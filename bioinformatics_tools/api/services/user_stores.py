"""
Each user's own copy of MARGIE's growing databases, on scratch, and their
backups on depot.

The databases a margie_sb run reads AND appends to -- the job database, the
OCC operon reference, the gene and operon fingerprint databases, the genome
pool -- used to be written in place on depot, or promoted to a per-user copy
beside the shared one, also on depot. Depot is shared and small; scratch is
not. So now:

  * The BASE copies on depot, under databases/margie-generated-databases/
    (margie-thesis-22-prokaryotes-base.db and the rest, below), are only ever
    read. Every new user starts from them.
  * On a user's first run, each base is copied to
        /scratch/<cluster>/<cluster-user>/margie-2026/<store>/<cluster-user>-<name>-v1<ext>
    and runs write there. (If the user already has a backup on depot -- say
    scratch was purged -- the newest backup is copied instead, as the next
    version, so no work is lost.)
  * "Back up" copies the working version vN to depot as
        databases/margie-generated-databases/<store>/<cluster-user>-<name>-vN<ext>
    beside the base -- only when the user says yes, after being told the size
    and how much room depot has left -- and the working copy on scratch
    becomes vN+1 (a rename:
    instant, nothing lost). Nothing is ever deleted: each backup is a new
    file, and earlier backups stay until someone removes them by hand.

The cluster username names everything -- it is what the files are owned by
on the cluster, and it does not change if someone makes a new web account.

The copies run as a SLURM job (charged to the account and partition in the
config), not on a login node: tens of GB of disk traffic is not login-node
work. Without an account they fall back to the login node, detached. Either
way it is one small bash script driven by a plan file, so the copy outlives
the request that started it. rsync reports its progress; GET /stores turns that into a
percentage for the progress bar and a log tail for the details under it.
"""

from __future__ import annotations

import json
import logging
import posixpath
import re
import shlex
import time

from bioinformatics_tools.utilities import ssh_sftp

LOGGER = logging.getLogger(__name__)

DEPOT = '/depot/lindems/data/margie'
# The base copies, one folder per database; each user's backups go beside them.
GENERATED = f'{DEPOT}/databases/margie-generated-databases'
# The folder under the user's scratch that holds every store (and the
# outputs that used to go to depot too).
ROOT_NAME = 'margie-2026'
# MARGIE's own bookkeeping inside that folder: what is where, and the copy under way.
STATE_DIR = '.margie'

# ---------------------------------------------------------------------------
# The stores. kind 'file': one file (plus companion files that travel with
# it). kind 'dir': a folder, copied whole.
#
# config: the config.yaml keys that point at the store -- for a 'dir', each
# maps to a file inside it ('' = the folder itself).
# ---------------------------------------------------------------------------
STORES: list[dict] = [
    {
        'id': 'job_db',
        'label': 'Job database',
        'note': 'Every run’s results and the job history, in SQLite.',
        'kind': 'file',
        'base': f'{GENERATED}/sqlite/margie-thesis-22-prokaryotes-base.db',
        'depot_sub': 'sqlite',
        'sub': 'sqlite',
        'name': 'margie-thesis-22-prokaryotes-base',
        'ext': '.db',
        'companions': [],
        'config': {'main_database': None},
    },
    {
        'id': 'occ',
        'label': 'Operon reference',
        'note': 'The cross-genome operon (OCC) reference scoring reads and adds each genome to.',
        'kind': 'file',
        'base': f'{GENERATED}/operon-database/occ_reference.pkl',
        'depot_sub': 'operon-database',
        'sub': 'operon-database',
        'name': 'occ_reference',
        'ext': '.pkl',
        # Written beside the reference by update-occ-reference-depot.py.
        'companions': ['.members.tsv', '.genome_stats.tsv'],
        'config': {'margie_sb.operon_database.occ_reference_pkl': None},
    },
    {
        'id': 'fingerprint',
        'label': 'Fingerprint databases',
        'note': 'The gene fingerprint database and the four operon fingerprint databases beside it.',
        'kind': 'dir',
        'base': f'{GENERATED}/fingerprint-database',
        # Only the databases themselves, not the backups and locks beside them.
        'base_files': [
            'fingerprint-database.tsv',
            'fingerprint-database-metadata.json',
            'operon-fingerprint-database-evidence-ordered.tsv',
            'operon-fingerprint-database-evidence-composition.tsv',
            'operon-fingerprint-database-label-ordered.tsv',
            'operon-fingerprint-database-label-composition.tsv',
            'operon-fingerprint-database-evidence_ordered-metadata.json',
            'operon-fingerprint-database-evidence_composition-metadata.json',
            'operon-fingerprint-database-label_ordered-metadata.json',
            'operon-fingerprint-database-label_composition-metadata.json',
        ],
        'depot_sub': 'fingerprint-database',
        'sub': 'fingerprint-database',
        'name': 'fingerprint-database',
        'config': {
            'margie_sb.fingerprint_database.path': 'fingerprint-database.tsv',
            # The report figures read the label-ordered operon fingerprints.
            'margie_sb.report_figures.operon_db': 'operon-fingerprint-database-label-ordered.tsv',
        },
    },
    {
        'id': 'genome_pool',
        'label': 'Genome pool',
        'note': 'Genomes ANI and AAI compare against; every annotated genome joins it.',
        'kind': 'dir',
        'base': f'{GENERATED}/genome-pool',
        # Its genomes only, not the backups kept beside them.
        'base_files': ['fna', 'faa'],
        'depot_sub': 'genome-pool',
        'sub': 'genome-pool',
        'name': 'genome-pool',
        'config': {'margie_sb.genome_pool.path': ''},
    },
]

# Written by runs but not databases: a folder each on scratch, never backed up.
OUTPUT_DIRS = {
    'margie_sb.scoring_results_historical.path': 'scoring-archive',
    'margie_sb.final_tables_depot.path': 'final-tables',
    'margie_sb.sqlite_pipeline_snapshot.path': 'sqlite/snapshots',
}

STORE_BY_ID = {s['id']: s for s in STORES}


def backup_root(cfg: dict | None) -> str:
    """Where backups go: margie_sb.backup_root if set (Settings, behind Unlock),
    else beside the bases on depot."""
    chosen = _cfg_get(cfg or {}, 'margie_sb.backup_root')
    return chosen.strip().rstrip('/') if isinstance(chosen, str) and chosen.strip() else GENERATED


def backup_dir(store: dict, cfg: dict | None) -> str:
    return f"{backup_root(cfg)}/{store['depot_sub']}"


class StoreError(Exception):
    """A request that cannot be done now; the message says why, for the page."""

    def __init__(self, message: str, status: int = 400):
        super().__init__(message)
        self.status = status


# ---------------------------------------------------------------------------
# Remote helpers
# ---------------------------------------------------------------------------

def _run(conn, command: str, timeout: float = 60.0) -> tuple[int, str]:
    ssh = conn.connect()
    _, stdout, stderr = ssh.exec_command(command, timeout=timeout)
    code = stdout.channel.recv_exit_status()
    out = stdout.read().decode(errors='replace') + stderr.read().decode(errors='replace')
    return code, out


def _cfg_get(cfg: dict, dotted: str):
    cur = cfg
    for part in dotted.split('.'):
        if not isinstance(cur, dict) or part not in cur:
            return None
        cur = cur[part]
    return cur


def _cfg_set(cfg: dict, dotted: str, value) -> None:
    parts = dotted.split('.')
    cur = cfg
    for part in parts[:-1]:
        if not isinstance(cur.get(part), dict):
            cur[part] = {}
        cur = cur[part]
    cur[parts[-1]] = value


def scratch_root(conn, cfg: dict) -> str:
    """margie_sb.stores_root if set, else /scratch/<cluster>/<user>/margie-2026.

    RCAC's scratch is /scratch/<cluster short name>/<user>; the short name is
    the login node's domain (login07.negishi.rcac.purdue.edu -> negishi).
    """
    chosen = _cfg_get(cfg, 'margie_sb.stores_root')
    if isinstance(chosen, str) and chosen.strip():
        return chosen.strip().rstrip('/')
    code, out = _run(conn, 'printf "%s %s" "$(hostname -d 2>/dev/null | cut -d. -f1)" "$USER"')
    cluster, _, user = out.strip().partition(' ')
    if code != 0 or not cluster or not user:
        raise StoreError('Could not work out your scratch folder on the cluster. Set margie_sb.stores_root in Settings.')
    return f'/scratch/{cluster}/{user}/{ROOT_NAME}'


# ---------------------------------------------------------------------------
# Names
# ---------------------------------------------------------------------------

def versioned(store: dict, user: str, version: int) -> str:
    """<user>-<name>-v<N><ext> -- the same form on scratch and on depot."""
    return f"{user}-{store['name']}-v{version}{store.get('ext', '') if store['kind'] == 'file' else ''}"


def _version_of(store: dict, user: str, name: str) -> int | None:
    ext = re.escape(store.get('ext', '')) if store['kind'] == 'file' else ''
    m = re.fullmatch(rf"{re.escape(user)}-{re.escape(store['name'])}-v(\d+){ext}", name)
    return int(m.group(1)) if m else None


def _backups(conn, store: dict, user: str, cfg: dict | None = None) -> list[int]:
    """Versions of this user's backups of the store on depot, oldest first."""
    try:
        entries = ssh_sftp.list_remote_dir(backup_dir(store, cfg), connection=conn)
    except Exception:
        return []
    want = 'file' if store['kind'] == 'file' else 'directory'
    found = []
    for e in entries:
        if e.get('type') != want:
            continue
        v = _version_of(store, user, e.get('name') or '')
        if v is not None:
            found.append(v)
    return sorted(found)


# ---------------------------------------------------------------------------
# State: <root>/.margie/stores.json holds each store's working version;
# <root>/.margie/op.* is the copy under way (or the last one).
# ---------------------------------------------------------------------------

def _state_path(root: str, name: str) -> str:
    return f'{root}/{STATE_DIR}/{name}'


def _read_text(conn, path: str) -> str:
    code, out = _run(conn, f'cat {shlex.quote(path)} 2>/dev/null')
    return out if code == 0 else ''


def read_manifest(conn, root: str) -> dict:
    try:
        data = json.loads(_read_text(conn, _state_path(root, 'stores.json')) or '{}')
        return data if isinstance(data, dict) else {}
    except json.JSONDecodeError:
        return {}


def _write_manifest(conn, root: str, manifest: dict) -> None:
    ssh_sftp.write_remote_text_file(_state_path(root, 'stores.json'), json.dumps(manifest, indent=2) + '\n', connection=conn)


def working_path(store: dict, root: str, user: str, version: int) -> str:
    return f"{root}/{store['sub']}/{versioned(store, user, version)}"


def _probe(conn, root: str, user: str, cfg: dict | None = None) -> dict:
    """Everything the page needs, in one trip to the cluster: which working
    copies exist, the backups on depot, and the copy under way (its status,
    the tail of its log, and whether its script is still alive)."""
    manifest_path = _state_path(root, 'stores.json')
    status_path = _state_path(root, 'op.status')
    log_path = _state_path(root, 'op.log')
    parts = [f'echo "@@manifest"; cat {shlex.quote(manifest_path)} 2>/dev/null',
             f'echo "@@status"; cat {shlex.quote(status_path)} 2>/dev/null',
             f'echo "@@log"; tail -c 16000 {shlex.quote(log_path)} 2>/dev/null | tr "\\r" "\\n"',
             'echo "@@host"; hostname',
             f'echo "@@squeue"; j=$(cat {shlex.quote(_state_path(root, "op.jobid"))} 2>/dev/null); '
             f'[ -n "$j" ] && {{ echo "job=$j"; squeue -h -j "$j" -o %T 2>/dev/null; }}',
             f'echo "@@slurmout"; tail -n 20 {shlex.quote(_state_path(root, "op.slurm.out"))} 2>/dev/null',
             f'echo "@@logage"; echo $(( $(date +%s) - $(stat -c %Y {shlex.quote(log_path)} 2>/dev/null || echo 0) ))']
    for st in STORES:
        parts.append(f'echo "@@depot {st["id"]}"; ls -1p {shlex.quote(backup_dir(st, cfg))} 2>/dev/null | grep "^{user}-{st["name"]}-v"')
    code, out = _run(conn, '; '.join(parts), timeout=60)
    sections: dict[str, list[str]] = {}
    current = None
    for line in out.splitlines():
        if line.startswith('@@'):
            current = line[2:].strip()
            sections[current] = []
        elif current is not None:
            sections[current].append(line)
    try:
        manifest = json.loads('\n'.join(sections.get('manifest', [])) or '{}')
    except json.JSONDecodeError:
        manifest = {}
    # The working copies, now that we know their versions (a second, small trip).
    wanted = {}
    for st in STORES:
        v = (manifest.get(st['id']) or {}).get('version') if isinstance(manifest, dict) else None
        if v:
            wanted[st['id']] = (int(v), working_path(st, root, user, int(v)))
    present = set()
    if wanted:
        cmd = '; '.join(f'test -e {shlex.quote(p)} && echo {sid}' for sid, (_, p) in wanted.items())
        _, out2 = _run(conn, cmd)
        present = set(out2.split())
    backups = {}
    for st in STORES:
        names = [n.rstrip('/') for n in sections.get(f'depot {st["id"]}', [])]
        is_dir = st['kind'] == 'dir'
        vs = sorted(v for n in names
                    if (v := _version_of(st, user, n)) is not None
                    and any(x == n + '/' for x in sections.get(f'depot {st["id"]}', [])) == is_dir)
        backups[st['id']] = vs
    return {
        'manifest': manifest if isinstance(manifest, dict) else {},
        'wanted': wanted,
        'present': present,
        'backups': backups,
        'op': _op_from(sections, root),
    }


def _op_from(sections: dict, root: str) -> dict:
    """The copy under way or last done: its status file, and progress read from its log."""
    op: dict = {}
    for line in sections.get('status', []):
        k, _, v = line.partition('=')
        if k:
            op[k] = v
    if not op:
        return {}
    # A 'running' op whose script is gone (the node rebooted, it was killed)
    # is failed, not running forever. Its pid means something only on the
    # login node it runs on, and this connection may have landed on another:
    # there, a log that has stopped growing for five minutes is the sign.
    squeue = [l.strip() for l in sections.get('squeue', []) if l.strip()]
    job = next((l[4:] for l in squeue if l.startswith('job=')), '')
    slurm_state = next((l for l in squeue if not l.startswith('job=')), '')
    if job:
        op['job'] = job
        op['slurm_state'] = slurm_state
    if op.get('state') in ('queued', 'running') and job:
        # A SLURM job: squeue knows whether it is still there.
        if not slurm_state:
            op['state'] = 'failed'
            op['message'] = op.get('message') or 'The copy job ended before it finished.'
            op['slurm_out'] = sections.get('slurmout', [])[-10:]
    elif op.get('state') == 'running':
        here = (sections.get('host') or [''])[0].strip()
        age = (sections.get('logage') or ['0'])[0].strip()
        if here and here == op.get('host') and op.get('pid', '').isdigit():
            op['_check_pid'] = op['pid']
        elif age.isdigit() and int(age) > 300:
            op['state'] = 'failed'
            op['message'] = op.get('message') or 'The copy stopped before it finished.'
    lines = [l.rstrip() for l in sections.get('log', []) if l.strip()]
    # rsync --info=progress2 prints "   1,234,567  42%  ..." for the file in hand.
    pct = 0.0
    for l in reversed(lines):
        m = re.search(r'\s(\d{1,3})%\s', f' {l} ')
        if m:
            pct = float(m.group(1))
            break
    total = float(op.get('total_bytes') or 0)
    done = float(op.get('done_bytes') or 0)
    item = float(op.get('item_bytes') or 0)
    if op.get('state') == 'done':
        percent = 100.0
    elif total > 0:
        percent = min(100.0, 100.0 * (done + item * pct / 100.0) / total)
    else:
        percent = pct
    op['percent'] = round(percent, 1)
    # rsync's own progress lines are the bar's business, not the log's.
    op['log'] = [l for l in lines if not re.match(r'^\s*[\d,]+\s+\d{1,3}%', l)][-40:]
    return op


def _read_op(conn, root: str) -> dict:
    op = _probe(conn, root, '__none__')['op']
    return _settle_pid(conn, op)


def _settle_pid(conn, op: dict) -> dict:
    """On the login node the copy runs on, its pid says whether it is alive."""
    pid = op.pop('_check_pid', None)
    if pid:
        code, _ = _run(conn, f'kill -0 {int(pid)} 2>/dev/null')
        if code != 0:
            op['state'] = 'failed'
            op['message'] = op.get('message') or 'The copy stopped before it finished.'
    return op


# ---------------------------------------------------------------------------
# The copy script: reads a plan, one step per line, and reports as it goes.
#   copy<TAB>label<TAB>src<TAB>dst<TAB>bytes   rsync src -> dst (via dst.partial)
#   files<TAB>label<TAB>srcdir<TAB>dst<TAB>bytes<TAB>name,name,...
#   move<TAB>src<TAB>dst                        rename (same filesystem: instant)
#   mkdir<TAB>path
#   link<TAB>src<TAB>dst                        symlink dst -> src, unless dst exists
#   pull<TAB>label<TAB>url<TAB>dst              apptainer pull url -> dst (via dst.partial)
#   repo<TAB>url                                where margie-build comes from, for build
#   build<TAB>label<TAB>mode<TAB>tool<TAB>sifdir<TAB>dbdir<TAB>statement
#                                               margie-build's build.sh --<mode> <tool>;
#                                               a statement accepts a gated tool's licence
# (link, pull, repo and build set up the tools: api/services/tool_assets.py.)
# ---------------------------------------------------------------------------
COPY_SCRIPT = r'''#!/bin/bash
set -u
dir="$1"; plan="$dir/op.plan"; status="$dir/op.status"; log="$dir/op.log"
total=$(awk -F'\t' '$1=="copy"||$1=="files"{s+=$5} END{print s+0}' "$plan")
done_bytes=0
say() { { printf 'state=%s\npid=%s\nhost=%s\nop=%s\nlabel=%s\ntotal_bytes=%s\ndone_bytes=%s\nitem_bytes=%s\nmessage=%s\n' \
          "$1" "$$" "$(hostname)" "$OP" "$2" "$total" "$done_bytes" "$3" "$4"; } > "$status.tmp" && mv "$status.tmp" "$status"; }
fail() { say failed "$1" 0 "$2"; echo "FAILED: $2" >> "$log"; exit 1; }
OP="$(head -n 1 "$dir/op.kind" 2>/dev/null)"
: > "$log"
say running "Starting" 0 ""
while IFS=$'\t' read -r kind a b c d e f; do
  case "$kind" in
    copy)
      label="$a"; src="$b"; dst="$c"; bytes="$d"
      say running "$label" "$bytes" ""
      echo "== $label: $src -> $dst" >> "$log"
      mkdir -p "$(dirname "$dst")" || fail "$label" "could not make $(dirname "$dst")"
      if [ -d "$src" ]; then
        rsync -a --info=progress2 --no-inc-recursive "$src/" "$dst.partial/" >> "$log" 2>&1 || fail "$label" "copy of $src failed"
      else
        rsync -a --info=progress2 "$src" "$dst.partial" >> "$log" 2>&1 || fail "$label" "copy of $src failed"
      fi
      mv -T "$dst.partial" "$dst" || fail "$label" "could not put $dst in place"
      done_bytes=$(( done_bytes + bytes ))
      ;;
    files)
      label="$a"; src="$b"; dst="$c"; bytes="$d"; names="$e"
      say running "$label" "$bytes" ""
      echo "== $label: $src/{$names} -> $dst" >> "$log"
      mkdir -p "$dst.partial" || fail "$label" "could not make $dst"
      list=(); IFS=',' read -r -a want <<< "$names"
      for n in "${want[@]}"; do [ -e "$src/$n" ] && list+=("$src/$n"); done
      if [ ${#list[@]} -gt 0 ]; then
        rsync -a --info=progress2 "${list[@]}" "$dst.partial/" >> "$log" 2>&1 || fail "$label" "copy from $src failed"
      fi
      mv -T "$dst.partial" "$dst" || fail "$label" "could not put $dst in place"
      done_bytes=$(( done_bytes + bytes ))
      ;;
    move)
      echo "== rename $a -> $b" >> "$log"
      mv -T "$a" "$b" || fail "Renaming" "could not rename $a"
      ;;
    mkdir)
      mkdir -p "$a" || fail "Folders" "could not make $a"
      ;;
    link)
      [ -e "$b" ] || [ -L "$b" ] || ln -s "$a" "$b" || fail "Links" "could not link $b"
      ;;
    pull)
      label="$a"; url="$b"; dst="$c"
      say running "$label" 0 ""
      echo "== $label: $url -> $dst" >> "$log"
      app=$(command -v apptainer || command -v singularity) || fail "$label" "apptainer is not installed here"
      mkdir -p "$(dirname "$dst")" || fail "$label" "could not make $(dirname "$dst")"
      # Its cache and scratch space beside this script, not in the small home folder.
      export APPTAINER_CACHEDIR="$dir/apptainer-cache" APPTAINER_TMPDIR="$dir/apptainer-tmp"
      mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
      "$app" pull --force "$dst.partial" "$url" >> "$log" 2>&1 || fail "$label" "could not pull $url"
      mv -T "$dst.partial" "$dst" || fail "$label" "could not put $dst in place"
      ;;
    repo)
      repo_url="$a"
      ;;
    build)
      label="$a"; mode="$b"; tool="$c"; sifdir="$d"; dbdir="$e"; statement="$f"
      say running "$label" 0 ""
      echo "== $label" >> "$log"
      repo="$dir/margie-build"
      if [ ! -d "$repo/.git" ]; then
        GIT_TERMINAL_PROMPT=0 git clone -q --depth 1 "${repo_url:-}" "$repo" >> "$log" 2>&1 \
          || fail "$label" "could not download margie-build from ${repo_url:-nowhere} (is it public?)"
      fi
      accept=()
      [ -n "$statement" ] && accept=(--accept-"$tool"-licence --licence-statement "$statement")
      mkdir -p "$sifdir" "$dbdir" "$dir/apptainer-cache" "$dir/apptainer-tmp" || fail "$label" "could not make $sifdir or $dbdir"
      ( cd "$repo" && APPTAINER_CACHEDIR="$dir/apptainer-cache" APPTAINER_TMPDIR="$dir/apptainer-tmp" \
          THREADS="${SLURM_CPUS_PER_TASK:-4}" ./build.sh --"$mode" "$tool" --sif-dir "$sifdir" --db-dir "$dbdir" \
          --audit-log "$dir/licence-acceptances.tsv" --no-color --quiet-notice ${accept[@]+"${accept[@]}"} ) >> "$log" 2>&1 \
        || fail "$label" "margie-build could not build $tool (the log above says why)"
      ;;
  esac
done < "$plan"
say done "Finished" 0 ""
echo "== finished" >> "$log"
'''


def _copy_cpus(cfg: dict | None) -> int:
    """Cores for the copy job: margie_sb.stores_copy_cpus, 4 by default, never
    more than 16 (the lab's limit for copying and moving). A copy is limited
    by the disks, not the processor, so more cores would mostly sit idle."""
    try:
        n = int(_cfg_get(cfg or {}, 'margie_sb.stores_copy_cpus') or 4)
    except (TypeError, ValueError):
        n = 4
    return max(1, min(16, n))


def _start(conn, root: str, kind: str, plan: list[list], after: dict, cfg: dict | None = None,
           walltime: str = '12:00:00', mem: str = '4G') -> None:
    """Write the plan and start the copy script: a SLURM job when the config
    names an account, detached on the login node otherwise."""
    state = f'{root}/{STATE_DIR}'
    code, out = _run(conn, f'mkdir -p {shlex.quote(state)} && chmod 700 {shlex.quote(state)}')
    if code != 0:
        raise StoreError(f'Could not make {state}: {out.strip()}')
    ssh_sftp.write_remote_text_file(f'{state}/op.sh', COPY_SCRIPT, connection=conn)
    ssh_sftp.write_remote_text_file(
        f'{state}/op.plan', ''.join('\t'.join(str(x) for x in step) + '\n' for step in plan), connection=conn)
    ssh_sftp.write_remote_text_file(f'{state}/op.kind', kind + '\n', connection=conn)
    # What to record once it has finished (see apply_finished).
    ssh_sftp.write_remote_text_file(f'{state}/op.after.json', json.dumps(after) + '\n', connection=conn)
    _run(conn, f'cd {shlex.quote(state)} && rm -f op.status op.log op.jobid op.slurm.out')
    account = str(_cfg_get(cfg or {}, 'compute.cluster_default.account') or '').strip()
    partition = str(_cfg_get(cfg or {}, 'compute.cluster_default.partition') or '').strip()
    q = shlex.quote
    if account:
        # Until the job starts, the page shows it waiting in the queue.
        ssh_sftp.write_remote_text_file(
            f'{state}/op.status', f'state=queued\nop={kind}\nlabel=Waiting for a SLURM slot\n', connection=conn)
        cmd = (f'cd {q(state)} && sbatch --parsable --job-name=margie-databases --account={q(account)} '
               + (f'--partition={q(partition)} ' if partition else '')
               + f'--time={q(walltime)} --ntasks=1 --cpus-per-task={_copy_cpus(cfg)} --mem={q(mem)} --output={q(state)}/op.slurm.out '
               + f'op.sh {q(state)}')
        code, out = _run(conn, cmd)
        job = out.strip().split(';')[0].split()[-1] if out.strip() else ''
        if code != 0 or not job.isdigit():
            raise StoreError(f'SLURM would not take the copy job: {out.strip()}')
        ssh_sftp.write_remote_text_file(f'{state}/op.jobid', job + '\n', connection=conn)
        return
    code, out = _run(conn, f'cd {q(state)} && setsid nohup bash op.sh {q(state)} '
                            f'> /dev/null 2>&1 < /dev/null & echo started')
    if 'started' not in out:
        raise StoreError(f'Could not start the copy: {out.strip()}')
    # The status file appears within a moment; wait for it so the page never
    # sees an op that is neither running nor finished.
    for _ in range(20):
        if _read_text(conn, f'{state}/op.status'):
            break
        time.sleep(0.25)


def _size(conn, paths: list[str]) -> int:
    if not paths:
        return 0
    code, out = _run(conn, 'du -sbc ' + ' '.join(shlex.quote(p) for p in paths) + ' 2>/dev/null | tail -1 | cut -f1',
                     timeout=300)
    try:
        return int(out.strip() or 0)
    except ValueError:
        return 0


# ---------------------------------------------------------------------------
# What the page asks
# ---------------------------------------------------------------------------

def _config_paths(store: dict, target: str) -> dict:
    """The config values for a store whose working copy is at target."""
    return {key: (posixpath.join(target, inner) if inner else target) if store['kind'] == 'dir' else target
            for key, inner in store['config'].items()}


def status(conn, cfg: dict, user: str) -> dict:
    root = scratch_root(conn, cfg)
    probe = _probe(conn, root, user, cfg)
    op = _settle_pid(conn, probe['op'])
    stores = []
    ready = True
    for s in STORES:
        v, path = probe['wanted'].get(s['id'], (None, ''))
        present = s['id'] in probe['present']
        ready = ready and present
        backups = probe['backups'].get(s['id'], [])
        stores.append({
            'id': s['id'],
            'label': s['label'],
            'note': s['note'],
            'version': v if present else None,
            'path': path if present else '',
            'base': s['base'],
            'backup': ({'version': backups[-1],
                        'path': posixpath.join(backup_dir(s, cfg), versioned(s, user, backups[-1]))}
                       if backups else None),
            'backups': len(backups),
        })
    return {'root': root, 'user': user, 'ready': ready, 'stores': stores, 'op': op or None}


def progress(conn, cfg: dict) -> dict | None:
    """Only the copy under way: what the progress bar polls, every second or two."""
    root = scratch_root(conn, cfg)
    return _read_op(conn, root) or None


def _busy(op: dict | None) -> bool:
    return bool(op) and op.get('state') in ('queued', 'running')


def start_setup(conn, cfg: dict, user: str) -> dict:
    """Copy whichever stores are missing on scratch: the newest backup if the
    user has one, otherwise the base."""
    st = status(conn, cfg, user)
    if _busy(st['op']):
        raise StoreError('A copy is already under way.', 409)
    root = st['root']
    plan: list[list] = [['mkdir', f'{root}/{sub}'] for sub in OUTPUT_DIRS.values()]
    after: dict = {'kind': 'setup', 'versions': {}}
    for s, info in zip(STORES, st['stores']):
        if info['version']:
            continue
        backups = _backups(conn, s, user, cfg)
        if backups:
            src = posixpath.join(backup_dir(s, cfg), versioned(s, user, backups[-1]))
            version = backups[-1] + 1
            label = f"{s['label']} (from your backup v{backups[-1]})"
        else:
            src = s['base']
            version = 1
            label = s['label']
        dst = working_path(s, root, user, version)
        if s['kind'] == 'file':
            plan.append(['copy', label, src, dst, _size(conn, [src])])
            for c in s['companions']:
                code, _ = _run(conn, f'test -e {shlex.quote(src + c)}')
                if code == 0:
                    plan.append(['copy', f"{s['label']} ({c.lstrip('.')})", src + c, dst + c, _size(conn, [src + c])])
        elif s.get('base_files') and src == s['base']:
            files = [posixpath.join(src, f) for f in s['base_files']]
            plan.append(['files', label, src, dst, _size(conn, files), ','.join(s['base_files'])])
        else:
            plan.append(['copy', label, src, dst, _size(conn, [src])])
        after['versions'][s['id']] = version
    _start(conn, root, 'setup', plan, after, cfg)
    return status(conn, cfg, user)


_UNITS = {'B': 1, 'KB': 1024, 'MB': 1024 ** 2, 'GB': 1024 ** 3, 'TB': 1024 ** 4, 'PB': 1024 ** 5,
          'K': 1024, 'M': 1024 ** 2, 'G': 1024 ** 3, 'T': 1024 ** 4, 'P': 1024 ** 5}


def _bytes(text: str) -> int | None:
    m = re.fullmatch(r'([\d.]+)\s*([KMGTP]?B?)', text.strip(), re.I)
    return int(float(m.group(1)) * _UNITS.get(m.group(2).upper() or 'B', 1)) if m else None


def depot_free(conn, where: str = DEPOT) -> int | None:
    """Bytes a backup can still use on depot: the smaller of what the group's
    quota leaves (RCAC's myquota) and what the filesystem itself has free."""
    group = where.split('/')[2] if where.startswith('/depot/') and len(where.split('/')) > 2 else ''
    # df on the nearest folder that exists: the backup folder may not yet.
    code, out = _run(conn, f'myquota 2>/dev/null; echo "@@df"; d={shlex.quote(where)}; '
                           f'while [ ! -e "$d" ] && [ "$d" != / ]; do d=$(dirname "$d"); done; '
                           f'df -B1 --output=avail "$d" 2>/dev/null | tail -1')
    quota_free = None
    head, _, df = out.partition('@@df')
    for line in head.splitlines():
        cols = line.split()
        if group and len(cols) >= 4 and cols[0] == 'depot' and cols[1] == group:
            used, limit = _bytes(cols[2]), _bytes(cols[3])
            if used is not None and limit is not None:
                quota_free = max(0, limit - used)
    try:
        fs_free = int(df.strip())
    except ValueError:
        fs_free = None
    known = [x for x in (quota_free, fs_free) if x is not None]
    return min(known) if known else None


def backup_check(conn, cfg: dict, user: str, store_id: str) -> dict:
    """What a backup would take and whether depot has room -- asked before the
    Yes / No, and again when Yes is pressed."""
    s = STORE_BY_ID.get(store_id)
    if not s:
        raise StoreError('No such database.', 404)
    st = status(conn, cfg, user)
    info = next(i for i in st['stores'] if i['id'] == store_id)
    if not info['version']:
        raise StoreError(f"{s['label']} is not set up on scratch yet.", 409)
    src = info['path']
    paths = [src] + [src + c for c in (s.get('companions') or [])]
    size = _size(conn, paths)
    free = depot_free(conn, backup_root(cfg))
    # Room to spare: a depot filled to the last byte breaks everyone's work.
    margin = 5 * 1024 ** 3
    fits = free is None or free - margin >= size
    return {
        'id': store_id,
        'label': s['label'],
        'version': info['version'],
        'size': size,
        'free': free,
        'fits': fits,
        'target': posixpath.join(backup_dir(s, cfg), versioned(s, user, info['version'])),
    }


def start_backup(conn, cfg: dict, user: str, store_id: str) -> dict:
    """Copy the working vN to depot as <user>-<name>-vN; the working copy becomes vN+1."""
    s = STORE_BY_ID.get(store_id)
    if not s:
        raise StoreError('No such database.', 404)
    st = status(conn, cfg, user)
    if _busy(st['op']):
        raise StoreError('A copy is already under way.', 409)
    check = backup_check(conn, cfg, user, store_id)
    if not check['fits']:
        gb = lambda n: f'{n / 1024 ** 3:.1f} GB'
        raise StoreError(f"Depot has {gb(check['free'])} left and this backup needs {gb(check['size'])} "
                         '(with 5 GB kept free): it would not fit.', 409)
    info = next(i for i in st['stores'] if i['id'] == store_id)
    v = info['version']
    if not v:
        raise StoreError(f"{s['label']} is not set up on scratch yet.", 409)
    root = st['root']
    src = working_path(s, root, user, v)
    dst = posixpath.join(backup_dir(s, cfg), versioned(s, user, v))
    code, _ = _run(conn, f'test -e {shlex.quote(dst)}')
    if code == 0:
        raise StoreError(f'{dst} already exists; nothing is overwritten.', 409)
    nxt = working_path(s, root, user, v + 1)
    plan: list[list] = [['copy', f"{s['label']} v{v} to depot", src, dst, _size(conn, [src])]]
    companions = s.get('companions', []) if s['kind'] == 'file' else []
    for c in companions:
        code, _ = _run(conn, f'test -e {shlex.quote(src + c)}')
        if code == 0:
            plan.append(['copy', f"{s['label']} ({c.lstrip('.')})", src + c, dst + c, _size(conn, [src + c])])
    # Only once the backup is in place: the working copy moves on to the next version.
    plan.append(['move', src, nxt])
    for c in companions:
        code, _ = _run(conn, f'test -e {shlex.quote(src + c)}')
        if code == 0:
            plan.append(['move', src + c, nxt + c])
    _start(conn, root, f'backup:{store_id}', plan, {'kind': 'backup', 'versions': {store_id: v + 1}}, cfg)
    return status(conn, cfg, user)


def apply_finished(conn, cfg: dict, user: str) -> bool:
    """After a copy has finished: record the new working versions, and point
    the config at them. Returns True when the config changed (the caller
    saves it). Safe to call any number of times."""
    root = scratch_root(conn, cfg)
    op = _read_op(conn, root)
    if op.get('state') != 'done':
        return False
    state = f'{root}/{STATE_DIR}'
    after_text = _read_text(conn, f'{state}/op.after.json')
    changed = False
    if after_text:
        try:
            after = json.loads(after_text)
        except json.JSONDecodeError:
            after = {}
        manifest = read_manifest(conn, root)
        for sid, version in (after.get('versions') or {}).items():
            manifest[sid] = {'version': int(version)}
        _write_manifest(conn, root, manifest)
        # Settings the copy moved (the containers' and reference databases' folders).
        for key, value in (after.get('config') or {}).items():
            if _cfg_get(cfg, key) != value:
                _cfg_set(cfg, key, value)
                changed = True
        _run(conn, f'mv -f {shlex.quote(state)}/op.after.json {shlex.quote(state)}/op.applied.json')
    if op.get('op') == 'assets':
        # Not the databases' copy: their settings stay as they are.
        return changed
    return point_config(conn, cfg, user, root) or changed


def point_config(conn, cfg: dict, user: str, root: str | None = None) -> bool:
    """Set every store's and output folder's config key to its scratch path."""
    root = root or scratch_root(conn, cfg)
    manifest = read_manifest(conn, root)
    changed = False
    wanted: dict = {}
    for s in STORES:
        v = (manifest.get(s['id']) or {}).get('version')
        if v:
            wanted.update(_config_paths(s, working_path(s, root, user, v)))
    for key, sub in OUTPUT_DIRS.items():
        wanted[key] = f'{root}/{sub}'
    # Remembered, so later requests need not ask the cluster where scratch is.
    wanted['margie_sb.stores_root'] = root
    for key, value in wanted.items():
        if _cfg_get(cfg, key) != value:
            _cfg_set(cfg, key, value)
            changed = True
    return changed


def is_ready(conn, cfg: dict, user: str) -> bool:
    return status(conn, cfg, user)['ready']
