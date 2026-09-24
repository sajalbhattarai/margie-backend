"""
The tools' containers and reference databases a margie_sb run reads: whether
each is where the account's config points, whether a run can start without
the missing ones, and setting up those that are not there.

Where a run looks (workflow_helpers.sif_path / db_path):
  * containers: <margie_sb.sif_path>/<tool>.sif
  * databases:  db.<tool> when the config names one, else <margie_sb.db_root>/<tool>
Both default to the lab's shared folders on depot (WORKFLOW_PATH_DEFAULTS).

A missing one comes from the first of these that has it:
  1. copy  -- the lab's shared folder, when the config points somewhere else
              and the shared copy is there (and readable from this account);
  2. pull  -- margie_sb.container_registry (e.g. docker://ghcr.io/<owner>),
              containers only, tag margie_sb.container_tag (default latest);
  3. build -- margie-build's recipes (build.sh), run on the cluster: phases
              1-8's containers, and the 15 reference databases. Its five
              licence-gated tools build only with the licence accepted, in
              words, for that run.
Nothing already there is replaced or deleted, and nothing is written under
the lab's depot folder. When the configured folder is there, or cannot be
written to, what is set up goes under the user's scratch root instead and
the config is pointed there: the containers folder as a whole (the ones
already there are linked into it), a database one db.<tool> at a time.

The work runs as the same SLURM copy job as the user's databases
(user_stores), with its progress bar: one copy at a time.
"""

from __future__ import annotations

import posixpath
import shlex

from bioinformatics_tools.api.services import user_stores
from bioinformatics_tools.api.services.user_stores import StoreError, _cfg_get, _run
from bioinformatics_tools.workflow_tools.workflow_helpers import WORKFLOW_PATH_DEFAULTS, genome_calls
from bioinformatics_tools.workflow_tools.workflow_registry import MARGIE_SB_PHASED_TOOLS

WORKFLOW = 'margie_sb'
SHARED_SIF = WORKFLOW_PATH_DEFAULTS[WORKFLOW]['sif_path']
SHARED_DB = WORKFLOW_PATH_DEFAULTS[WORKFLOW]['db_root']

# Folders under the user's scratch root for what cannot go where the config points.
SIF_SUB = 'sif'
DB_SUB = 'reference-databases'

# Run only when chosen, and heavy: never what stops a run that did not choose them.
OPTIONAL = {'gtdbtk', 'llm'}

# Not a phase of its own: run_prodigal calls the genes of a genome without a
# domain or genetic code while GTDB-Tk is off (workflow_helpers.genome_calls).
GENE_CALLER = {'key': 'prodigal', 'label': 'Prodigal', 'sif': 'prodigal.sif'}

# What margie-build (build.sh) can make.
BUILD_REPO = 'https://github.com/sajalbhattarai/margie-build.git'
BUILDABLE_IMAGES = {'quast', 'gtdbtk', 'prodigal', 'rasttk', 'cog', 'dbcan', 'eggnog', 'geneprop', 'interpro', 'kegg',
                    'merops', 'pfam', 'pgap', 'tcdb', 'tigrfam', 'uniprot', 'operon', 'phobius', 'tmbed', 'envelope',
                    'deepsig', 'psortb'}
BUILDABLE_DATABASES = {'cog', 'dbcan', 'eggnog', 'geneprop', 'gtdbtk', 'interpro', 'kegg', 'merops', 'pfam', 'pgap',
                       'rasttk', 'tcdb', 'tigrfam', 'tmbed', 'uniprot'}
# build.sh refuses these without --accept-<tool>-licence and this statement, word for word.
GATED = {'interpro', 'merops', 'tcdb', 'tmbed', 'phobius'}
LICENCE_STATEMENT = ('I accept that I am using these tools for non-commercial purposes '
                     'and have received all permissions from the upstream developers.')

# Room kept free on scratch, as a backup keeps on depot.
MARGIN = 5 * 1024 ** 3


def _expand(path: str, home: str) -> str:
    """~ is the account's home on the cluster, not this server's."""
    path = str(path).strip()
    if path == '~' or path.startswith('~/'):
        path = home.rstrip('/') + path[1:]
    return path.rstrip('/') or '/'


def _in_lab(path: str) -> bool:
    """Under the lab's depot folder: only ever read, even by an account that
    could write there (depot is shared, and nearly full)."""
    return path == user_stores.DEPOT or path.startswith(user_stores.DEPOT + '/')


def folders(cfg: dict, home: str) -> tuple[str, str]:
    """The containers folder and the databases root the run will use."""
    sif_dir = _cfg_get(cfg, f'{WORKFLOW}.sif_path') or SHARED_SIF
    db_root = _cfg_get(cfg, f'{WORKFLOW}.db_root') or SHARED_DB
    return _expand(sif_dir, home), _expand(db_root, home)


def inventory(cfg: dict, home: str) -> list[dict]:
    """Every container and reference database margie_sb can use, and where the run looks for it."""
    sif_dir, db_root = folders(cfg, home)
    rows: list[dict] = []

    def image(key: str, label: str, sif: str) -> None:
        rows.append({'id': f'image:{key}', 'kind': 'image', 'tool': key, 'label': label, 'name': sif,
                     'path': f'{sif_dir}/{sif}', 'shared': f'{SHARED_SIF}/{sif}', 'optional': key in OPTIONAL})

    for t in MARGIE_SB_PHASED_TOOLS:
        if t.get('uses_container', True):
            image(t['key'], t['label'], t['sif'])
        if t['key'] == 'rasttk':
            image(GENE_CALLER['key'], GENE_CALLER['label'], GENE_CALLER['sif'])
    for t in MARGIE_SB_PHASED_TOOLS:
        if not t.get('database'):
            continue
        key = t['key']
        explicit = _cfg_get(cfg, f'db.{key}')
        rows.append({'id': f'database:{key}', 'kind': 'database', 'tool': key, 'label': t['label'], 'name': key,
                     'path': _expand(explicit, home) if explicit else f'{db_root}/{key}',
                     'shared': f'{SHARED_DB}/{key}', 'optional': key in OPTIONAL})
    return rows


# One trip to the cluster: for each row, whether it is there (ok / missing /
# unknown, when the folder cannot be read from this account), whether the
# shared copy is, and whether the folder it belongs in can be written to.
_PROBE = r'''
near() { d="$1"; while [ ! -e "$d" ] && [ "$d" != / ]; do d=$(dirname "$d"); done; printf '%s' "$d"; }
has() {
  if [ -d "$1" ]; then
    if [ -r "$1" ] && [ -x "$1" ]; then [ -n "$(ls -A "$1" 2>/dev/null | head -c 1)" ] && echo ok || echo missing; else echo unknown; fi
  elif [ -e "$1" ]; then [ -s "$1" ] && echo ok || echo missing
  else d=$(near "$(dirname "$1")"); [ -r "$d" ] && [ -x "$d" ] && echo missing || echo unknown; fi
}
writable() { d=$(near "$1"); [ -w "$d" ] && [ -x "$d" ] && echo w || echo r; }
'''


def status(conn, cfg: dict, home: str) -> dict:
    rows = inventory(cfg, home)
    sif_dir, db_root = folders(cfg, home)
    q = shlex.quote
    lines = [_PROBE, f'echo "sif $(writable {q(sif_dir)})"', f'echo "db $(writable {q(db_root)})"']
    for i, r in enumerate(rows):
        shared = f'$(has {q(r["shared"])})' if r['shared'] != r['path'] else 'same'
        lines.append(f'echo "{i} $(has {q(r["path"])}) {shared} $(writable {q(posixpath.dirname(r["path"]))})"')
    code, out = _run(conn, '\n'.join(lines), timeout=120)
    if code != 0 and not out.strip():
        raise StoreError('Could not look at the containers and databases folders on the cluster.', 502)
    seen: dict[str, list[str]] = {}
    for line in out.splitlines():
        parts = line.split()
        if parts:
            seen[parts[0]] = parts[1:]
    registry = str(_cfg_get(cfg, f'{WORKFLOW}.container_registry') or '').strip().rstrip('/')
    tag = str(_cfg_get(cfg, f'{WORKFLOW}.container_tag') or 'latest').strip()
    for i, r in enumerate(rows):
        got = seen.get(str(i), ['unknown', 'unknown', 'r'])
        r['status'] = got[0] if got[0] in ('ok', 'missing', 'unknown') else 'unknown'
        shared_ok = len(got) > 1 and got[1] == 'ok'
        r['writable'] = len(got) > 2 and got[2] == 'w' and not _in_lab(r['path'])
        stem = r['name'][:-4] if r['name'].endswith('.sif') else r['name']
        r['source'] = None
        r['url'] = ''
        if r['status'] != 'ok':
            if shared_ok:
                r['source'] = 'copy'
            elif r['kind'] == 'image' and registry:
                r['source'] = 'pull'
                r['url'] = f'{registry}/{stem.lower()}:{tag}'
            elif r['tool'] in (BUILDABLE_IMAGES if r['kind'] == 'image' else BUILDABLE_DATABASES):
                r['source'] = 'build'
        r['gated'] = r['source'] == 'build' and r['tool'] in GATED
    return {
        'sif_dir': sif_dir,
        'db_root': db_root,
        'sif_dir_writable': seen.get('sif', ['r'])[0] == 'w' and not _in_lab(sif_dir),
        'db_root_writable': seen.get('db', ['r'])[0] == 'w' and not _in_lab(db_root),
        'rows': rows,
    }


# ---------------------------------------------------------------------------
# Before a run
# ---------------------------------------------------------------------------

def needed(rows: list[dict], cfg: dict, tools: set[str], genome_names: list[str] | None) -> set[str]:
    """The ids of the rows a run of these tools reads.

    RASTtk is always part of a run (phase 3's gate). Which gene caller a genome
    gets is decided from the config (genome_calls): with GTDB-Tk on, RASTtk
    calls them all; with it off, a genome without a domain or genetic code goes
    to Prodigal. Without the genomes' names, only RASTtk is asked for.
    """
    tools = set(tools) | {'rasttk'}
    if 'gtdbtk' in tools:
        callers = {'rasttk'}
    elif genome_names is None:
        callers = {'rasttk'}
    else:
        calls = genome_calls({n: n for n in genome_names}, {**cfg, 'run_gtdbtk': False})
        callers = {c['gene_caller'] for c in calls.values()}
    out = set()
    for r in rows:
        if r['tool'] in ('rasttk', 'prodigal'):
            if r['tool'] in callers:
                out.add(r['id'])
        elif r['tool'] in tools:
            out.add(r['id'])
    return out


def missing_for_run(conn, cfg: dict, home: str, tools: set[str], genome_names: list[str] | None) -> list[dict]:
    """What this run needs and is not there. "Cannot tell" (a folder this
    account cannot read) is not missing: the run is the judge of that."""
    st = status(conn, cfg, home)
    want = needed(st['rows'], cfg, tools, genome_names)
    return [r for r in st['rows'] if r['id'] in want and r['status'] == 'missing']


def describe_missing(rows: list[dict]) -> str:
    what = '; '.join(f"{r['label']} {'container' if r['kind'] == 'image' else 'database'} ({r['path']})" for r in rows)
    return (f'Not where your config points, and this run needs them: {what}. '
            'Set them up on Install, under "Tools and reference data", or leave those tools out of this run.')


# ---------------------------------------------------------------------------
# Setting up
# ---------------------------------------------------------------------------

def _sizes(conn, paths: list[str]) -> dict[str, int]:
    if not paths:
        return {}
    code, out = _run(conn, 'du -sbL ' + ' '.join(shlex.quote(p) for p in paths) + ' 2>/dev/null', timeout=600)
    sizes = {}
    for line in out.splitlines():
        size, _, path = line.partition('\t')
        if size.isdigit():
            sizes[path] = int(size)
    return sizes


def plan(conn, cfg: dict, home: str, include_builds: bool = True, accept: set[str] | None = None,
         include_optional: bool = False) -> dict:
    """What setting up would do, item by item, and what it cannot do and why.
    GTDB-Tk and the LLM layer (hundreds of GB between them) only when asked.

    The same answer the Yes / No question shows and the job then runs, so it
    is worked out again when Yes is pressed."""
    accept = accept or set()
    st = status(conn, cfg, home)
    root = user_stores.scratch_root(conn, cfg)
    rows = st['rows']
    todo = [r for r in rows if r['status'] != 'ok']
    items: list[dict] = []
    skipped: list[dict] = []
    config: dict = {}

    def skip(r: dict, reason: str) -> None:
        skipped.append({'id': r['id'], 'label': r['label'], 'kind': r['kind'], 'path': r['path'], 'reason': reason})

    def usable(r: dict) -> bool:
        if r['optional'] and not include_optional:
            skip(r, 'only needed when chosen for a run, and left out here')
            return False
        if not r['source']:
            skip(r, 'nothing to copy it from: not in the lab folder, no container registry set, and margie-build has no recipe for it')
            return False
        if r['source'] == 'build' and not include_builds:
            skip(r, 'only margie-build can make it, and building was left out')
            return False
        if r['source'] == 'build' and r['gated'] and r['tool'] not in accept:
            skip(r, 'its licence was not accepted, and margie-build builds it only with that')
            return False
        return True

    # Containers: one folder for all of them.
    images = [r for r in todo if r['kind'] == 'image' and usable(r)]
    sif_target = st['sif_dir']
    links: list[dict] = []
    if images and not st['sif_dir_writable']:
        sif_target = f'{root}/{SIF_SUB}'
        config[f'{WORKFLOW}.sif_path'] = sif_target
        fetching = {r['id'] for r in images}
        links = [r for r in rows if r['kind'] == 'image' and r['id'] not in fetching and r['status'] != 'missing']
    for r in images:
        items.append({**_item(r), 'dst': f"{sif_target}/{r['name']}"})
    ready_images = {r['tool'] for r in rows if r['kind'] == 'image' and r['status'] == 'ok'} | {r['tool'] for r in images}

    # Databases: each where the config points, or under scratch with db.<tool>.
    for r in todo:
        if r['kind'] != 'database' or not usable(r):
            continue
        dst = r['path']
        movable = r['writable'] and (r['source'] != 'build' or posixpath.basename(dst) == r['tool'])
        if not movable:
            dst = f"{root}/{DB_SUB}/{r['tool']}"
            config[f"db.{r['tool']}"] = dst
        if r['source'] == 'build' and r['tool'] not in ready_images:
            # build.sh indexes a database inside its tool's container.
            skip(r, 'margie-build needs its container first, and that cannot be set up either')
            config.pop(f"db.{r['tool']}", None)
            continue
        items.append({**_item(r), 'dst': dst})

    sizes = _sizes(conn, [i['src'] for i in items if i['action'] == 'copy'])
    for i in items:
        i['bytes'] = sizes.get(i['src'], 0) if i['action'] == 'copy' else None
    total = sum(i['bytes'] or 0 for i in items)
    free = user_stores.depot_free(conn, root) if items else None
    return {
        'root': root,
        'sif_dir': sif_target,
        'items': items,
        'links': [{'src': r['path'], 'dst': f"{sif_target}/{r['name']}"} for r in links],
        'skipped': skipped,
        'config': config,
        'copy_bytes': total,
        'free': free,
        'fits': free is None or free - MARGIN >= total,
        'builds': sum(1 for i in items if i['action'] == 'build'),
        'gated': sorted({r['tool'] for r in todo if r['source'] == 'build' and r['gated']}),
        'statement': LICENCE_STATEMENT,
    }


def _item(r: dict) -> dict:
    return {'id': r['id'], 'label': r['label'], 'kind': r['kind'], 'tool': r['tool'], 'action': r['source'],
            'src': r['shared'] if r['source'] == 'copy' else r['url'], 'gated': r['gated']}


def start(conn, cfg: dict, home: str, include_builds: bool = True, accept: set[str] | None = None,
          include_optional: bool = False) -> dict:
    """Set up what the plan can: one SLURM job, the copy job's progress bar."""
    op = user_stores.progress(conn, cfg)
    if op and op.get('state') in ('queued', 'running'):
        raise StoreError('A copy is already under way.', 409)
    p = plan(conn, cfg, home, include_builds, accept, include_optional)
    if not p['items']:
        raise StoreError('Nothing here can be set up: ' + '; '.join(f"{s['label']}: {s['reason']}" for s in p['skipped'])
                         if p['skipped'] else 'Everything is already there.', 409)
    if not p['fits']:
        gb = lambda n: f'{n / 1024 ** 3:.1f} GB'
        raise StoreError(f"Scratch has {gb(p['free'])} left and the copies need {gb(p['copy_bytes'])} "
                         '(with 5 GB kept free): they would not fit.', 409)
    account = str(_cfg_get(cfg, 'compute.cluster_default.account') or '').strip()
    if p['builds'] and not account:
        raise StoreError('Building needs a SLURM account (Settings): it is too heavy for a login node.', 409)

    steps: list[list] = []
    if p['config'].get(f'{WORKFLOW}.sif_path'):
        steps.append(['mkdir', p['sif_dir']])
    for link in p['links']:
        steps.append(['link', link['src'], link['dst']])
    if p['builds']:
        steps.append(['repo', str(_cfg_get(cfg, f'{WORKFLOW}.build_repo') or BUILD_REPO).strip()])
    for i in p['items']:
        what = 'container' if i['kind'] == 'image' else 'database'
        if i['action'] == 'copy':
            steps.append(['copy', f"{i['label']} {what}, from the lab folder", i['src'], i['dst'], i['bytes'] or 0])
        elif i['action'] == 'pull':
            steps.append(['pull', f"{i['label']} container", i['src'], i['dst']])
        else:
            mode = 'containers' if i['kind'] == 'image' else 'databases'
            name = i['dst'].rsplit('/', 1)[-1]
            tool = name[:-4] if i['kind'] == 'image' else i['tool']
            steps.append(['build', f"Building the {i['label']} {what}", mode, tool, p['sif_dir'],
                          posixpath.dirname(i['dst']) if i['kind'] == 'database' else f"{p['root']}/{DB_SUB}",
                          LICENCE_STATEMENT if i['gated'] else ''])
    # Downloading and indexing a database can take most of a day, and memory.
    heavy = p['builds'] > 0
    user_stores._start(conn, p['root'], 'assets', steps, {'kind': 'assets', 'config': p['config']}, cfg,
                       walltime=str(_cfg_get(cfg, f'{WORKFLOW}.setup_walltime') or ('24:00:00' if heavy else '12:00:00')),
                       mem='32G' if heavy else '4G')
    return {'op': user_stores.progress(conn, cfg), 'plan': p}
