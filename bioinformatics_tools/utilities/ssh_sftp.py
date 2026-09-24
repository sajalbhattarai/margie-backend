"""
SFTP file operations over SSH.

Provides remote directory listing, file streaming, and YAML config
read/write via paramiko SFTP, plus a paginated line-range reader (via
SSH exec, not SFTP) for the file viewer.

All functions are API-layer only. Pass a per-user SSHConnection built
with make_user_connection() for every call.
"""
import logging
import re
import shlex
import stat
import time

import yaml

from bioinformatics_tools.utilities.ssh_connection import SSHConnection

LOGGER = logging.getLogger(__name__)

_PAGE_SENTINEL = "___MARGIE_PAGE_SENTINEL___"


def list_remote_dir(
    remote_path: str,
    connection: SSHConnection,
) -> list[dict]:
    """List files and directories in a remote path via SFTP.

    Returns a list of dicts: {name, type, size, mtime}; mtime is seconds
    since the epoch, or None when the server does not report it.
    """
    ssh = connection.connect()
    sftp = ssh.open_sftp()
    entries = []
    for attr in sftp.listdir_attr(remote_path):
        entry_type = 'directory' if stat.S_ISDIR(attr.st_mode) else 'file'
        entries.append({
            'name': attr.filename,
            'type': entry_type,
            'size': attr.st_size,
            'mtime': attr.st_mtime,
        })
    sftp.close()
    pass  # pooled client: closing it would defeat SSHConnection's pool
    return entries


def list_remote_dir_checked(
    remote_path: str,
    connection: SSHConnection,
) -> list[dict]:
    """list_remote_dir, with the check that the path is a directory made in
    the same SFTP session (one session per listing, not two: each one is a
    round trip and a new sftp-server on the login node). Raises
    FileNotFoundError when the path is not there and NotADirectoryError when
    it is a file."""
    ssh = connection.connect()
    sftp = ssh.open_sftp()
    try:
        try:
            attr = sftp.stat(remote_path)
        except FileNotFoundError:
            raise FileNotFoundError(f'Path not found on cluster: {remote_path}')
        if not stat.S_ISDIR(attr.st_mode):
            raise NotADirectoryError(remote_path)
        return [{
            'name': a.filename,
            'type': 'directory' if stat.S_ISDIR(a.st_mode) else 'file',
            'size': a.st_size,
            'mtime': a.st_mtime,
        } for a in sftp.listdir_attr(remote_path)]
    finally:
        sftp.close()


def stream_remote_file(
    remote_path: str,
    connection: SSHConnection,
):
    """Open a remote file eagerly and return a generator that streams it in 8KB chunks.

    The SFTP open happens at call time (not during iteration), so callers can
    catch FileNotFoundError / IOError before wrapping the result in
    StreamingResponse.  Without this, all SFTP errors happen inside the
    StreamingResponse generator after headers are sent, which drops the
    connection and causes "fetch failed" in the browser.
    """
    ssh = connection.connect()
    sftp = ssh.open_sftp()
    f = sftp.open(remote_path, 'rb')   # raises FileNotFoundError/IOError here if absent

    def _chunks():
        try:
            while True:
                chunk = f.read(8192)
                if not chunk:
                    break
                yield chunk
        finally:
            f.close()
            sftp.close()
            pass  # pooled client: closing it would defeat SSHConnection's pool

    return _chunks()


def read_remote_yaml(
    remote_path: str,
    connection: SSHConnection,
) -> dict:
    """Read and parse a YAML file from the remote cluster.

    Returns the parsed dict, or an empty dict if the file does not exist.
    """
    ssh = connection.connect()
    sftp = ssh.open_sftp()
    try:
        with sftp.open(remote_path, 'r') as f:
            content = f.read().decode('utf-8')
        return yaml.safe_load(content) or {}
    except FileNotFoundError:
        LOGGER.info('Remote config not found at %s — returning empty dict', remote_path)
        return {}
    finally:
        sftp.close()
        pass  # pooled client: closing it would defeat SSHConnection's pool


def check_remote_file(
    path: str,
    connection: SSHConnection,
) -> None:
    """Verify a remote file exists and is a regular file via SFTP.

    Raises FileNotFoundError if the path does not exist on the cluster,
    or IsADirectoryError if it resolves to a directory rather than a file.
    """
    ssh = connection.connect()
    sftp = ssh.open_sftp()
    try:
        attr = sftp.stat(path)
        if stat.S_ISDIR(attr.st_mode):
            raise IsADirectoryError(f'Path is a directory, not a file: {path}')
    except FileNotFoundError:
        raise FileNotFoundError(f'File not found on cluster: {path}')
    finally:
        sftp.close()
        pass  # pooled client: closing it would defeat SSHConnection's pool


def stat_remote_file(path: str, connection: SSHConnection) -> tuple[float, int]:
    """Return (mtime, size) for a remote file via a single SFTP stat call.

    Used to key a cache entry against a remote file's identity, cheaply
    (metadata only, no content read) -- output files don't change once a
    job completes, so (mtime, size) staying the same means a cached
    result (e.g. a line count) is still valid.

    Raises FileNotFoundError if the path does not exist on the cluster.
    """
    ssh = connection.connect()
    sftp = ssh.open_sftp()
    try:
        attr = sftp.stat(path)
    except FileNotFoundError:
        raise FileNotFoundError(f'File not found on cluster: {path}')
    finally:
        sftp.close()
        pass  # pooled client: closing it would defeat SSHConnection's pool
    return (attr.st_mtime, attr.st_size)


def check_remote_path_kind(path: str, connection: SSHConnection) -> str:
    """Check whether a remote path exists and whether it's a file or directory.

    Returns 'file' or 'directory'. Raises FileNotFoundError if the path
    does not exist on the cluster.
    """
    ssh = connection.connect()
    sftp = ssh.open_sftp()
    try:
        attr = sftp.stat(path)
    except FileNotFoundError:
        raise FileNotFoundError(f'Path not found on cluster: {path}')
    finally:
        sftp.close()
        pass  # pooled client: closing it would defeat SSHConnection's pool
    return 'directory' if stat.S_ISDIR(attr.st_mode) else 'file'


def write_remote_yaml(
    remote_path: str,
    data: dict,
    connection: SSHConnection,
) -> None:
    """Write a dict as YAML to a remote path via SFTP.

    Creates parent directories on the remote if they do not exist.
    """
    ssh = connection.connect()

    # Ensure parent directory exists
    parent = remote_path.rsplit('/', 1)[0]
    if parent:
        ssh.exec_command(f'mkdir -p {parent}')

    sftp = ssh.open_sftp()
    content = yaml.dump(data, default_flow_style=False, allow_unicode=True)
    with sftp.open(remote_path, 'w') as f:
        f.write(content)
    sftp.close()
    pass  # pooled client: closing it would defeat SSHConnection's pool
    LOGGER.info('Wrote remote config to %s', remote_path)


def write_remote_text_file(
    remote_path: str,
    content: str,
    connection: SSHConnection,
) -> None:
    """Write raw text to a remote path via SFTP, creating parent directories
    if needed. Unlike write_remote_yaml, writes the content verbatim -- for
    the file explorer's Save action on arbitrary text files (config.yaml,
    scripts, etc.), not just structured config.
    """
    ssh = connection.connect()

    parent = remote_path.rsplit('/', 1)[0]
    if parent:
        ssh.exec_command(f'mkdir -p {shlex.quote(parent)}')

    sftp = ssh.open_sftp()
    try:
        with sftp.open(remote_path, 'w') as f:
            f.write(content)
    finally:
        sftp.close()
        pass  # pooled client: closing it would defeat SSHConnection's pool
    LOGGER.info('Wrote remote text file to %s', remote_path)


def copy_remote_directory(
    src_path: str,
    dest_path: str,
    connection: SSHConnection,
) -> None:
    """Copy a remote directory tree into a new path entirely on the cluster
    filesystem via a single SSH exec_command -- the API process never reads
    or writes the file bytes itself, so this scales to however large a
    job's output_dir is without taxing dane-api's own memory/bandwidth.

    Excludes .snakemake/ (Snakemake's own bookkeeping -- its metadata
    filenames and JSON content embed the OLD absolute output_dir path,
    which would be wrong after the copy; Snakemake regenerates this fully
    fresh on the next invocation against dest_path) and any tool's
    original_container_outputs/*/stage/ subdirectory (pure rule-local
    scratch space, re-staged fresh by the rule itself every run -- copying
    it just wastes bytes on genome FASTA copies that get blown away
    immediately).

    Creates dest_path if needed. Raises RuntimeError if the remote rsync
    command exits non-zero.
    """
    ssh = connection.connect()
    try:
        quoted_src = shlex.quote(src_path.rstrip('/') + '/')
        quoted_dest = shlex.quote(dest_path)
        cmd = (
            f'mkdir -p {quoted_dest} && '
            f'rsync -a '
            f"--exclude='.snakemake' "
            f"--exclude='original_container_outputs/*/stage' "
            f'{quoted_src} {quoted_dest}'
        )
        LOGGER.info('Copying remote directory %s -> %s', src_path, dest_path)
        _, stdout, stderr = ssh.exec_command(cmd)
        exit_code = stdout.channel.recv_exit_status()
        if exit_code != 0:
            err = stderr.read().decode('utf-8', errors='replace')
            raise RuntimeError(f'rsync failed (exit {exit_code}) copying {src_path} to {dest_path}: {err}')
    finally:
        pass  # pooled client: closing it would defeat SSHConnection's pool


def stage_selected_genomes(
    source_dir: str,
    names: list[str],
    connection: SSHConnection,
    label: str = '',
) -> str:
    """Build a folder on the cluster holding only the chosen genomes, and
    return its absolute path.

    A workflow is pointed at a folder and annotates everything in it; there is
    no "run only these" flag anywhere in the chain. So to run a subset, the
    subset is given a folder of its own.

    The entries are symlinks, not copies: a genome is tens of megabytes, the
    run only ever reads them, and hundreds of them would otherwise be
    duplicated on scratch for no reason. Creating them is a single
    exec_command, so the API process never touches the bytes.

    The folder lives under the user's home, which is always writable and
    always present -- unlike scratch, whose path differs between clusters.
    Old selections are left where they are: they are a handful of symlinks
    each, and they record what a finished job was actually given.

    Raises ValueError for a name that is not a plain file name, and
    RuntimeError if the remote command fails.
    """
    if not names:
        raise ValueError('No genomes were selected.')
    for name in names:
        # These come from a browser, so they are checked rather than trusted:
        # anything with a slash in it could reach outside source_dir.
        if not name or '/' in name or name in ('.', '..'):
            raise ValueError(f'Not a genome file name: {name!r}')

    stamp = time.strftime('%Y%m%d-%H%M%S')
    suffix = f"-{re.sub(r'[^A-Za-z0-9]+', '-', label).strip('-')}" if label else ''
    rel = f'.margie/selections/{stamp}{suffix}'

    src = shlex.quote(source_dir.rstrip('/'))
    links = ' && '.join(
        f'ln -sfn {src}/{shlex.quote(n)} "$d"/{shlex.quote(n)}' for n in names
    )
    # printf at the end so the absolute path comes back without a trailing
    # newline to strip off guesswork later.
    cmd = f'd="$HOME"/{shlex.quote(rel)} && mkdir -p "$d" && {links} && printf %s "$d"'

    ssh = connection.connect()
    LOGGER.info('Staging %d selected genome(s) from %s', len(names), source_dir)
    _, stdout, stderr = ssh.exec_command(cmd)
    out = stdout.read().decode('utf-8', errors='replace').strip()
    exit_code = stdout.channel.recv_exit_status()
    if exit_code != 0 or not out:
        err = stderr.read().decode('utf-8', errors='replace')
        raise RuntimeError(
            f'Could not stage the selected genomes on the cluster (exit {exit_code}): {err or "no path returned"}'
        )
    return out


def _build_path_rewrite_script(directory: str, old_path: str, new_path: str) -> str:
    """Buildsd the embedded Python source run remotely by
    rewrite_path_references(). Separated out so tests can run this exact
    script locally against a real temp directory, proving the
    find-and-replace logic itself works, not just that some command got
    sent over SSH.

    Values are embedded via repr() (not f-string interpolation of the raw
    string) so paths containing quotes/backslashes/etc. round-trip as
    correct Python source -- the script itself then does a pure literal
    str.replace(), so the path values never need regex/sed escaping at all.
    """
    return (
        "import os\n"
        f"old = {old_path!r}\n"
        f"new = {new_path!r}\n"
        "modified = 0\n"
        f"for root, _, files in os.walk({directory!r}):\n"
        "    for name in files:\n"
        "        path = os.path.join(root, name)\n"
        "        try:\n"
        "            with open(path, 'r', encoding='utf-8') as f:\n"
        "                content = f.read()\n"
        "        except (UnicodeDecodeError, OSError):\n"
        "            continue\n"
        "        if old not in content:\n"
        "            continue\n"
        "        try:\n"
        "            with open(path, 'w', encoding='utf-8') as f:\n"
        "                f.write(content.replace(old, new))\n"
        "            modified += 1\n"
        "        except OSError:\n"
        "            continue\n"
        "print(modified)\n"
    )


def rewrite_path_references(
    directory: str,
    old_path: str,
    new_path: str,
    connection: SSHConnection,
) -> int:
    """Find-and-replace every literal occurrence of old_path with new_path
    across every text file under directory, in place.

    Used after copy_remote_directory() during Resume to fix up provenance
    metadata that some tools (GTDB-Tk, KEGG, ...) stamp their own
    invocation path into their results.tsv content -- confirmed cosmetic
    only (nothing downstream reads those columns), but cheap to correct so
    a resumed run's output doesn't carry stale path strings pointing at a
    directory that no longer matches where the file actually lives.

    Runs as a small embedded Python script over SSH exec rather than shell
    grep/sed specifically to avoid escaping headaches -- paths can contain
    '.', '-', and other characters that are regex/sed metacharacters;
    Python's str.replace() is a pure literal substitution with no such
    concerns. Binary files (anything that fails UTF-8 decoding) are
    skipped, not corrupted. Best-effort: a single file's read/write error
    is skipped, not fatal to the whole pass -- this is a cosmetic cleanup,
    not something that should ever block a resumed job from launching.

    Returns the number of files modified.
    """
    ssh = connection.connect()
    try:
        script = _build_path_rewrite_script(directory, old_path, new_path)
        remote_python = '~/bioinformatics-tools/.venv/bin/python'
        cmd = f'{remote_python} -c {shlex.quote(script)}'
        _, stdout, stderr = ssh.exec_command(cmd)
        out = stdout.read().decode().strip()
        exit_code = stdout.channel.recv_exit_status()
        if exit_code != 0:
            err = stderr.read().decode('utf-8', errors='replace')
            LOGGER.warning('rewrite_path_references failed (exit %d) for %s: %s', exit_code, directory, err)
            return 0
        return int(out) if out.isdigit() else 0
    finally:
        pass  # pooled client: closing it would defeat SSHConnection's pool


def read_remote_file_page(
    remote_path: str,
    start_row: int,
    end_row: int,
    connection: SSHConnection,
    known_total_lines: int | None = None,
) -> dict:
    """Reads a 1-indexed inclusive line range [start_row, end_row] from a
    remote text file, plus its header line (line 1) and total line count,
    in a single SSH exec_command.

    There is no connection pooling in this codebase (connect() opens a
    fresh SSH session every call), so this combines the wc -l / header /
    range reads into one remote shell invocation instead of three, to
    avoid paying the handshake cost more than once per page request.

    Uses SSH exec (sed/wc) rather than SFTP byte-seeking: paramiko's SFTP
    has no line-counting primitive, and remote sed/wc scan a file far
    faster than reading it chunk-by-chunk over SFTP to count newlines
    ourselves would. The range read quits as soon as it passes end_row
    (sed's `q` command), so cost scales with how far into the file a page
    is, not with total file size.

    If known_total_lines is given, skips the wc -l pass entirely (caller
    already has a cached, still-valid count for this file).

    Returns {"total_lines": int, "header": str, "lines": list[str]}.
    Raises FileNotFoundError if remote_path does not exist, or
    IsADirectoryError if it resolves to a directory.
    """
    check_remote_file(remote_path, connection)

    quoted = shlex.quote(remote_path)
    parts = []
    if known_total_lines is None:
        parts.append(f'wc -l < {quoted}')
    # `1{p;q}` prints line 1 then quits immediately, regardless of file size.
    parts.append(f"sed -n '1{{p;q}}' {quoted}")
    # Print MUST come before quit: sed's `q` exits immediately (skipping any
    # later command in that cycle), so `end_row p` before `end_row q` is
    # required or the last row of the page would be silently dropped.
    parts.append(f"sed -n '{start_row},{end_row}p;{end_row}q' {quoted}")
    script = f'; echo {_PAGE_SENTINEL}; '.join(parts)

    ssh = connection.connect()
    try:
        _, stdout, stderr = ssh.exec_command(script)
        output = stdout.read().decode('utf-8', errors='replace')
        err = stderr.read().decode('utf-8', errors='replace').strip()
        if err:
            LOGGER.warning('read_remote_file_page stderr for %s: %s', remote_path, err)
    finally:
        pass  # pooled client: closing it would defeat SSHConnection's pool

    sections = output.split(f'{_PAGE_SENTINEL}\n')
    if known_total_lines is None:
        total_lines = int(sections[0].strip() or 0)
        header = sections[1].rstrip('\n')
        data_block = sections[2]
    else:
        total_lines = known_total_lines
        header = sections[0].rstrip('\n')
        data_block = sections[1]

    lines = data_block.split('\n')
    if lines and lines[-1] == '':
        lines.pop()

    return {"total_lines": total_lines, "header": header, "lines": lines}
