"""
Local mode: the API runs on the user's own computer and works there instead
of on an HPC cluster over SSH.

The GUI's "Use MARGIE locally" option starts this API with BSP_LOCAL_MODE=1,
bound to 127.0.0.1. Nothing else sets it, so a hosted deployment keeps
requiring cluster credentials.
"""
import getpass
import os

# Stored as cluster_host for local accounts, so they are recognizable in the
# users table (whose cluster columns are NOT NULL).
LOCAL_CLUSTER_HOST = 'localhost'


def is_local_mode() -> bool:
    return os.getenv('BSP_LOCAL_MODE') == '1'


def local_account() -> dict:
    """The cluster-side fields for a local account: this machine and user."""
    return {
        'cluster_host': LOCAL_CLUSTER_HOST,
        'cluster_username': getpass.getuser(),
        'home_dir': os.path.expanduser('~'),
    }
