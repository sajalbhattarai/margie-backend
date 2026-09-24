"""
What makes two genomes, or two proteins, the same for MARGIE's caches.

A genome used to be identified by the SHA-256 of its FASTA file's bytes, so
the same assembly saved with a different line width, in lower case, with
Windows line endings or with longer header descriptions was a different
genome, and was annotated again. genome_hash() is taken over the sequence
instead: each record's id (the header's first word) and its bases, upper
case, with every line break and space removed, in file order.

The ids and their order stay in on purpose. RASTtk numbers its features by
contig order and names them after the contig ids, and a whole-genome cache
hit restores RASTtk-numbered tables as they are -- a renamed or reordered
genome must not match them. Such a genome is still recognised protein by
protein (protein_cache.py), where no numbering is involved.

legacy_file_hash() is the old byte hash, kept so that everything already
cached or recorded under it is still found (see genome_hashes()).

protein_hash() is a protein's identity for the per-protein cache: the
SHA-256 of its residues, upper case, whitespace and a final stop (*) removed.

Only the standard library is used.
"""
from __future__ import annotations

import gzip
import hashlib
from pathlib import Path
from typing import IO, Iterator

# Marks a sequence hash, so it is never mistaken for a byte hash of a file.
GENOME_HASH_PREFIX = 'fa1:'


def _open_text(path: str | Path) -> IO[str]:
    p = str(path)
    if p.endswith('.gz'):
        return gzip.open(p, 'rt', encoding='utf-8', errors='replace', newline=None)
    return open(p, 'r', encoding='utf-8', errors='replace', newline=None)


def genome_hash(path: str | Path) -> str:
    """The genome's identity: its records' ids and bases, whatever the file's formatting."""
    sha = hashlib.sha256()
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith('>'):
                words = line[1:].split()
                sha.update(b'>' + (words[0] if words else '').encode() + b'\n')
            else:
                seq = ''.join(line.split()).upper()
                if seq:
                    sha.update(seq.encode())
    return GENOME_HASH_PREFIX + sha.hexdigest()


def legacy_file_hash(path: str | Path) -> str:
    """The SHA-256 of the file's bytes (how genomes were identified before genome_hash)."""
    sha = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            sha.update(chunk)
    return sha.hexdigest()


def genome_hashes(path: str | Path) -> list[str]:
    """[genome_hash, legacy_file_hash]: look records up under either, write new ones under the first."""
    return [genome_hash(path), legacy_file_hash(path)]


def normalise_protein(seq: str) -> str:
    return ''.join(seq.split()).upper().rstrip('*')


def protein_hash(seq: str) -> str:
    return hashlib.sha256(normalise_protein(seq).encode()).hexdigest()


def read_fasta(path: str | Path) -> Iterator[tuple[str, str, str]]:
    """(id, header line without '>', sequence) for each record."""
    header = None
    parts: list[str] = []
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if header is not None:
                    words = header.split()
                    yield (words[0] if words else ''), header, ''.join(parts)
                header = line[1:].rstrip('\n')
                parts = []
            else:
                parts.append(''.join(line.split()))
    if header is not None:
        words = header.split()
        yield (words[0] if words else ''), header, ''.join(parts)
