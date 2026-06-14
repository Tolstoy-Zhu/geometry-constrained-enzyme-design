"""Foldseek command-line wrapper.

Foldseek is not a Python package here. This module calls the Foldseek executable
installed in Linux, WSL2, or a cloud server.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Iterable, List

from .io_utils import ensure_dir

# These columns are requested from Foldseek. If your Foldseek version does not
# support one field, remove it from this list and rerun.
FOLDSEEK_FIELDS: List[str] = [
    "query",
    "target",
    "evalue",
    "bits",
    "prob",
    "alntmscore",
    "qtmscore",
    "ttmscore",
    "rmsd",
    "fident",
    "alnlen",
    "qlen",
    "tlen",
    "qcov",
    "tcov",
]


def check_foldseek_executable(executable: str = "foldseek") -> bool:
    """Return True if Foldseek executable is available."""
    return shutil.which(executable) is not None


def run_command(cmd: Iterable[str]) -> None:
    """Run a shell command and raise an error if it fails."""
    cmd = [str(x) for x in cmd]
    print("[CMD]", " ".join(cmd))
    subprocess.run(cmd, check=True)


def run_easy_search(
    query_pdb: str | Path,
    target_db: str | Path,
    out_tsv: str | Path,
    tmp_dir: str | Path,
    executable: str = "foldseek",
    max_seqs: int = 1000,
    threads: int = 8,
) -> Path:
    """Run Foldseek easy-search.

    Parameters
    ----------
    query_pdb:
        Input PDB/mmCIF file, e.g. Bur01.pdb.
    target_db:
        Foldseek database path, e.g. databases/pdb100.
    out_tsv:
        Output tab-separated result file.
    tmp_dir:
        Temporary folder.
    """
    query_pdb = Path(query_pdb)
    target_db = Path(target_db)
    out_tsv = Path(out_tsv)
    tmp_dir = Path(tmp_dir)

    ensure_dir(out_tsv.parent)
    ensure_dir(tmp_dir)

    if not query_pdb.exists():
        raise FileNotFoundError(f"Query PDB not found: {query_pdb}")

    if not check_foldseek_executable(executable):
        raise RuntimeError(
            f"Foldseek executable not found: {executable}. "
            "Install Foldseek first or run in demo mode."
        )

    format_output = ",".join(FOLDSEEK_FIELDS)
    cmd = [
        executable,
        "easy-search",
        query_pdb,
        target_db,
        out_tsv,
        tmp_dir,
        "--max-seqs",
        str(max_seqs),
        "--threads",
        str(threads),
        "--format-output",
        format_output,
    ]
    run_command(cmd)
    return out_tsv
