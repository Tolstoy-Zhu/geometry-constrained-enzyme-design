"""Input/output helper functions."""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Tuple

import pandas as pd


def ensure_dir(path: str | Path) -> Path:
    """Create a directory if it does not exist and return it as Path."""
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def write_fasta(records: Iterable[Tuple[str, str]], out_file: str | Path) -> None:
    """Write FASTA records.

    Parameters
    ----------
    records:
        Iterable of (record_id, sequence).
    out_file:
        Output FASTA path.
    """
    out_file = Path(out_file)
    ensure_dir(out_file.parent)
    with out_file.open("w", encoding="utf-8") as f:
        for rec_id, seq in records:
            seq = str(seq).replace(" ", "").replace("\n", "").strip()
            if not rec_id or not seq or seq.lower() == "nan":
                continue
            f.write(f">{rec_id}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")


def write_excel_safely(df: pd.DataFrame, out_file: str | Path) -> None:
    """Write an Excel file and create its parent directory."""
    out_file = Path(out_file)
    ensure_dir(out_file.parent)
    df.to_excel(out_file, index=False)
