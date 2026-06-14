"""Demo script that does not require Foldseek or internet.

It only extracts the amino-acid sequence from data/query/Bur01m.pdb and writes
results/demo_query_sequence.fasta. This is useful for checking whether PyCharm
and dependencies are installed correctly.
"""
from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from bur01_structmine.io_utils import write_fasta  # noqa: E402
from bur01_structmine.pdb_utils import extract_sequence_from_pdb  # noqa: E402


def main() -> None:
    pdb_file = PROJECT_ROOT / "data" / "query" / "Bur01m.pdb"
    out_file = PROJECT_ROOT / "results" / "demo_query_sequence.fasta"

    seqs = extract_sequence_from_pdb(pdb_file)
    records = [(f"Bur01m_chain_{chain}", seq) for chain, seq in seqs.items()]
    write_fasta(records, out_file)

    print("[OK] Demo finished.")
    print(f"Extracted chains: {list(seqs.keys())}")
    print(f"Output FASTA: {out_file}")


if __name__ == "__main__":
    main()
