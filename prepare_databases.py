"""Prepare search databases configured for AquaEsterNet-StructMine.

Examples:
  python prepare_databases.py --dry-run
  python prepare_databases.py --only alphafold_swissprot
  python prepare_databases.py --only ncbi_prostt5
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from bur01_structmine.pipeline import prepare_search_databases  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare configured Foldseek databases")
    parser.add_argument("--config", default="config.yaml", help="Path to config.yaml")
    parser.add_argument(
        "--only",
        nargs="+",
        default=None,
        help="Prepare only these enabled search database names",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned database preparation commands without running Foldseek",
    )
    args = parser.parse_args()

    selected = prepare_search_databases(
        config_path=args.config,
        only=set(args.only) if args.only else None,
        dry_run=args.dry_run,
    )
    names = ", ".join(source.name for source in selected)
    if args.dry_run:
        print(f"[OK] Dry-run inspected: {names}")
    else:
        print(f"[OK] Prepared: {names}")


if __name__ == "__main__":
    main()
