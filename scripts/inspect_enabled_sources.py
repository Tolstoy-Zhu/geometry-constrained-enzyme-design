"""Inspect enabled AquaEsterNet-StructMine search sources in a config file."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from bur01_structmine.config import load_config  # noqa: E402
from bur01_structmine.datasources import load_search_databases  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Inspect enabled search sources")
    parser.add_argument("--config", default="config.yaml", help="Path to config.yaml or a preset config")
    parser.add_argument("--require-source", action="append", default=[], help="Required enabled source_db name")
    parser.add_argument("--json", action="store_true", help="Write machine-readable JSON")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    sources = load_search_databases(cfg)
    rows = [
        {
            "name": source.name,
            "kind": source.kind,
            "annotation_source": source.annotation_source,
            "target_db": str(source.target_db or ""),
            "search_engine": source.search_engine,
            "foldseek_database": source.foldseek_database,
            "input_fasta": str(source.input_fasta or ""),
            "structure_dir": str(source.structure_dir or ""),
        }
        for source in sources
    ]

    names = {row["name"] for row in rows}
    missing = [name for name in args.require_source if name not in names]
    if missing:
        raise SystemExit(f"Missing required enabled source(s): {', '.join(missing)}")

    if args.json:
        print(json.dumps(rows, indent=2, ensure_ascii=False))
    else:
        for row in rows:
            print(
                f"{row['name']}\t{row['kind']}\t{row['annotation_source']}\t"
                f"{row['search_engine']}\t{row['target_db']}\t"
                f"{row['foldseek_database'] or row['input_fasta'] or row['structure_dir']}"
            )
        print(f"[OK] Enabled sources: {', '.join(sorted(names))}")


if __name__ == "__main__":
    main()
