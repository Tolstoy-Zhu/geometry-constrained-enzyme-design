"""Command-line entry point for AquaEsterNet-StructMine.

Run examples:
  python run_pipeline.py --config config.yaml
  python run_pipeline.py --config config.yaml --skip-foldseek --input-tsv data/example_foldseek_result.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Make the local src/ folder importable in PyCharm and command line.
PROJECT_ROOT = Path(__file__).resolve().parent
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from bur01_structmine.pipeline import run_pipeline  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description="AquaEsterNet-StructMine pipeline")
    parser.add_argument("--config", default="config.yaml", help="Path to config.yaml")
    parser.add_argument(
        "--skip-foldseek",
        action="store_true",
        help="Skip Foldseek and use an existing TSV result file",
    )
    parser.add_argument(
        "--input-tsv",
        default=None,
        help="Existing Foldseek TSV file when --skip-foldseek is used",
    )
    parser.add_argument("--run-id", default=None, help="Optional stable run identifier")
    parser.add_argument("--output-root", default=None, help="Versioned run directory root")
    parser.add_argument(
        "--legacy-output",
        action="store_true",
        help="Also use the configured fixed output prefix; intended only for compatibility checks",
    )
    args = parser.parse_args()

    summary = run_pipeline(
        config_path=args.config,
        skip_foldseek=args.skip_foldseek,
        input_tsv=args.input_tsv,
        run_id=args.run_id,
        output_root=args.output_root,
        legacy_output=args.legacy_output,
    )
    print(f"Manifest: {summary.manifest_path}")


if __name__ == "__main__":
    main()
