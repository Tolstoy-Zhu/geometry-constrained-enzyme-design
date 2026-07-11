"""Compatibility launcher for common AquaEsterNet-StructMine actions."""
from __future__ import annotations

import argparse

from bur01_structmine.cli import main as cli_main


def _overview() -> int:
    print("AquaEsterNet-StructMine v0.3")
    print("Local AlphaFold2 prediction and Bur01-like structure mining with Foldseek, TM-align, and WASE-screen.")
    print()
    print("Common actions:")
    print("  python start_here.py --action doctor")
    print("  python start_here.py --action ui")
    print("  python start_here.py --action benchmark")
    print("  python start_here.py --action production-run")
    print("  python start_here.py --action validate")
    print("  python start_here.py --action check-files")
    print("  python -m bur01_structmine.cli predict --fasta query.fasta --preset standard --confirm-public-msa")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="AquaEsterNet-StructMine launcher")
    parser.add_argument(
        "--action",
        choices=(
            "overview",
            "check-files",
            "list-sources",
            "quick-test",
            "doctor",
            "ui",
            "benchmark",
            "target-run",
            "production-run",
            "validate",
            "validate-live",
        ),
        default="overview",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.action == "overview":
        return _overview()
    mapping = {
        "check-files": ["check-files"],
        "list-sources": ["list-sources"],
        "quick-test": ["validate"],
        "doctor": ["doctor"],
        "ui": ["ui"],
        "benchmark": ["benchmark", "--engine", "foldseek"],
        "target-run": ["benchmark", "--engine", "foldseek"],
        "production-run": ["production-run"],
        "validate": ["validate"],
        "validate-live": ["validate-live"],
    }
    command = mapping[args.action]
    if args.dry_run:
        print("[DRY-RUN] aquaesternet-structmine " + " ".join(command))
        return 0
    return cli_main(command)


if __name__ == "__main__":
    raise SystemExit(main())
