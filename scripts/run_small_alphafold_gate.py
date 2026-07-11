"""Compatibility wrapper for the fixed Bur01-like benchmark.

The former P00520-only gate was removed in v0.2. P00520 is now a remote
negative control and cannot serve as evidence that the workflow succeeded.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the fixed Bur01-like benchmark")
    parser.add_argument("--engine", choices=("foldseek", "tmtools"), default="foldseek")
    parser.add_argument("--dry-run", action="store_true")
    args, legacy = parser.parse_known_args()
    if legacy:
        print(f"[WARN] Ignoring removed v0.1 arguments: {' '.join(legacy)}")
    command = [sys.executable, "scripts/run_benchmark.py", "--engine", args.engine]
    if args.dry_run:
        print("[DRY-RUN] " + " ".join(command))
        return 0
    return subprocess.run(command, cwd=PROJECT_ROOT).returncode


if __name__ == "__main__":
    raise SystemExit(main())
