"""Compatibility wrapper for isolated v0.2 tests and optional live doctor checks."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from bur01_structmine.doctor import print_doctor, run_doctor

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--with-ncbi-network", action="store_true")
    parser.add_argument("--with-alphafold-network", action="store_true")
    args = parser.parse_args()
    status = subprocess.run([sys.executable, "-m", "pytest", "-q"], cwd=PROJECT_ROOT).returncode
    if args.with_ncbi_network or args.with_alphafold_network:
        report = run_doctor("config.yaml", check_network=True)
        print_doctor(report)
        status = status or (0 if report.passed else 1)
    return status


if __name__ == "__main__":
    raise SystemExit(main())
