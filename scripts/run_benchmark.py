"""Run the fixed benchmark with Foldseek or the TM-align preview engine."""
from __future__ import annotations

import argparse
import json

from bur01_structmine.benchmark import run_fixed_benchmark


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--engine", choices=("foldseek", "tmtools"), default="foldseek")
    parser.add_argument("--output-root")
    parser.add_argument("--no-prepare", action="store_true")
    args = parser.parse_args()
    _, report, report_path = run_fixed_benchmark(
        engine=args.engine,
        output_root=args.output_root,
        prepare=not args.no_prepare,
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))
    print(f"Benchmark report: {report_path}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
