"""Unified command-line interface for AquaEsterNet-StructMine."""
from __future__ import annotations

import argparse
import json
import socket
import subprocess
import sys
from pathlib import Path

from .benchmark import run_fixed_benchmark
from .config import load_config, project_root
from .datasources import load_search_databases
from .doctor import print_doctor, run_doctor
from .models import PredictionRequest, RunRequest
from .pipeline import run_pipeline
from .prediction import run_prediction, run_prediction_and_mine

PROJECT_NAME = "AquaEsterNet-StructMine"
REQUIRED_FILES = (
    "pyproject.toml",
    "requirements.txt",
    "requirements-dev.txt",
    "requirements.lock.txt",
    "README_FIRST.txt",
    "README.md",
    "README_CN.md",
    "MODEL_GUIDE_CN.md",
    "MODEL_GUIDE_EN.md",
    "CHANGELOG.md",
    "THIRD_PARTY_LICENSES.md",
    "docs/RELEASE_NOTES_CN.md",
    "docs/RELEASE_NOTES_EN.md",
    "docs/validation/v0.3_validation_summary.json",
    "docs/screenshots/af_prediction_desktop.png",
    "docs/screenshots/af_prediction_mobile.png",
    "config.yaml",
    "configs/benchmark_foldseek.yaml",
    "configs/prediction_runtime.lock.yaml",
    "configs/colabfold-pip-lock.txt",
    ".github/workflows/ci.yml",
    ".streamlit/config.toml",
    "data/query/Bur01m.pdb",
    "data/query/Bur01m_seed.fasta",
    "data/benchmark/benchmark_manifest.yaml",
    "data/benchmark/README.md",
    "run_pipeline.py",
    "start_here.py",
    "start_ui.bat",
    "ui/app.py",
    "src/bur01_structmine/pipeline.py",
    "src/bur01_structmine/structure_search.py",
    "src/bur01_structmine/scoring.py",
    "src/bur01_structmine/doctor.py",
    "src/bur01_structmine/benchmark.py",
    "src/bur01_structmine/models.py",
    "src/bur01_structmine/prediction.py",
    "src/bur01_structmine/prediction_sequence.py",
    "src/bur01_structmine/colabfold_runtime.py",
    "scripts/prepare_benchmark.py",
    "scripts/run_benchmark.py",
    "scripts/setup_wsl_foldseek.ps1",
    "scripts/run_production_acceptance.ps1",
    "scripts/write_database_manifest.py",
    "scripts/build_verified_archive.py",
    "scripts/install_colabfold_wsl.sh",
    "scripts/patch_colabfold_result_timeout.py",
    "scripts/install_alphafold_params.sh",
    "scripts/run_af2_acceptance.py",
    "scripts/build_public_release.py",
    "scripts/validate_ui_playwright.py",
    "tests/fixtures/fake_colabfold.py",
    "tests/test_prediction.py",
    "tests/test_prediction_sequence.py",
)


def _free_port(preferred: int = 8501) -> int:
    for port in range(preferred, preferred + 100):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            try:
                sock.bind(("127.0.0.1", port))
            except OSError:
                continue
            return port
    raise RuntimeError("No free local port was found in the 8501-8600 range")


def check_files() -> int:
    root = project_root()
    missing = [path for path in REQUIRED_FILES if not (root / path).exists()]
    print(f"{PROJECT_NAME} file check")
    print(f"Project root: {root}")
    print(f"Required files: {len(REQUIRED_FILES)}")
    if missing:
        print("[FAIL] Missing files:")
        for path in missing:
            print(f"  - {path}")
        return 1
    print("[OK] Required project files are present.")
    return 0


def list_sources(config_paths: list[str]) -> int:
    for config_path in config_paths:
        cfg = load_config(config_path)
        print(f"\n{config_path}")
        for source in load_search_databases(cfg):
            path = source.target_db or source.structure_dir or source.input_fasta or source.result_tsv or ""
            print(f"  {source.name}: engine={source.search_engine}, kind={source.kind}, path={path}")
    return 0


def launch_ui(port: int | None = None, no_browser: bool = False) -> int:
    root = project_root()
    app = root / "ui" / "app.py"
    if not app.exists():
        raise FileNotFoundError(f"Streamlit app not found: {app}")
    selected_port = port or _free_port()
    url = f"http://127.0.0.1:{selected_port}"
    print(f"Starting {PROJECT_NAME} at {url}")
    command = [
        sys.executable,
        "-m",
        "streamlit",
        "run",
        str(app),
        "--server.address",
        "127.0.0.1",
        "--server.port",
        str(selected_port),
        "--browser.gatherUsageStats",
        "false",
        "--server.headless",
        "true" if no_browser else "false",
    ]
    return subprocess.run(command, cwd=root).returncode


def _run_tests(extra_args: list[str] | None = None) -> int:
    command = [sys.executable, "-m", "pytest", "-q", *(extra_args or [])]
    return subprocess.run(command, cwd=project_root()).returncode


def _add_prediction_arguments(parser: argparse.ArgumentParser) -> None:
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--fasta", help="Path to one 20-1000 aa monomer FASTA")
    source.add_argument("--candidate-run", help="Run ID containing an annotated candidate sequence")
    parser.add_argument("--accession", help="Candidate accession; required with --candidate-run")
    parser.add_argument("--preset", choices=("quick", "standard", "high"), default="standard")
    parser.add_argument("--templates", action="store_true", help="Enable public template search")
    parser.add_argument(
        "--confirm-public-msa",
        action="store_true",
        help="Confirm that a cache miss may send the sequence to the public ColabFold MSA service",
    )
    parser.add_argument("--prediction-id")
    parser.add_argument("--prediction-output-root", default="results/predictions")
    parser.add_argument("--runtime-backend", choices=("wsl", "native"), default="wsl")
    parser.add_argument("--wsl-distribution", default="Ubuntu-Foldseek")
    parser.add_argument("--colabfold-executable", default="/opt/colabfold/bin/colabfold_batch")
    parser.add_argument(
        "--catalytic-site",
        action="append",
        default=[],
        metavar="LABEL=POSITION",
        help="Override an auto-mapped site, for example S111=100",
    )


def _prediction_request(args: argparse.Namespace) -> PredictionRequest:
    if args.candidate_run and not args.accession:
        raise ValueError("--accession is required with --candidate-run")
    catalytic_override: dict[str, int] = {}
    for value in args.catalytic_site:
        try:
            label, position = value.split("=", 1)
            catalytic_override[label.strip().upper()] = int(position)
        except (ValueError, TypeError) as exc:
            raise ValueError(f"Invalid --catalytic-site value: {value!r}") from exc
    return PredictionRequest(
        fasta_path=args.fasta,
        candidate_run_id=args.candidate_run,
        accession=args.accession,
        preset=args.preset,
        use_templates=args.templates,
        confirm_public_msa=args.confirm_public_msa,
        prediction_id=args.prediction_id,
        output_root=args.prediction_output_root,
        runtime_backend=args.runtime_backend,
        wsl_distribution=args.wsl_distribution,
        colabfold_executable=args.colabfold_executable,
        catalytic_override=catalytic_override or None,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="aquaesternet-structmine", description=PROJECT_NAME)
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("check-files", help="Verify the project delivery file set")
    sources = subparsers.add_parser("list-sources", help="Inspect enabled data sources")
    sources.add_argument("--config", action="append", default=[])

    doctor = subparsers.add_parser("doctor", help="Diagnose the local runtime and data environment")
    doctor.add_argument("--config", default="config.yaml")
    doctor.add_argument("--json", action="store_true", dest="as_json")
    doctor.add_argument("--no-network", action="store_true")
    doctor.add_argument("--strict", action="store_true", help="Return nonzero when required checks fail")

    run = subparsers.add_parser("run", help="Run the shared structure-mining pipeline")
    run.add_argument("--config", default="config.yaml")
    run.add_argument("--run-id")
    run.add_argument("--output-root")
    run.add_argument("--skip-search", action="store_true")
    run.add_argument("--input-tsv")
    run.add_argument("--legacy-output", action="store_true")
    run.add_argument("--query-structure")
    run.add_argument("--query-chain")
    run.add_argument("--parent-prediction-id")

    benchmark = subparsers.add_parser("benchmark", help="Run the fixed Bur01-like acceptance benchmark")
    benchmark.add_argument("--engine", choices=("foldseek", "tmtools"), default="foldseek")
    benchmark.add_argument("--output-root")
    benchmark.add_argument("--no-prepare", action="store_true")
    benchmark.add_argument("--query-structure")
    benchmark.add_argument("--query-chain")
    benchmark.add_argument("--parent-prediction-id")

    production = subparsers.add_parser("production-run", help="Run the full AlphaFold/Swiss-Prot Foldseek workflow")
    production.add_argument("--config", default="configs/alphafold_swissprot.yaml")
    production.add_argument("--run-id")
    production.add_argument("--output-root")
    production.add_argument("--query-structure")
    production.add_argument("--query-chain")
    production.add_argument("--parent-prediction-id")

    predict = subparsers.add_parser("predict", help="Predict one monomer with local AlphaFold2-PTM")
    _add_prediction_arguments(predict)

    predict_and_mine = subparsers.add_parser(
        "predict-and-mine",
        help="Predict rank_001 locally, then use it as the structure-mining query",
    )
    _add_prediction_arguments(predict_and_mine)
    predict_and_mine.add_argument("--mining-config", default="configs/alphafold_swissprot.yaml")
    predict_and_mine.add_argument("--mining-output-root")

    ui = subparsers.add_parser("ui", help="Start the bilingual local web workbench")
    ui.add_argument("--port", type=int)
    ui.add_argument("--no-browser", action="store_true")

    validate = subparsers.add_parser("validate", help="Run isolated automated tests")
    validate.add_argument("pytest_args", nargs="*")
    validate_live = subparsers.add_parser("validate-live", help="Run tests plus live environment diagnostics")
    validate_live.add_argument("--config", default="config.yaml")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "check-files":
        return check_files()
    if args.command == "list-sources":
        configs = args.config or [
            "config.yaml",
            "configs/benchmark_foldseek.yaml",
            "configs/alphafold_swissprot.yaml",
        ]
        return list_sources(configs)
    if args.command == "doctor":
        report = run_doctor(args.config, check_network=not args.no_network)
        print(report.to_json() if args.as_json else "", end="" if args.as_json else "")
        if not args.as_json:
            print_doctor(report)
        elif report.to_json():
            print()
        return 1 if args.strict and not report.passed else 0
    if args.command == "run":
        request = RunRequest(
            config_path=args.config,
            run_id=args.run_id,
            output_root=args.output_root,
            skip_search=args.skip_search,
            input_tsv=args.input_tsv,
            query_structure=args.query_structure,
            query_chain=args.query_chain,
            parent_prediction_id=args.parent_prediction_id,
        )
        summary = run_pipeline(
            config_path=request.config_path,
            run_id=request.run_id,
            output_root=request.output_root,
            skip_foldseek=request.skip_search,
            input_tsv=request.input_tsv,
            legacy_output=args.legacy_output,
            query_structure=request.query_structure,
            query_chain=request.query_chain,
            parent_prediction_id=request.parent_prediction_id,
        )
        print(summary.model_dump_json(indent=2))
        return 0
    if args.command == "benchmark":
        _, report, report_path = run_fixed_benchmark(
            engine=args.engine,
            output_root=args.output_root,
            prepare=not args.no_prepare,
            query_structure=args.query_structure,
            query_chain=args.query_chain,
            parent_prediction_id=args.parent_prediction_id,
        )
        print(json.dumps(report, indent=2, ensure_ascii=False))
        print(f"Benchmark report: {report_path}")
        return 0 if report["passed"] else 1
    if args.command == "production-run":
        summary = run_pipeline(
            args.config,
            run_id=args.run_id,
            output_root=args.output_root,
            query_structure=args.query_structure,
            query_chain=args.query_chain,
            parent_prediction_id=args.parent_prediction_id,
        )
        print(summary.model_dump_json(indent=2))
        return 0
    if args.command == "predict":
        summary = run_prediction(_prediction_request(args))
        print(summary.model_dump_json(indent=2))
        return 0
    if args.command == "predict-and-mine":
        prediction, mining = run_prediction_and_mine(
            _prediction_request(args),
            mining_config=args.mining_config,
            mining_output_root=args.mining_output_root,
        )
        print(
            json.dumps(
                {
                    "prediction": json.loads(prediction.model_dump_json()),
                    "mining": json.loads(mining.model_dump_json()),
                },
                indent=2,
                ensure_ascii=False,
            )
        )
        return 0
    if args.command == "ui":
        return launch_ui(args.port, args.no_browser)
    if args.command == "validate":
        return _run_tests(args.pytest_args)
    if args.command == "validate-live":
        test_status = _run_tests()
        report = run_doctor(args.config, check_network=True)
        print_doctor(report)
        return test_status or (0 if report.passed else 1)
    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
