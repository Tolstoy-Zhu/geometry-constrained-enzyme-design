"""Validate the real v0.3 AF2 predictions and their linked mining runs."""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from bur01_structmine.config import project_root
from bur01_structmine.models import PredictionManifest, RunManifest
from bur01_structmine.structure_search import compare_ca_traces, load_ca_trace


def _check(name: str, passed: bool, observed: Any, required: str) -> dict[str, Any]:
    return {"name": name, "passed": bool(passed), "observed": observed, "required": required}


def _prediction_checks(
    label: str,
    manifest: PredictionManifest,
    reference_path: Path,
    require_candidate_lineage: bool,
) -> list[dict[str, Any]]:
    checks: list[dict[str, Any]] = []
    rank_path = Path(manifest.rank_001_path)
    query = load_ca_trace(rank_path, chain_id="A")
    reference = load_ca_trace(reference_path, chain_id="A")
    comparison = compare_ca_traces(query, reference, catalytic_query_residues=())
    minimum_tm = min(float(comparison["qtmscore"]), float(comparison["ttmscore"]))
    checks.extend(
        [
            _check(f"{label}_status", manifest.status == "completed", manifest.status, "completed"),
            _check(f"{label}_model_count", len(manifest.models) == 5, len(manifest.models), "5"),
            _check(
                f"{label}_mean_plddt",
                manifest.mean_plddt is not None and manifest.mean_plddt >= 70.0,
                manifest.mean_plddt,
                ">= 70",
            ),
            _check(
                f"{label}_ptm",
                manifest.ptm is not None and 0.0 < manifest.ptm <= 1.0,
                manifest.ptm,
                "0 < pTM <= 1",
            ),
            _check(f"{label}_minimum_tm", minimum_tm >= 0.75, minimum_tm, ">= 0.75"),
            _check(
                f"{label}_parameter_hashes",
                bool(manifest.model_parameter_checksums),
                len(manifest.model_parameter_checksums),
                "non-empty",
            ),
            _check(
                f"{label}_parameter_archive_hash",
                bool(manifest.model_parameter_provenance.get("archive_sha256")),
                manifest.model_parameter_provenance.get("archive_sha256"),
                "recorded SHA-256",
            ),
        ]
    )
    if require_candidate_lineage:
        checks.extend(
            [
                _check(
                    f"{label}_source_run",
                    bool(manifest.source_run_id),
                    manifest.source_run_id,
                    "traceable run ID",
                ),
                _check(
                    f"{label}_source_accession",
                    manifest.source_accession.upper() == "Q9FVW3",
                    manifest.source_accession,
                    "Q9FVW3",
                ),
            ]
        )
    return checks


def _linked_run_manifests(manifest: PredictionManifest) -> list[RunManifest]:
    root = project_root() / "results" / "runs"
    linked: list[RunManifest] = []
    for run_id in manifest.child_run_ids:
        path = root / run_id / "run_manifest.json"
        if path.is_file():
            linked.append(RunManifest.model_validate_json(path.read_text(encoding="utf-8")))
    return linked


def _is_swissprot_run(manifest: RunManifest) -> bool:
    sources = manifest.config_snapshot.get("search_databases") or []
    if isinstance(sources, dict):
        sources = list(sources.values())
    return any(source.get("name") == "alphafold_swissprot" for source in sources if isinstance(source, dict))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prediction-8hfw", required=True)
    parser.add_argument("--prediction-q9fvw3", required=True)
    parser.add_argument("--output", default="results/validation/af2_acceptance_report.json")
    args = parser.parse_args()

    root = project_root()
    prediction_root = root / "results" / "predictions"
    eight_path = prediction_root / args.prediction_8hfw / "prediction_manifest.json"
    q9_path = prediction_root / args.prediction_q9fvw3 / "prediction_manifest.json"
    eight = PredictionManifest.model_validate_json(eight_path.read_text(encoding="utf-8"))
    q9 = PredictionManifest.model_validate_json(q9_path.read_text(encoding="utf-8"))
    checks = _prediction_checks(
        "8hfw",
        eight,
        root / "data" / "benchmark" / "structures" / "8HFW_A.pdb",
        require_candidate_lineage=False,
    )
    checks.extend(
        _prediction_checks(
            "q9fvw3",
            q9,
            root / "data" / "benchmark" / "structures" / "AF-Q9FVW3-F1-model_v6.pdb",
            require_candidate_lineage=True,
        )
    )

    linked = _linked_run_manifests(eight)
    benchmark_reports = []
    for child in linked:
        report_path = root / "results" / "runs" / child.run_id / "benchmark_report.json"
        if report_path.is_file():
            benchmark_reports.append(json.loads(report_path.read_text(encoding="utf-8")))
    passing_predicted_benchmark = next(
        (
            report
            for report in benchmark_reports
            if report.get("passed") and report.get("query_mode") == "predicted"
        ),
        None,
    )
    checks.append(
        _check(
            "predicted_8hfw_benchmark",
            passing_predicted_benchmark is not None,
            passing_predicted_benchmark.get("run_id") if passing_predicted_benchmark else None,
            "passing Foldseek predicted-query benchmark",
        )
    )

    production = next(
        (
            manifest
            for manifest in linked
            if manifest.status == "completed"
            and _is_swissprot_run(manifest)
            and manifest.raw_hit_count > 0
            and manifest.annotated_hit_count > 0
            and not manifest.error
        ),
        None,
    )
    checks.append(
        _check(
            "predicted_8hfw_swissprot_run",
            production is not None,
            {
                "run_id": production.run_id,
                "raw": production.raw_hit_count,
                "filtered": production.filtered_hit_count,
                "annotated": production.annotated_hit_count,
            }
            if production
            else None,
            "completed, non-empty, annotated, no error",
        )
    )

    report = {
        "project_version": "0.3.0",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "prediction_8hfw": eight.prediction_id,
        "prediction_q9fvw3": q9.prediction_id,
        "passed": all(check["passed"] for check in checks),
        "checks": checks,
    }
    output = Path(args.output)
    if not output.is_absolute():
        output = root / output
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = output.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    tmp.replace(output)
    print(json.dumps(report, indent=2, ensure_ascii=False))
    print(f"Report: {output}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
