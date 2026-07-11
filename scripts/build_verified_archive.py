"""Build a verified v0.2 archive only after formal acceptance evidence exists."""
from __future__ import annotations

import csv
import json
import re
import zipfile
from datetime import datetime
from pathlib import Path

from bur01_structmine.config import project_root
from bur01_structmine.run_context import file_sha256

EXCLUDED_PARTS = {
    ".cache",
    ".pytest_cache",
    "__pycache__",
    "tmp",
    "databases",
    "results",
}


def _latest_formal_benchmark(run_root: Path) -> tuple[Path, dict]:
    accepted = []
    for report_path in run_root.glob("*/benchmark_report.json"):
        report = json.loads(report_path.read_text(encoding="utf-8"))
        if report.get("passed") and report.get("search_engine") == "foldseek":
            accepted.append((report_path, report))
    if not accepted:
        raise RuntimeError("No passing formal Foldseek benchmark report was found")
    return max(accepted, key=lambda item: item[0].stat().st_mtime)


def _is_production_manifest(manifest: dict) -> bool:
    entries = manifest.get("config_snapshot", {}).get("search_databases") or []
    if isinstance(entries, dict):
        entries = list(entries.values())
    names = {entry.get("name") for entry in entries if isinstance(entry, dict)}
    return "alphafold_swissprot" in names


def _latest_valid_production(run_root: Path) -> tuple[Path, dict]:
    accepted = []
    for manifest_path in run_root.glob("*/run_manifest.json"):
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("status") != "completed" or not _is_production_manifest(manifest):
            continue
        if int(manifest.get("raw_hit_count", 0)) <= 0 or int(manifest.get("annotated_hit_count", 0)) <= 0:
            continue
        annotated = Path(manifest.get("artifacts", {}).get("annotated_tsv", ""))
        if not annotated.is_file():
            continue
        with annotated.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        valid_rows = [
            row
            for row in rows
            if row.get("final_protein_sequence", "").strip()
            and row.get("annotation_status", "").strip() in {"complete", "partial"}
        ]
        if valid_rows:
            accepted.append((manifest_path, manifest))
    if not accepted:
        raise RuntimeError("No completed AlphaFold/Swiss-Prot production run with valid annotations was found")
    return max(accepted, key=lambda item: item[0].stat().st_mtime)


def _source_files(root: Path):
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(root)
        if any(part in EXCLUDED_PARTS or part.endswith(".egg-info") for part in relative.parts):
            continue
        yield path, relative


def _unique_archive_path(parent: Path) -> Path:
    stem = f"AquaEsterNet-StructMine_{datetime.now().strftime('%Y%m%d')}_v0.2_verified"
    candidate = parent / f"{stem}.zip"
    if not candidate.exists():
        return candidate
    return parent / f"{stem}_{datetime.now().strftime('%H%M%S')}.zip"


def main() -> int:
    root = project_root()
    run_root = root / "results" / "runs"
    benchmark_path, benchmark = _latest_formal_benchmark(run_root)
    production_path, production = _latest_valid_production(run_root)
    audit_text = (root / "COMPLETION_AUDIT.md").read_text(encoding="utf-8")
    if "| PENDING |" in audit_text or "not marked complete" in audit_text:
        raise RuntimeError("Completion audit still contains pending requirements")
    validation_text = (root / "VALIDATION_REPORT.md").read_text(encoding="utf-8")
    for required_run_id in (benchmark.get("run_id"), production.get("run_id")):
        if required_run_id and required_run_id not in validation_text:
            raise RuntimeError(f"Validation report does not cite required run: {required_run_id}")
    environment_dir = root / "databases" / "environment"
    required_environment = [
        environment_dir / "foldseek_install_manifest.json",
        environment_dir / "afdb_swissprot_download_provenance.json",
        environment_dir / "alphafold_swissprot_database_manifest.json",
    ]
    missing_environment = [str(path) for path in required_environment if not path.is_file()]
    if missing_environment:
        raise RuntimeError("Missing environment manifest(s): " + ", ".join(missing_environment))

    archive = _unique_archive_path(root.parent)
    evidence_paths = [
        benchmark_path,
        benchmark_path.parent / "run_manifest.json",
        production_path,
        *[Path(value) for value in production.get("artifacts", {}).values()],
        *required_environment,
    ]
    screenshots = [path for path in (root / "results").glob("ui_*.png") if path.is_file()]
    evidence_paths.extend(screenshots)
    evidence_paths = list(dict.fromkeys(path.resolve() for path in evidence_paths if path.is_file()))

    verification_index = {
        "project_version": "0.2.0",
        "created_at": datetime.now().astimezone().isoformat(),
        "formal_benchmark_run": benchmark.get("run_id"),
        "production_run": production.get("run_id"),
        "files": [
            {"path": str(path), "sha256": file_sha256(path), "size_bytes": path.stat().st_size}
            for path in evidence_paths
        ],
    }

    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as bundle:
        for path, relative in _source_files(root):
            bundle.write(path, Path("AquaEsterNet-StructMine") / relative)
        for path in evidence_paths:
            if path.is_relative_to(root):
                relative = path.relative_to(root)
            else:
                relative = Path("external_evidence") / re.sub(r"[^A-Za-z0-9_.-]+", "_", path.name)
            bundle.write(path, Path("AquaEsterNet-StructMine") / relative)
        bundle.writestr(
            "AquaEsterNet-StructMine/VERIFICATION_INDEX.json",
            json.dumps(verification_index, indent=2, ensure_ascii=False),
        )
    print(f"Archive: {archive}")
    print(f"SHA-256: {file_sha256(archive)}")
    print(f"Formal benchmark: {benchmark.get('run_id')}")
    print(f"Production run: {production.get('run_id')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
