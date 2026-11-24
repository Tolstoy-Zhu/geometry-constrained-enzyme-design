"""Build the public read-only bundle only from a passing real Foldseek run."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from bur01_structmine.benchmark import _label_for_target, validate_benchmark_tables
from bur01_structmine.config import project_root
from bur01_structmine.experiments import simulated_experiments, template_bytes
from bur01_structmine.rescoring import enrich_local_catalytic_evidence, local_candidate_path
from bur01_structmine.run_context import file_sha256
from bur01_structmine.scoring import add_candidate_scores


def build(run_dir: Path) -> Path:
    root = project_root()
    manifest_path = run_dir / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest["status"] != "completed" or "foldseek" not in manifest["software_versions"]:
        raise ValueError("A completed real Foldseek run is required")
    raw = pd.read_csv(manifest["artifacts"]["raw_tsv"], sep="\t")
    filtered = pd.read_csv(manifest["artifacts"]["filtered_tsv"], sep="\t")
    report = validate_benchmark_tables(raw, filtered)
    if not report["passed"]:
        raise ValueError("Fixed benchmark acceptance failed")
    raw["benchmark_label"] = raw.target.map(_label_for_target)
    raw = raw.sort_values("alntmscore", ascending=False).drop_duplicates("benchmark_label").copy()
    annotations = pd.read_csv(manifest["artifacts"]["annotated_tsv"], sep="\t").set_index("target")
    for col in annotations:
        if col not in raw:
            raw[col] = raw.target.map(annotations[col])
    raw["default_filter_pass"] = raw.target.isin(filtered.target)
    thresholds = manifest["config_snapshot"]["filter"]
    def exclusions(row):
        failed = []
        for field, threshold_key, minimum in [("alntmscore", "min_alntmscore", True), ("rmsd", "max_rmsd", False),
                                             ("qcov", "min_qcov", True), ("tcov", "min_tcov", True),
                                             ("alnlen", "min_alnlen", True), ("fident", "max_sequence_identity", False)]:
            if (row[field] < thresholds[threshold_key]) if minimum else (row[field] > thresholds[threshold_key]):
                failed.append(f"{field}{'<' if minimum else '>'}{thresholds[threshold_key]}")
        return ";".join(failed)
    raw["filter_exclusion_reasons"] = raw.apply(exclusions, axis=1)
    query = root / "data" / "query" / "Bur01m.pdb"
    raw = enrich_local_catalytic_evidence(raw, query, "A", {"S111": 111, "D241": 241, "H274": 274})
    raw["demo_structure_path"] = raw.apply(lambda row: str(local_candidate_path(row, root).relative_to(root)).replace("\\", "/"), axis=1)
    table = add_candidate_scores(raw)
    table["source_run_id"] = manifest["run_id"]
    # Never bundle local paths, production candidates, logs or remote credentials.
    table = table.drop(columns=[col for col in table if col.endswith("_path") and col != "demo_structure_path"] + ["target_database"], errors="ignore")
    destination = root / "data" / "demo"
    destination.mkdir(parents=True, exist_ok=True)
    if (destination / "demo_manifest.json").exists():
        raise FileExistsError("Demo bundle already exists; preserve it before rebuilding")
    table.to_csv(destination / "candidates.tsv", sep="\t", index=False)
    examples = root / "data" / "examples" / "simulated"
    examples.mkdir(parents=True, exist_ok=True)
    simulated_experiments().to_csv(examples / "experimental_measurements.csv", index=False, encoding="utf-8-sig")
    templates = root / "data" / "templates"
    templates.mkdir(parents=True, exist_ok=True)
    (templates / "experiment_template.xlsx").write_bytes(template_bytes())
    files = [query, destination / "candidates.tsv", root / "ui" / "assets" / "3Dmol-2.5.5-min.js"]
    files += [root / relative for relative in table.demo_structure_path]
    provenance = {"demo_version": "0.5.0.dev0", "source_run_id": manifest["run_id"],
                  "source_manifest_sha256": file_sha256(manifest_path), "source_raw_sha256": file_sha256(manifest["artifacts"]["raw_tsv"]),
                  "software_versions": manifest["software_versions"], "benchmark_checks": report["checks"],
                  "data_kind": "real_computational_records", "experimental_example": "simulated_separate_not_Bur01_results",
                  "files": {str(path.relative_to(root)).replace("\\", "/"): file_sha256(path) for path in files}}
    (destination / "demo_manifest.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    return destination


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    print(build(parser.parse_args().run_dir))
