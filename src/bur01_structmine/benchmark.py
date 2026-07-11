"""Fixed Bur01-like benchmark preparation, execution, and acceptance checks."""
from __future__ import annotations

import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from .alphafold_client import AlphaFoldClient
from .config import load_config, project_root, resolve_path
from .datasources import load_search_databases
from .foldseek_runner import create_foldseek_database
from .http_utils import CachedSession
from .pipeline import LOCAL_ENGINES, run_pipeline
from .run_context import file_sha256
from .structure_search import load_ca_traces

SELF_LABEL = "8HFW"
POSITIVE_LABELS = ("8HFW", "1XKL", "6TJ2", "Q6RYA0", "Q9FVW3")
NEGATIVE_LABELS = ("P60174", "P00520")
AF_ACCESSIONS = ("Q6RYA0", "Q9FVW3", "P60174", "P00520")
PDB_IDS = ("1XKL", "6TJ2")


def _atomic_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(content, encoding="utf-8")
    tmp.replace(path)


def _validate_structure_file(path: Path) -> dict[str, Any]:
    traces = load_ca_traces(path)
    return {
        "chains": [trace.chain_id for trace in traces],
        "residues_by_chain": {trace.chain_id: len(trace.residues) for trace in traces},
    }


def prepare_benchmark_structures(force: bool = False) -> Path:
    """Download the fixed controls and write a checksummed provenance manifest."""
    root = project_root()
    structure_dir = root / "data" / "benchmark" / "structures"
    structure_dir.mkdir(parents=True, exist_ok=True)
    http = CachedSession(root / ".cache" / "downloads" / "benchmark", timeout=60)
    records: list[dict[str, Any]] = []

    query_source = root / "data" / "query" / "Bur01m.pdb"
    self_path = structure_dir / "8HFW_A.pdb"
    if force or not self_path.exists() or file_sha256(self_path) != file_sha256(query_source):
        shutil.copy2(query_source, self_path)
    self_info = _validate_structure_file(self_path)
    records.append(
        {
            "label": SELF_LABEL,
            "role": "self_positive",
            "path": str(self_path.relative_to(root)),
            "source_url": "https://www.rcsb.org/structure/8HFW",
            "sha256": file_sha256(self_path),
            **self_info,
        }
    )

    for pdb_id in PDB_IDS:
        path = structure_dir / f"{pdb_id}.pdb"
        source_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        if force or not path.exists():
            text = http.get_text(source_url)
            if not text or "ATOM" not in text:
                raise RuntimeError(f"RCSB download failed for {pdb_id}: {http.last_error or 'no ATOM records'}")
            _atomic_text(path, text)
        records.append(
            {
                "label": pdb_id,
                "role": "related_positive",
                "path": str(path.relative_to(root)),
                "source_url": source_url,
                "sha256": file_sha256(path),
                **_validate_structure_file(path),
            }
        )

    af_client = AlphaFoldClient(0.15, cache_dir=root / ".cache" / "annotations" / "alphafold")
    af_records = af_client.fetch_many(AF_ACCESSIONS)
    for accession in AF_ACCESSIONS:
        model = af_records.get(accession)
        if model is None or not model.pdb_url:
            raise RuntimeError(f"AlphaFold model metadata failed for {accession}: {af_client.errors.get(accession, 'no PDB URL')}")
        version = model.latest_version or "latest"
        path = structure_dir / f"AF-{accession}-F1-model_v{version}.pdb"
        if force or not path.exists():
            text = http.get_text(model.pdb_url)
            if not text or "ATOM" not in text:
                raise RuntimeError(f"AlphaFold structure download failed for {accession}: {http.last_error or 'no ATOM records'}")
            _atomic_text(path, text)
        records.append(
            {
                "label": accession,
                "role": "related_positive" if accession in POSITIVE_LABELS else ("length_matched_negative" if accession == "P60174" else "remote_negative"),
                "path": str(path.relative_to(root)),
                "source_url": model.pdb_url,
                "model_version": version,
                "global_plddt": model.global_plddt,
                "sha256": file_sha256(path),
                **_validate_structure_file(path),
            }
        )

    manifest = {
        "benchmark_id": "bur01_like_v0.1",
        "prepared_at": datetime.now(timezone.utc).isoformat(),
        "records": records,
    }
    manifest_path = root / "data" / "benchmark" / "download_manifest.json"
    _atomic_text(manifest_path, json.dumps(manifest, indent=2, ensure_ascii=False))
    return manifest_path


def _label_for_target(target: str) -> str | None:
    upper = str(target).upper()
    for label in (*POSITIVE_LABELS, *NEGATIVE_LABELS):
        if label in upper:
            return label
    return None


def validate_benchmark_tables(
    raw_df: pd.DataFrame,
    filtered_df: pd.DataFrame,
    self_tm_threshold: float = 0.98,
    require_self_rmsd: bool = True,
) -> dict[str, Any]:
    """Apply the fixed scientific acceptance criteria to search outputs."""
    table = raw_df.copy()
    table["benchmark_label"] = table["target"].map(_label_for_target)
    table = table[table["benchmark_label"].notna()].copy()
    best = (
        table.sort_values("alntmscore", ascending=False)
        .drop_duplicates("benchmark_label")
        .set_index("benchmark_label")
    )
    missing = [label for label in (*POSITIVE_LABELS, *NEGATIVE_LABELS) if label not in best.index]
    checks: list[dict[str, Any]] = []

    if SELF_LABEL in best.index:
        self_row = best.loc[SELF_LABEL]
        self_tm = float(self_row.get("alntmscore", 0.0))
        self_rmsd = float(self_row.get("rmsd", 999.0))
        checks.append(
            {
                "name": "self_tm_score",
                "passed": self_tm >= self_tm_threshold,
                "observed": self_tm,
                "required": f">= {self_tm_threshold:.2f}",
            }
        )
        checks.append(
            {
                "name": "self_rmsd",
                "passed": self_rmsd < 0.1 if require_self_rmsd else True,
                "observed": self_rmsd,
                "required": "< 0.1 A" if require_self_rmsd else "reported, not identity-gated for predicted query",
            }
        )
    else:
        checks.extend(
            [
                {"name": "self_tm_score", "passed": False, "observed": None, "required": "> 0.98"},
                {"name": "self_rmsd", "passed": False, "observed": None, "required": "< 0.1 A"},
            ]
        )

    if not missing:
        positive_floor = min(float(best.loc[label, "alntmscore"]) for label in POSITIVE_LABELS)
        negative_ceiling = max(float(best.loc[label, "alntmscore"]) for label in NEGATIVE_LABELS)
        checks.append(
            {
                "name": "positive_ranking",
                "passed": positive_floor > negative_ceiling,
                "observed": {"lowest_positive": positive_floor, "highest_negative": negative_ceiling},
                "required": "every positive above both negatives",
            }
        )
    else:
        checks.append({"name": "positive_ranking", "passed": False, "observed": {"missing": missing}, "required": "all controls observed"})

    filtered_labels = {_label_for_target(target) for target in filtered_df.get("target", pd.Series(dtype=str))}
    p00520_filtered = "P00520" in filtered_labels
    checks.append(
        {
            "name": "p00520_default_filter",
            "passed": not p00520_filtered,
            "observed": "passed" if p00520_filtered else "rejected",
            "required": "rejected",
        }
    )
    ranked = [
        {
            "label": label,
            "target": str(best.loc[label, "target"]),
            "alntmscore": float(best.loc[label, "alntmscore"]),
            "qtmscore": float(best.loc[label, "qtmscore"]) if pd.notna(best.loc[label].get("qtmscore")) else None,
            "ttmscore": float(best.loc[label, "ttmscore"]) if pd.notna(best.loc[label].get("ttmscore")) else None,
            "rmsd": float(best.loc[label, "rmsd"]) if pd.notna(best.loc[label].get("rmsd")) else None,
        }
        for label in best.sort_values("alntmscore", ascending=False).index
    ]
    return {
        "benchmark_id": "bur01_like_v0.1",
        "passed": not missing and all(check["passed"] for check in checks),
        "missing_controls": missing,
        "checks": checks,
        "ranking": ranked,
    }


def run_fixed_benchmark(
    engine: str = "foldseek",
    output_root: str | Path | None = None,
    prepare: bool = True,
    query_structure: str | None = None,
    query_chain: str | None = None,
    parent_prediction_id: str | None = None,
    catalytic_reference_sites: dict[str, int] | None = None,
) -> tuple[Any, dict[str, Any], Path]:
    """Run the fixed benchmark through the shared pipeline service."""
    engine = engine.strip().lower()
    if engine not in {"foldseek", "tmtools"}:
        raise ValueError("Benchmark engine must be 'foldseek' or 'tmtools'")
    if prepare:
        prepare_benchmark_structures()
    config_path = f"configs/benchmark_{engine}.yaml"
    cfg = load_config(config_path)
    source = load_search_databases(cfg)[0]
    if source.search_engine not in LOCAL_ENGINES:
        if not source.structure_dir or not source.target_db:
            raise ValueError("Foldseek benchmark requires structure_dir and target_db")
        if not (source.target_db.exists() or Path(f"{source.target_db}.dbtype").exists()):
            fs_cfg = cfg.get("foldseek", {})
            create_foldseek_database(
                source.structure_dir,
                source.target_db,
                executable=str(fs_cfg.get("executable", "foldseek")),
                backend=str(fs_cfg.get("backend", "auto")),
                wsl_distribution=str(fs_cfg.get("wsl_distribution", "Ubuntu-Foldseek")),
                wsl_executable=str(fs_cfg.get("wsl_executable", "/opt/foldseek/bin/foldseek")),
            )
    summary = run_pipeline(
        config_path=config_path,
        output_root=str(output_root) if output_root else None,
        query_structure=query_structure,
        query_chain=query_chain,
        parent_prediction_id=parent_prediction_id,
        catalytic_reference_sites=catalytic_reference_sites,
    )
    raw_df = pd.read_csv(summary.artifacts["raw_tsv"], sep="\t")
    filtered_df = pd.read_csv(summary.artifacts["filtered_tsv"], sep="\t")
    report = validate_benchmark_tables(
        raw_df,
        filtered_df,
        self_tm_threshold=0.75 if query_structure else 0.98,
        require_self_rmsd=not bool(query_structure),
    )
    report["search_engine"] = engine
    report["query_mode"] = "predicted" if query_structure else "experimental_self"
    report["run_id"] = summary.run_id
    report["run_manifest"] = str(summary.manifest_path.resolve())
    report_path = summary.run_dir / "benchmark_report.json"
    _atomic_text(report_path, json.dumps(report, indent=2, ensure_ascii=False))
    return summary, report, report_path
