"""Copy-on-write WASE rescoring with source hashes and explicit evidence gaps."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .config import project_root
from .experiments import safe_table
from .run_context import file_sha256, make_run_id
from .results_parser import parse_pdb_chain
from .scoring import WASE_VERSION, add_candidate_scores
from .structure_search import compare_ca_traces, load_ca_trace, load_ca_traces


def local_candidate_path(row: pd.Series, root: Path) -> Path | None:
    path_text = str(row.get("target_structure_path", ""))
    if path_text and path_text != "nan" and Path(path_text).is_file():
        return Path(path_text)
    target = str(row.get("target", ""))
    candidates = [root / ".cache" / "ui_structures" / f"{target}.pdb"]
    for path in (root / "data" / "benchmark" / "structures").glob("*.pdb"):
        if any(label in target.upper() and label in path.name.upper() for label in ("8HFW", "1XKL", "6TJ2", "Q6RYA0", "Q9FVW3", "P60174", "P00520")):
            candidates.append(path)
    return next((path for path in candidates if path.is_file()), None)


def enrich_local_catalytic_evidence(frame: pd.DataFrame, query_path: Path, query_chain: str,
                                   sites: dict[str, int], root: Path | None = None) -> pd.DataFrame:
    root = root or project_root()
    out = frame.copy()
    for column in out:
        if column.startswith(("mapped_catalytic", "catalytic_", "query_catalytic")):
            out[column] = out[column].astype(object)
    query = load_ca_trace(query_path, chain_id=query_chain)
    reference = tuple((number, label[0]) for label, number in sites.items())
    for index, row in out.iterrows():
        path = local_candidate_path(row, root)
        if not path:
            out.at[index, "catalytic_evidence_status"] = "local_coordinates_unavailable"
            continue
        if not reference:
            out.at[index, "catalytic_evidence_status"] = "reference_sites_not_mapped"
            continue
        try:
            chain = str(row.get("target_chain", row.get("chain", "")))
            if chain in {"", "nan", "None"}:
                _, chain = parse_pdb_chain(str(row.get("target", ""))) if not str(row.get("target", "")).startswith("AF-") else (None, None)
            target = load_ca_trace(path, chain_id=chain)
            evidence = compare_ca_traces(query, target, reference)
            # Preserve Foldseek metrics; TM-align is only the provenance for this mapping.
            for field, value in evidence.items():
                if field.startswith(("mapped_catalytic", "catalytic_", "query_catalytic")):
                    if field not in out:
                        out[field] = pd.Series(index=out.index, dtype=object)
                    out.at[index, field] = value
            out.at[index, "catalytic_evidence_status"] = "local_coordinates_assessed"
            out.at[index, "catalytic_mapping_source"] = "tmtools-0.3.0/TM-align-20210224"
            out.at[index, "catalytic_structure_sha256"] = file_sha256(path)
            out.at[index, "target_chain"] = target.chain_id
            out.at[index, "query_chain"] = query.chain_id
        except (ValueError, OSError) as exc:
            out.at[index, "catalytic_evidence_status"] = f"not_assessable:{exc}"
    return out


def rescore_run(run_dir: str | Path, output_root: str | Path | None = None) -> Path:
    source = Path(run_dir).resolve()
    manifest_path = source / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    table_path = Path(manifest["artifacts"]["annotated_tsv"])
    if not table_path.is_file():
        table_path = source / table_path.name
    original_hash = file_sha256(table_path)
    frame = pd.read_csv(table_path, sep="\t")
    sites = manifest.get("catalytic_reference_sites")
    if not sites and not manifest.get("parent_prediction_id"):
        sites = {"S111": 111, "D241": 241, "H274": 274} if manifest.get("query_pdb_id", "").upper() == "8HFW" else {}
    query = Path(manifest["query_path"])
    enriched = enrich_local_catalytic_evidence(frame, query, manifest.get("query_chain") or "A", sites) if query.is_file() else frame.copy()
    output = add_candidate_scores(enriched)
    output["source_run_id"] = manifest["run_id"]
    output["previous_wase_version"] = frame.get("wase_version", "unknown")
    output["previous_wase_screen_score"] = frame.get("wase_screen_score", float("nan"))
    output["wase_score_delta"] = output.wase_screen_score - pd.to_numeric(output.previous_wase_screen_score, errors="coerce")
    destination = Path(output_root or project_root() / "results" / "rescores") / make_run_id("rescore")
    destination.mkdir(parents=True, exist_ok=False)
    safe_table(output).to_csv(destination / "rescored.tsv", sep="\t", index=False)
    safe_table(output).to_excel(destination / "rescored.xlsx", index=False)
    metadata = {"rescore_id": destination.name, "source_run_id": manifest["run_id"], "source_manifest_sha256": file_sha256(manifest_path),
                "source_table_sha256": original_hash, "wase_version": WASE_VERSION, "rows": len(output),
                "coordinate_policy": "local_only_no_download", "artifact": "rescored.tsv", "source_preserved": file_sha256(table_path) == original_hash}
    (destination / "rescore_manifest.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return destination
