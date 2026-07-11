from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import yaml

from bur01_structmine.config import project_root
from bur01_structmine.pipeline import run_pipeline


def test_pipeline_writes_versioned_run_to_temporary_directory(tmp_path: Path) -> None:
    root = project_root()
    config = {
        "query": {
            "pdb_file": str(root / "data" / "query" / "Bur01m.pdb"),
            "fasta_file": str(root / "data" / "query" / "Bur01m_seed.fasta"),
            "pdb_id": "8HFW",
            "chain": "A",
            "residue_range": "12-294",
        },
        "foldseek": {"prepare_databases": False},
        "search_databases": [{"name": "fixture", "enabled": True, "result_tsv": str(root / "data" / "example_foldseek_result.tsv")}],
        "filter": {
            "min_alntmscore": 0.45,
            "max_rmsd": 4,
            "min_qcov": 0.6,
            "min_tcov": 0.4,
            "min_alnlen": 180,
            "max_sequence_identity": 1.0,
        },
        "annotation": {
            "enable_rcsb": False,
            "enable_uniprot": False,
            "enable_ncbi": False,
            "enable_alphafold": False,
            "cache_dir": str(tmp_path / "cache"),
            "max_records": 20,
        },
        "output": {"out_prefix": "results/test_run", "run_root": str(tmp_path / "runs")},
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    summary = run_pipeline(
        config_path=str(config_path),
        skip_foldseek=True,
        input_tsv=str(root / "data" / "example_foldseek_result.tsv"),
        output_root=str(tmp_path / "runs"),
        query_structure=str(root / "data" / "query" / "Bur01m.pdb"),
        query_chain="A",
        parent_prediction_id="prediction-test",
        catalytic_reference_sites={"S111": 111, "D241": 241, "H274": 274},
    )
    assert summary.status == "completed"
    assert summary.raw_hit_count == 6
    assert summary.filtered_hit_count == 5
    assert summary.run_dir.parent == tmp_path / "runs"
    manifest = json.loads(summary.manifest_path.read_text(encoding="utf-8"))
    assert manifest["query_sha256"]
    assert manifest["query_chain"] == "A"
    assert manifest["parent_prediction_id"] == "prediction-test"
    assert manifest["catalytic_reference_sites"] == {"S111": 111, "D241": 241, "H274": 274}
    assert Path(manifest["artifacts"]["events"]).exists()
    annotated = pd.read_csv(summary.artifacts["annotated_tsv"], sep="\t")
    assert set(annotated["annotation_status"]) == {"disabled"}
    assert "wase_screen_uncertainty" in annotated.columns


def test_formal_results_are_not_touched_by_tests(tmp_path: Path) -> None:
    assert not (tmp_path / "results").exists()
