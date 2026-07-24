from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest
import pandas as pd
import yaml

from bur01_structmine.colabfold_runtime import (
    ColabFoldRuntime,
    PredictionCancelled,
    ProcessFileLock,
    ResourceBusyError,
)
from bur01_structmine.models import NrRunManifest, PredictionManifest, PredictionRequest
from bur01_structmine.config import project_root
from bur01_structmine.prediction import PrivacyConsentRequired, run_prediction, run_prediction_and_mine


FAKE_COLABFOLD = Path(__file__).parent / "fixtures" / "fake_colabfold.py"
SEQUENCE = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY"


def _request(tmp_path: Path, prediction_id: str, consent: bool = True, preset: str = "quick") -> PredictionRequest:
    fasta = tmp_path / "query.fasta"
    fasta.write_text(f">test_query\n{SEQUENCE}\n", encoding="ascii")
    return PredictionRequest(
        fasta_path=str(fasta),
        preset=preset,
        confirm_public_msa=consent,
        prediction_id=prediction_id,
        output_root=str(tmp_path / "predictions"),
        msa_cache_dir=str(tmp_path / "msa_cache"),
        runtime_backend="native",
        colabfold_executable=str(FAKE_COLABFOLD),
    )


@pytest.fixture(autouse=True)
def fake_gpu_environment(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("AQUA_SKIP_GPU_CHECK", "1")
    monkeypatch.setenv("AQUA_FAKE_COLABFOLD_VERSION", "1.6.1")
    monkeypatch.setenv("AQUA_MSA_RETRY_BASE_SECONDS", "0")


def test_fake_colabfold_end_to_end_and_msa_cache(tmp_path: Path) -> None:
    first = run_prediction(_request(tmp_path, "fake-first", consent=True, preset="standard"))
    assert first.status == "completed"
    assert first.model_count == 5
    assert first.rank_001_path and first.rank_001_path.is_file()
    manifest = PredictionManifest.model_validate_json(first.manifest_path.read_text(encoding="utf-8"))
    assert manifest.msa_cache_hit is False
    assert manifest.models[0].relaxed is True
    assert manifest.model_parameter_provenance["archive_sha256"] == "test-only"
    assert manifest.model_parameter_checksums == {"fake_parameters": "test-only"}
    assert Path(manifest.artifacts["pae"]).is_file()
    assert Path(manifest.artifacts["plddt"]).is_file()

    cached = run_prediction(_request(tmp_path, "fake-cached", consent=False, preset="quick"))
    cached_manifest = PredictionManifest.model_validate_json(cached.manifest_path.read_text(encoding="utf-8"))
    assert cached.status == "completed"
    assert cached_manifest.msa_cache_hit is True


def test_public_msa_requires_explicit_consent_on_cache_miss(tmp_path: Path) -> None:
    with pytest.raises(PrivacyConsentRequired):
        run_prediction(_request(tmp_path, "no-consent", consent=False))
    manifest_path = tmp_path / "predictions" / "no-consent" / "prediction_manifest.json"
    manifest = PredictionManifest.model_validate_json(manifest_path.read_text(encoding="utf-8"))
    assert manifest.status == "failed"
    assert "PrivacyConsentRequired" in manifest.error


def test_public_msa_retries_then_succeeds(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    counter = tmp_path / "msa_attempts.txt"
    monkeypatch.setenv("FAKE_COLABFOLD_FAIL_MSA_ATTEMPTS", "2")
    monkeypatch.setenv("FAKE_COLABFOLD_COUNTER", str(counter))
    summary = run_prediction(_request(tmp_path, "retry-success"))
    assert summary.status == "completed"
    assert counter.read_text(encoding="ascii") == "3"


def test_process_lock_recovers_stale_owner_and_rejects_live_owner(tmp_path: Path) -> None:
    path = tmp_path / "gpu.lock"
    path.write_text(json.dumps({"pid": 99999999}), encoding="utf-8")
    with ProcessFileLock(path, "GPU"):
        assert json.loads(path.read_text(encoding="utf-8"))["pid"] == os.getpid()
        with pytest.raises(ResourceBusyError):
            ProcessFileLock(path, "GPU").acquire()
    assert not path.exists()


def test_runtime_honors_existing_cancel_marker(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    fasta = tmp_path / "query.fasta"
    fasta.write_text(f">cancel\n{SEQUENCE}\n", encoding="ascii")
    output = tmp_path / "output"
    cancel = tmp_path / "cancel.requested"
    cancel.write_text("cancel", encoding="ascii")
    monkeypatch.setenv("FAKE_COLABFOLD_SLEEP", "3")
    runtime = ColabFoldRuntime(backend="native", executable=str(FAKE_COLABFOLD))
    with pytest.raises(PredictionCancelled):
        runtime.run(
            [str(fasta), str(output), "--msa-only"],
            tmp_path / "run.log",
            cancel,
            tmp_path / "pid",
            timeout=20,
        )


def test_fake_prediction_and_local_mining_are_bidirectionally_linked(tmp_path: Path) -> None:
    config = {
        "query": {"pdb_file": str(project_root() / "data" / "query" / "Bur01m.pdb"), "chain": "A"},
        "foldseek": {"prepare_databases": False},
        "search_databases": [
            {
                "name": "benchmark_structures",
                "enabled": True,
                "kind": "custom_structure",
                "structure_dir": str(project_root() / "data" / "benchmark" / "structures"),
                "search_engine": "tmtools_tmalign",
                "annotation_source": "none",
            }
        ],
        "filter": {
            "min_alntmscore": 0.0,
            "max_rmsd": 999.0,
            "min_qcov": 0.0,
            "min_tcov": 0.0,
            "min_alnlen": 1,
            "max_sequence_identity": 1.0,
        },
        "annotation": {
            "enable_rcsb": False,
            "enable_uniprot": False,
            "enable_ncbi": False,
            "enable_alphafold": False,
            "cache_dir": str(tmp_path / "annotation_cache"),
            "max_records": 50,
        },
        "output": {"out_prefix": "results/predict_and_mine", "run_root": str(tmp_path / "runs")},
    }
    config_path = tmp_path / "mining.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    prediction, mining = run_prediction_and_mine(
        _request(tmp_path, "linked-prediction"),
        mining_config=str(config_path),
        mining_output_root=str(tmp_path / "runs"),
    )
    assert prediction.child_run_ids == [mining.run_id]
    prediction_manifest = PredictionManifest.model_validate_json(prediction.manifest_path.read_text(encoding="utf-8"))
    assert mining.run_id in prediction_manifest.child_run_ids
    mining_manifest = json.loads(mining.manifest_path.read_text(encoding="utf-8"))
    assert mining_manifest["parent_prediction_id"] == prediction.prediction_id
    assert mining.raw_hit_count > 0


def test_nr_candidate_prediction_keeps_bidirectional_provenance(
    tmp_path: Path,
) -> None:
    nr_run_id = "nr-candidate-source"
    nr_run_dir = tmp_path / "nr_runs" / nr_run_id
    nr_run_dir.mkdir(parents=True)
    candidates_path = nr_run_dir / "candidates.tsv"
    pd.DataFrame(
        [
            {
                "target": "ref|WP_TEST_NR.1|",
                "target_accession": "WP_TEST_NR.1",
                "sequence_id": "nr_0000001",
                "target_sequence": SEQUENCE,
                "af2_ready": True,
                "wase_tier": "Not evaluated - structure required",
            }
        ]
    ).to_csv(candidates_path, sep="\t", index=False)
    nr_manifest_path = nr_run_dir / "nr_run_manifest.json"
    nr_manifest_path.write_text(
        NrRunManifest(
            nr_run_id=nr_run_id,
            status="completed",
            parent_structure_run_id="structure-run",
            database_scope="nr_remote_snapshot",
            artifacts={"candidates_tsv": str(candidates_path.resolve())},
        ).model_dump_json(indent=2),
        encoding="utf-8",
    )
    request = PredictionRequest(
        candidate_nr_run_id=nr_run_id,
        accession="WP_TEST_NR.1",
        preset="quick",
        confirm_public_msa=True,
        prediction_id="nr-linked-prediction",
        output_root=str(tmp_path / "predictions"),
        msa_cache_dir=str(tmp_path / "msa_cache"),
        candidate_nr_run_root=str(tmp_path / "nr_runs"),
        runtime_backend="native",
        colabfold_executable=str(FAKE_COLABFOLD),
    )

    summary = run_prediction(request)

    prediction = PredictionManifest.model_validate_json(
        summary.manifest_path.read_text(encoding="utf-8")
    )
    linked_nr = NrRunManifest.model_validate_json(
        nr_manifest_path.read_text(encoding="utf-8")
    )
    assert prediction.source_type == "nr_candidate"
    assert prediction.source_nr_run_id == nr_run_id
    assert prediction.source_accession == "WP_TEST_NR.1"
    assert prediction.prediction_id in linked_nr.child_prediction_ids
