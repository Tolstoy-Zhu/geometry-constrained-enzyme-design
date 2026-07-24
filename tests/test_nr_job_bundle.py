from __future__ import annotations

import json
import subprocess
import sys
import zipfile
from pathlib import Path

import pandas as pd
import pytest

from bur01_structmine.config import project_root
from bur01_structmine.mmseqs_runner import MMSEQS_FIELDS
from bur01_structmine.models import NrExpansionRequest, RunManifest
from bur01_structmine.nr_job_bundle import (
    _safe_extract,
    _validate_database_manifest,
    _validated_run_id,
    export_nr_job_bundle,
    import_nr_job_bundle,
)
from bur01_structmine.prediction_sequence import read_single_fasta


def _structure_run(tmp_path: Path) -> tuple[Path, str]:
    structure_root = tmp_path / "runs"
    run_id = "source-structure-run"
    run_dir = structure_root / run_id
    run_dir.mkdir(parents=True)
    query = run_dir / "query.fasta"
    query.write_text(
        (project_root() / "data" / "query" / "Bur01m_seed.fasta").read_text(
            encoding="utf-8"
        ),
        encoding="utf-8",
    )
    manifest = RunManifest(
        run_id=run_id,
        status="completed",
        config_path=str(project_root() / "config.yaml"),
        query_pdb_id="8HFW",
        query_chain="A",
        artifacts={"query_fasta": str(query.resolve())},
    )
    (run_dir / "run_manifest.json").write_text(
        manifest.model_dump_json(indent=2),
        encoding="utf-8",
    )
    return structure_root, run_id


def test_exported_full_nr_bundle_is_credential_free_and_fixed(
    tmp_path: Path,
) -> None:
    structure_root, source_run = _structure_run(tmp_path)
    summary = export_nr_job_bundle(
        NrExpansionRequest(
            source_run_id=source_run,
            config_path=str(project_root() / "configs" / "nr_full_server.yaml"),
            nr_run_id="nr-full-test",
            output_root=str(tmp_path / "nr_runs"),
            structure_run_root=str(structure_root),
            database_scope="full_nr",
            execution_mode="job_bundle",
        )
    )

    archive = summary.artifacts["job_bundle_zip"]
    with zipfile.ZipFile(archive) as bundle:
        names = set(bundle.namelist())
        run_script = bundle.read("run_nr_pipeline.sh").decode("utf-8")
        environment = bundle.read("job.env").decode("utf-8")
        prepare_script = bundle.read("prepare_nr_database.sh").decode("utf-8")
    assert {
        "seeds.fasta",
        "bur01_query.fasta",
        "job_metadata.json",
        "resource_precheck.sh",
        "prepare_nr_database.sh",
        "run_nr_pipeline.sh",
        "submit.slurm",
        "submit.pbs",
        "BUNDLE_SHA256SUMS",
    } <= names
    assert "--num-iterations 2" in run_script
    assert "--max-seqs 2500" in run_script
    assert "--search-type 1" in run_script
    assert "--min-seq-id 0.5" in run_script
    assert 'DIRECT_MAX_SEQS="$REP_COUNT"' in run_script
    assert '--max-seqs "$DIRECT_MAX_SEQS"' in run_script
    assert 'MMSEQS" databases NR' in prepare_script
    assert "password=" not in environment.lower()
    assert "token=" not in environment.lower()
    assert summary.warnings
    assert summary.artifacts["run_checksums"].is_file()


def _write_server_result(
    destination: Path,
    nr_run_id: str,
    parent_run_id: str,
    fake_mmseqs_executable: Path,
) -> Path:
    destination.mkdir(parents=True)
    query = read_single_fasta(project_root() / "data" / "query" / "Bur01m_seed.fasta")
    representative_id = "nr_000000001"
    values = {
        "query": "Bur01_8HFW_A",
        "target": "WP_FULL_NR.1",
        "pident": 100.0,
        "fident": 1.0,
        "nident": len(query.sequence),
        "alnlen": len(query.sequence),
        "qlen": len(query.sequence),
        "tlen": len(query.sequence),
        "qstart": 1,
        "qend": len(query.sequence),
        "tstart": 1,
        "tend": len(query.sequence),
        "qcov": 1.0,
        "tcov": 1.0,
        "bits": 900.0,
        "evalue": 1e-100,
        "qaln": query.sequence,
        "taln": query.sequence,
        "tseq": query.sequence,
        "theader": "WP_FULL_NR.1 Bur01-like hydrolase [Full NR species]",
    }
    raw = destination / "nr_hits.raw.tsv"
    raw.write_text(
        "\t".join(str(values[field]) for field in MMSEQS_FIELDS) + "\n",
        encoding="utf-8",
    )
    pd.DataFrame([{"sequence_id": representative_id, **values}]).to_csv(
        destination / "nr_hits.deduplicated.tsv",
        sep="\t",
        index=False,
    )
    (destination / "nr_s50_rep_seq.fasta").write_text(
        f">{representative_id}\n{query.sequence}\n",
        encoding="ascii",
    )
    (destination / "nr_s50_cluster.tsv").write_text(
        f"{representative_id}\t{representative_id}\n",
        encoding="ascii",
    )
    direct_values = {
        **values,
        "target": representative_id,
        "theader": representative_id,
    }
    (destination / "bur01_vs_representatives.tsv").write_text(
        "\t".join(str(direct_values[field]) for field in MMSEQS_FIELDS) + "\n",
        encoding="utf-8",
    )
    (destination / "job_metadata.json").write_text(
        json.dumps(
            {
                "project_version": "0.4.0",
                "nr_run_id": nr_run_id,
                "parent_structure_run_id": parent_run_id,
                "database_scope": "full_nr",
                "paper_preset": "hong2025",
                "parameters": {
                    "iterations": 2,
                    "max_seqs": 2500,
                    "search_type": 1,
                    "cluster_min_sequence_identity": 0.5,
                },
            }
        ),
        encoding="utf-8",
    )
    (destination / "database_manifest.json").write_text(
        json.dumps(
            {
                "database_scope": "full_nr",
                "database_name": "NCBI NR",
                "database_snapshot_date": "2026-07-24",
                "database_prefix_name": "nr",
                "mmseqs_version": "18.8cc5c",
                "file_count": 12,
                "total_bytes": 123456789,
                "dbtype_sha256": "0" * 64,
                "source_command": "mmseqs databases NR",
            }
        ),
        encoding="utf-8",
    )
    script = project_root() / "scripts" / "nr_server_postprocess.py"
    subprocess.run(
        [
            sys.executable,
            str(script),
            "finalize",
            str(destination),
            nr_run_id,
            "NCBI NR full database",
            "2026-07-24",
            "--mmseqs",
            str(fake_mmseqs_executable),
        ],
        check=True,
    )
    archive = destination.parent / f"{nr_run_id}_server_result.zip"
    subprocess.run(
        [
            sys.executable,
            str(script),
            "archive",
            str(destination),
            str(archive),
        ],
        check=True,
    )
    return archive


def test_imported_full_nr_result_is_verified_and_sequence_only(
    tmp_path: Path,
    fake_mmseqs_executable: Path,
) -> None:
    structure_root, source_run = _structure_run(tmp_path)
    output_root = tmp_path / "nr_runs"
    nr_run_id = "nr-full-import-test"
    export_nr_job_bundle(
        NrExpansionRequest(
            source_run_id=source_run,
            config_path=str(project_root() / "configs" / "nr_full_server.yaml"),
            nr_run_id=nr_run_id,
            output_root=str(output_root),
            structure_run_root=str(structure_root),
            database_scope="full_nr",
            execution_mode="job_bundle",
        )
    )
    archive = _write_server_result(
        tmp_path / "server_output",
        nr_run_id,
        source_run,
        fake_mmseqs_executable,
    )

    summary = import_nr_job_bundle(
        archive,
        output_root=output_root,
        structure_run_root=structure_root,
    )

    assert summary.status == "completed"
    assert summary.candidate_count == 1
    candidates = pd.read_csv(summary.artifacts["candidates_tsv"], sep="\t")
    assert candidates.iloc[0]["wase_tier"] == "Not evaluated - structure required"
    assert candidates.iloc[0]["wase_scope"] == "not_evaluated_sequence_only"
    manifest = json.loads(summary.manifest_path.read_text(encoding="utf-8"))
    assert manifest["database_scope"] == "full_nr"
    assert manifest["database_snapshot"]["validated"] is True
    assert manifest["database_snapshot"]["snapshot_date"] == "2026-07-24"
    assert Path(manifest["artifacts"]["run_checksums"]).is_file()


def test_result_archive_rejects_path_traversal(tmp_path: Path) -> None:
    archive = tmp_path / "unsafe.zip"
    with zipfile.ZipFile(archive, "w") as output:
        output.writestr("../escape.txt", "unsafe")

    with pytest.raises(ValueError, match="Unsafe path"):
        _safe_extract(archive, tmp_path / "extract")


def test_result_archive_rejects_windows_path_traversal(tmp_path: Path) -> None:
    archive = tmp_path / "unsafe-windows.zip"
    with zipfile.ZipFile(archive, "w") as output:
        output.writestr(r"..\escape.txt", "unsafe")

    with pytest.raises(ValueError, match="Unsafe path"):
        _safe_extract(archive, tmp_path / "extract")


def test_import_identifiers_and_database_manifest_are_strict(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Invalid nr_run_id"):
        _validated_run_id("../outside")

    database = tmp_path / "database_manifest.json"
    database.write_text(
        json.dumps(
            {
                "database_scope": "full_nr",
                "database_snapshot_date": "2026-07-24",
                "source_command": "mmseqs databases NR",
                "file_count": 1,
                "total_bytes": 100,
                "dbtype_sha256": "0" * 64,
                "mmseqs_version": "different",
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="different MMseqs2 versions"):
        _validate_database_manifest(
            database,
            {
                "database_snapshot_date": "2026-07-24",
                "mmseqs_version": "18.8cc5c",
            },
        )
