from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
import yaml

from bur01_structmine.config import project_root
from bur01_structmine.models import (
    NrExpansionRequest,
    NrNoveltyRequest,
    PredictionManifest,
    RunManifest,
)
from bur01_structmine.nr_pipeline import (
    build_nr_demo_snapshot,
    deduplicate_hits,
    load_nr_config,
    map_catalytic_alignment,
    run_nr_expansion,
    run_nr_novelty,
)
from bur01_structmine.ncbi_blast import NcbiBlastHit
from bur01_structmine.prediction_sequence import read_single_fasta


def _write_structure_run(tmp_path: Path) -> tuple[Path, str]:
    run_root = tmp_path / "runs"
    run_id = "structure-seed-run"
    run_dir = run_root / run_id
    run_dir.mkdir(parents=True)
    query_fasta = run_dir / "query.fasta"
    query_fasta.write_text(
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
        artifacts={"query_fasta": str(query_fasta.resolve())},
    )
    (run_dir / "run_manifest.json").write_text(
        manifest.model_dump_json(indent=2),
        encoding="utf-8",
    )
    return run_root, run_id


def _write_nr_database(tmp_path: Path) -> Path:
    query = read_single_fasta(project_root() / "data" / "query" / "Bur01m_seed.fasta")
    negative = "A" * len(query.sequence)
    path = tmp_path / "demo_nr.fasta"
    path.write_text(
        (
            f">WP_POSITIVE.1 Bur01-like esterase [Positive species]\n{query.sequence}\n"
            f">P00520 distant negative control [Negative species]\n{negative}\n"
        ),
        encoding="ascii",
    )
    return path


def _write_config(
    tmp_path: Path,
    executable: Path,
    database_path: Path,
) -> Path:
    config = {
        "database": {
            "scope": "nr_remote_snapshot",
            "name": "fixed fake NR snapshot",
            "path": str(database_path),
        },
        "runtime": {
            "backend": "native",
            "executable": str(executable),
        },
        "search": {
            "paper_preset": "hong2025",
            "iterations": 2,
            "max_seqs": 2500,
            "search_type": 1,
            "evalue": 1e-5,
            "threads": 2,
            "cluster_min_sequence_identity": 0.5,
        },
        "filter": {
            "max_evalue": 1e-5,
            "min_qcov": 0.60,
            "min_tcov": 0.40,
            "min_length": 180,
            "max_length": 500,
            "max_x_fraction": 0.05,
            "required_catalytic_identities": 3,
        },
    }
    path = tmp_path / "nr_test.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return path


def test_nr_pipeline_fake_executable_end_to_end(
    tmp_path: Path,
    fake_mmseqs_executable: Path,
) -> None:
    structure_root, source_run_id = _write_structure_run(tmp_path)
    database = _write_nr_database(tmp_path)
    config = _write_config(tmp_path, fake_mmseqs_executable, database)
    output_root = tmp_path / "nr_runs"

    summary = run_nr_expansion(
        NrExpansionRequest(
            source_run_id=source_run_id,
            config_path=str(config),
            output_root=str(output_root),
            structure_run_root=str(structure_root),
            database_path=str(database),
            database_scope="nr_remote_snapshot",
        )
    )

    assert summary.status == "completed"
    assert summary.candidate_count == 2
    candidates = pd.read_csv(summary.artifacts["candidates_tsv"], sep="\t")
    positive = candidates[candidates["target_accession"] == "WP_POSITIVE.1"].iloc[0]
    negative = candidates[candidates["target_accession"] == "P00520"].iloc[0]
    assert bool(positive["af2_ready"])
    assert positive["mapped_catalytic_identity_count"] == 3
    assert not bool(negative["af2_ready"])
    assert "catalytic_triad" in negative["sequence_screen_status"]
    assert candidates["wase_screen_score"].isna().all()
    assert set(candidates["wase_scope"]) == {"not_evaluated_sequence_only"}

    parent = RunManifest.model_validate_json(
        (structure_root / source_run_id / "run_manifest.json").read_text(
            encoding="utf-8"
        )
    )
    assert summary.nr_run_id in parent.child_nr_run_ids
    manifest = json.loads(summary.manifest_path.read_text(encoding="utf-8"))
    assert manifest["database_scope"] == "nr_remote_snapshot"
    assert manifest["resolved_parameters"]["iterations"] == 2
    assert manifest["software_versions"]["mmseqs"] == "18.8cc5c"
    checksum_path = Path(manifest["artifacts"]["run_checksums"])
    assert checksum_path.is_file()
    for line in checksum_path.read_text(encoding="ascii").splitlines():
        digest, relative = line.split(maxsplit=1)
        assert len(digest) == 64
        assert (summary.run_dir / relative.strip()).is_file()


def test_deduplicate_prefers_best_evalue_then_bits() -> None:
    raw = pd.DataFrame(
        [
            {"target": "ref|WP_1.1|", "evalue": 1e-5, "bits": 100, "qcov": 0.8},
            {"target": "WP_1.1", "evalue": 1e-20, "bits": 80, "qcov": 0.7},
            {"target": "WP_2.1", "evalue": 1e-10, "bits": 70, "qcov": 0.7},
        ]
    )

    result = deduplicate_hits(raw)

    assert len(result) == 2
    assert result[result["target_accession"] == "WP_1.1"].iloc[0]["evalue"] == 1e-20


def test_deduplicate_merges_identical_sequences_across_accessions() -> None:
    raw = pd.DataFrame(
        [
            {"target": "WP_1.1", "tseq": "MSTAG", "evalue": 1e-20, "bits": 80},
            {"target": "WP_ALIAS.1", "tseq": "MSTAG", "evalue": 1e-10, "bits": 70},
        ]
    )

    result = deduplicate_hits(raw)

    assert len(result) == 1
    assert result.iloc[0]["target_accession"] == "WP_1.1"
    assert len(result.iloc[0]["sequence_sha256"]) == 64


def test_catalytic_mapping_reports_identity_mismatch() -> None:
    query = read_single_fasta(project_root() / "data" / "query" / "Bur01m_seed.fasta")
    baseline = map_catalytic_alignment(
        pd.Series(
            {
                "qstart": 1,
                "tstart": 1,
                "qaln": query.sequence,
                "taln": query.sequence,
            }
        )
    )
    serine_position = baseline["mapping"]["S111"]["target_position"]
    target = list(query.sequence)
    target[serine_position - 1] = "A"
    row = pd.Series(
        {
            "qstart": 1,
            "tstart": 1,
            "qaln": query.sequence,
            "taln": "".join(target),
        }
    )

    mapping = map_catalytic_alignment(row)

    assert mapping["identity_count"] == 2
    assert mapping["mapping"]["S111"]["status"] == "identity_mismatch"
    assert mapping["complete"] is False


def test_old_manifests_gain_compatible_nr_defaults() -> None:
    run = RunManifest.model_validate(
        {"run_id": "old", "config_path": "config.yaml"}
    )
    prediction = PredictionManifest.model_validate(
        {"prediction_id": "old-prediction"}
    )

    assert run.child_nr_run_ids == []
    assert prediction.source_nr_run_id == ""


def test_hong2025_config_rejects_parameter_drift(tmp_path: Path) -> None:
    config = yaml.safe_load(
        (project_root() / "configs" / "nr_local_demo.yaml").read_text(
            encoding="utf-8"
        )
    )
    config["search"]["max_seqs"] = 100
    path = tmp_path / "drifted.yaml"
    path.write_text(yaml.safe_dump(config), encoding="utf-8")

    with pytest.raises(ValueError, match="hong2025 preset parameters are fixed"):
        load_nr_config(path)


def test_nr_run_ids_cannot_escape_or_overwrite_output(
    tmp_path: Path,
    fake_mmseqs_executable: Path,
) -> None:
    with pytest.raises(ValueError, match="unsupported path"):
        NrExpansionRequest(source_run_id="../outside")

    structure_root, source_run_id = _write_structure_run(tmp_path)
    database = _write_nr_database(tmp_path)
    config = _write_config(tmp_path, fake_mmseqs_executable, database)
    output_root = tmp_path / "nr_runs"
    request = NrExpansionRequest(
        source_run_id=source_run_id,
        config_path=str(config),
        nr_run_id="fixed-run",
        output_root=str(output_root),
        structure_run_root=str(structure_root),
        database_path=str(database),
        database_scope="nr_remote_snapshot",
    )
    run_nr_expansion(request)

    with pytest.raises(FileExistsError):
        run_nr_expansion(request)


def test_nr_snapshot_preserves_sequence_fetch_errors(tmp_path: Path) -> None:
    class StubClient:
        def search(self, *args, **kwargs):
            query = read_single_fasta(
                project_root() / "data" / "query" / "Bur01m_seed.fasta"
            )

            def hit(accession: str) -> NcbiBlastHit:
                return NcbiBlastHit(
                    query=query.name,
                    target=accession,
                    target_accession=accession,
                    target_title=f"{accession} test",
                    pident=80.0,
                    fident=0.8,
                    nident=200,
                    alnlen=250,
                    qlen=len(query.sequence),
                    tlen=None,
                    qstart=1,
                    qend=250,
                    tstart=1,
                    tend=250,
                    qcov=250 / len(query.sequence),
                    tcov=None,
                    bits=300.0,
                    evalue=1e-50,
                    qaln=query.sequence[:250],
                    taln=query.sequence[:250],
                )

            return [hit("WP_GOOD.1"), hit("WP_MISSING.1")], {"database": "nr"}

        def fetch_target_sequences(self, hits, **kwargs):
            hits[0].target_sequence = "MSTAG" * 50
            hits[0].tlen = len(hits[0].target_sequence)
            hits[0].tcov = hits[0].alnlen / hits[0].tlen
            hits[1].annotation_error = "efetch did not return a sequence"
            return hits

    manifest_path = build_nr_demo_snapshot(
        project_root() / "data" / "query" / "Bur01m_seed.fasta",
        contact_email="researcher@example.org",
        confirm_ncbi=True,
        output_root=tmp_path / "snapshots",
        client=StubClient(),  # type: ignore[arg-type]
    )

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    hits = pd.read_csv(manifest["hits_path"], sep="\t")
    assert manifest["raw_hit_count"] == 2
    assert manifest["sequence_count"] == 1
    assert manifest["sequence_fetch_error_count"] == 1
    assert hits.loc[hits["target_accession"] == "WP_MISSING.1", "annotation_error"].iloc[0]


def test_nr_novelty_mmseqs_end_to_end_is_traceable(
    tmp_path: Path,
    fake_mmseqs_executable: Path,
) -> None:
    query_fasta = project_root() / "data" / "query" / "Bur01m_seed.fasta"
    database = _write_nr_database(tmp_path)
    manifest, manifest_path = run_nr_novelty(
        NrNoveltyRequest(
            fasta_path=str(query_fasta),
            novelty_id="novelty-test",
            output_root=str(tmp_path / "novelty"),
            backend="mmseqs",
            database_path=str(database),
            database_scope="custom_mmseqs",
            runtime_backend="native",
            mmseqs_executable=str(fake_mmseqs_executable),
        )
    )

    assert manifest.status == "completed"
    assert manifest.nearest_nr_accession == "WP_POSITIVE.1"
    assert manifest.nearest_nr_identity == pytest.approx(1.0)
    assert manifest.self_hit_warning
    assert Path(manifest.artifacts["run_checksums"]).is_file()
    persisted = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert persisted["database_scope"] == "custom_mmseqs"
    assert persisted["software_versions"]["mmseqs"] == "18.8cc5c"
