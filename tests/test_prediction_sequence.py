from __future__ import annotations

from pathlib import Path

import pytest

from bur01_structmine.colabfold_runtime import (
    ColabFoldRuntime,
    build_msa_arguments,
    build_prediction_arguments,
    preset_parameters,
    to_wsl_path,
)
from bur01_structmine.config import project_root
from bur01_structmine.prediction_sequence import (
    SequenceValidationError,
    map_8hfw_catalytic_sites,
    parse_single_fasta,
    read_single_fasta,
)


def test_single_fasta_validation_collects_input_errors() -> None:
    with pytest.raises(SequenceValidationError) as exc_info:
        parse_single_fasta(">bad\nXXXXXXBZ\n")
    message = str(exc_info.value)
    assert "20-1000" in message
    assert "unsupported residue" in message
    assert "5%" in message


def test_multiple_fasta_records_are_rejected() -> None:
    with pytest.raises(SequenceValidationError, match="exactly one FASTA"):
        parse_single_fasta(">one\nACDEFGHIKLMNPQRSTVWY\n>two\nACDEFGHIKLMNPQRSTVWY\n")


def test_8hfw_catalytic_sites_map_to_seed_sequence() -> None:
    record = read_single_fasta(project_root() / "data" / "query" / "Bur01m_seed.fasta")
    mapping = map_8hfw_catalytic_sites(record.sequence)
    assert mapping["complete"] is True
    assert set(mapping["valid_positions"]) == {"S111", "D241", "H274"}
    assert all(site["status"] == "mapped" for site in mapping["sites"].values())


def test_manual_catalytic_mismatch_is_explicit() -> None:
    record = read_single_fasta(project_root() / "data" / "query" / "Bur01m_seed.fasta")
    mapping = map_8hfw_catalytic_sites(record.sequence, override={"S111": 1})
    assert mapping["complete"] is False
    assert mapping["sites"]["S111"]["source"] == "manual"
    assert mapping["sites"]["S111"]["status"] == "identity_mismatch"


def test_wsl_path_and_preset_command_construction() -> None:
    converted = to_wsl_path(Path("E:/Structure mining/query.fasta"))
    assert converted == "/mnt/e/Structure mining/query.fasta"
    standard = preset_parameters("standard")
    assert standard.num_models == 5
    assert standard.num_recycle == 3
    assert standard.num_seeds == 1
    assert standard.num_relax == 1
    msa = build_msa_arguments("input.fasta", "out", use_templates=True)
    assert "--msa-only" in msa
    assert "--templates" in msa
    command = build_prediction_arguments("input.fasta", "out", standard, use_templates=False)
    assert command[command.index("--num-models") + 1] == "5"
    assert command[command.index("--num-relax") + 1] == "1"
    assert "--disable-unified-memory" in command
    assert "--use-gpu-relax" not in command


def test_parameter_checksums_reuse_verified_provenance(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        ColabFoldRuntime,
        "parameter_provenance",
        lambda self: {"parameter_files": {"params_model_1_ptm.npz": "abc123"}},
    )
    assert ColabFoldRuntime().parameter_checksums() == {
        "params_model_1_ptm.npz": "abc123"
    }
