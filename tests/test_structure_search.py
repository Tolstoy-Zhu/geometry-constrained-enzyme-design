from __future__ import annotations

import math
from pathlib import Path

import pytest

from bur01_structmine.config import project_root
from bur01_structmine.structure_search import compare_ca_traces, load_ca_trace, load_ca_traces
from bur01_structmine.visualization import align_structure_overlay


QUERY = project_root() / "data" / "query" / "Bur01m.pdb"


def test_self_alignment_is_identity_and_preserves_chain() -> None:
    trace = load_ca_trace(QUERY, chain_id="A")
    row = compare_ca_traces(trace, trace)
    assert row["query_chain"] == "A"
    assert row["target_chain"] == "A"
    assert row["qtmscore"] > 0.999
    assert row["ttmscore"] > 0.999
    assert row["rmsd"] < 0.001
    assert row["mapped_catalytic_residues"] == "S111;D241;H274"
    assert row["mapped_catalytic_identity_count"] == 3
    assert row["search_engine"] == "tmtools_tmalign"
    assert math.isnan(row["evalue"])
    assert math.isnan(row["bits"])
    assert math.isnan(row["prob"])


def test_multichain_structure_requires_explicit_selection() -> None:
    multichain = project_root() / "data" / "benchmark" / "structures" / "1XKL.pdb"
    if not multichain.exists():
        pytest.skip("benchmark structures have not been prepared")
    traces = load_ca_traces(multichain)
    assert {trace.chain_id for trace in traces} == {"A", "B", "C", "D"}
    with pytest.raises(ValueError, match="Multiple protein chains"):
        load_ca_trace(multichain)
    selected = load_ca_trace(multichain, chain_id="C")
    assert selected.chain_id == "C"
    assert len(selected.residues) == 248


def test_unknown_query_chain_is_rejected() -> None:
    with pytest.raises(ValueError, match="was not found"):
        load_ca_trace(QUERY, chain_id="Z")


def test_structure_overlay_serializes_transformed_models() -> None:
    overlay = align_structure_overlay(QUERY, QUERY, query_chain="A", target_chain="A")
    assert overlay.qtmscore > 0.999
    assert overlay.rmsd < 0.001
    assert "ATOM" in overlay.query_pdb
    assert "ATOM" in overlay.target_pdb
    assert overlay.mapped_target_residue_numbers == (111, 241, 274)
