from __future__ import annotations

import copy
import json

import pandas as pd
import pytest

from bur01_structmine.config import project_root
from bur01_structmine.rescoring import rescore_run
from bur01_structmine.scoring import read_score_display, score_wase_screen
from bur01_structmine.structure_search import _map_catalytic_residues, load_ca_trace, DEFAULT_CATALYTIC_QUERY_RESIDUES
from bur01_structmine.results_parser import explain_filter_hits, filter_hits


def test_unrelated_motif_no_catalytic_credit():
    row = pd.Series({"final_protein_sequence": "GASAGAAAAAAAAA", "mapped_catalytic_ser_motif": "AASAA"})
    assert score_wase_screen(row)["catalytic_core_score"] == 0
    row["mapped_catalytic_ser_motif"] = "GASAG"
    assert score_wase_screen(row)["catalytic_core_score"] == 5


def test_unmapped_motif_and_legacy_contacts_are_missing_not_negative():
    result = score_wase_screen(pd.Series({"final_protein_sequence": "GASAG", "catalytic_geometry_contacts": 0}))
    assert "mapped_ser_motif" in result["wase_missing_evidence"]
    assert "catalytic_ser_his_distance" in result["wase_missing_evidence"]
    assert result["wase_prediction_error"] is None


def test_missing_sidechain_and_out_of_range_are_distinct():
    trace = load_ca_trace(project_root() / "data/query/Bur01m.pdb", chain_id="A")
    mapping = dict(enumerate(range(len(trace.residues))))
    target = copy.deepcopy(trace)
    serine = next(res for res in target.residues if res.residue_number == 111)
    del serine.atoms["OG"]
    missing = _map_catalytic_residues(trace, target, mapping, DEFAULT_CATALYTIC_QUERY_RESIDUES)
    assert missing["catalytic_geometry_status"] == "not_assessable"
    assert missing["catalytic_geometry_contacts"] is None
    assert missing["catalytic_ser_his_distance"] is None
    assert "Ser:missing_sidechain_atoms" in missing["catalytic_geometry_missing"]
    target = copy.deepcopy(trace)
    next(res for res in target.residues if res.residue_number == 111).atoms["OG"] += 100
    bad = _map_catalytic_residues(trace, target, mapping, DEFAULT_CATALYTIC_QUERY_RESIDUES)
    assert bad["catalytic_geometry_status"] == "out_of_range"
    assert bad["catalytic_ser_his_distance"] > 5


def test_history_compatibility_does_not_mutate_row():
    original = {"wase_version": "0.1.0", "wase_screen_uncertainty": 35}
    assert read_score_display(original)["information_missingness"] == 35
    assert original == {"wase_version": "0.1.0", "wase_screen_uncertainty": 35}


def test_filter_explanations_match_existing_filter():
    raw = pd.DataFrame([{"target": "ok", "alntmscore": .8, "rmsd": 1, "qcov": .8, "tcov": .9, "alnlen": 250, "fident": .5},
                        {"target": "negative", "alntmscore": .2, "rmsd": 10, "qcov": .4, "tcov": .2, "alnlen": 50, "fident": .1}])
    explanation = explain_filter_hits(raw)
    assert explanation.loc[explanation.structural_filter_pass, "target"].tolist() == filter_hits(raw).target.tolist()
    assert "alntmscore<0.45" in explanation.iloc[1].filter_exclusion_reasons


def test_rescore_creates_child_with_hash_and_deltas(tmp_path):
    source = tmp_path / "old"
    source.mkdir()
    query = project_root() / "data/query/Bur01m.pdb"
    table_path = source / "annotated.tsv"
    pd.DataFrame([{"target": "8HFW_A", "wase_version": "0.1.0", "wase_screen_score": 50,
                   "qtmscore": 1, "ttmscore": 1, "qcov": 1, "tcov": 1, "rmsd": 0}]).to_csv(table_path, sep="\t", index=False)
    manifest_path = source / "run_manifest.json"
    manifest_path.write_text(json.dumps({"run_id": "old", "query_path": str(query), "query_pdb_id": "8HFW", "query_chain": "A",
                                         "artifacts": {"annotated_tsv": str(table_path)}}))
    before = {p.name: p.read_bytes() for p in source.iterdir()}
    destination = rescore_run(source, tmp_path / "rescores")
    assert before == {p.name: p.read_bytes() for p in source.iterdir()}
    row = pd.read_csv(destination / "rescored.tsv", sep="\t").iloc[0]
    assert row.source_run_id == "old"
    assert row.wase_version == "0.2.0"
    assert row.wase_score_delta == pytest.approx(row.wase_screen_score - 50)
    assert pd.isna(row.wase_screen_uncertainty)
