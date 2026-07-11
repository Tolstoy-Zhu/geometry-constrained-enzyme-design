from __future__ import annotations

import math

import pandas as pd

from bur01_structmine.scoring import FORMAL_WASE_MISSING, score_wase_screen


def test_wase_screen_uses_documented_formula_and_high_priority_gate() -> None:
    row = pd.Series(
        {
            "qtmscore": 0.8,
            "ttmscore": 0.6,
            "qcov": 0.9,
            "tcov": 0.4,
            "rmsd": 2.0,
            "final_protein_sequence": "MAGASAGAA",
            "final_protein_length": 250,
            "mapped_catalytic_identity_count": 3,
            "catalytic_geometry_contacts": 2,
            "uniprot_protein_name": "carboxylesterase",
            "uniprot_is_reviewed": True,
            "alphafold_global_plddt": 90,
            "target": "Q6RYA0",
        }
    )
    result = score_wase_screen(row)
    expected_structure = 25 * 0.6 + 10 * math.sqrt(0.9 * 0.4) + 10 * (1 - 2 / 5)
    assert result["structure_fit_score"] == round(expected_structure, 3)
    assert result["catalytic_core_score"] == 30.0
    assert result["annotation_confidence_score"] == 15.0
    assert result["experimental_readiness_score"] == 8.0
    assert result["wase_screen_score"] == round(expected_structure + 53, 3)
    assert result["wase_evidence_coverage"] == 100.0
    assert result["wase_tier"] == "High priority"
    assert result["wase_version"] == "0.1.0"
    for missing in FORMAL_WASE_MISSING:
        assert missing in result["wase_missing_evidence"]


def test_high_priority_requires_75_percent_evidence_coverage() -> None:
    row = pd.Series(
        {
            "qtmscore": 1.0,
            "ttmscore": 1.0,
            "qcov": 1.0,
            "tcov": 1.0,
            "rmsd": 0.0,
            "mapped_catalytic_identity_count": 3,
            "catalytic_geometry_contacts": 2,
            "target": "candidate",
        }
    )
    result = score_wase_screen(row)
    assert result["wase_screen_score"] >= 70
    assert result["wase_evidence_coverage"] < 75
    assert result["wase_tier"] != "High priority"
