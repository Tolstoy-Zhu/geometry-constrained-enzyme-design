"""Transparent candidate triage and WASE-screen v0.1 scoring."""
from __future__ import annotations

import math
import re
from typing import Any

import pandas as pd

WASE_VERSION = "0.1.0"
FORMAL_WASE_MISSING = [
    "tunnel_specialization",
    "alcohol_retention",
    "water_exclusion",
    "product_escape",
    "trap_channel_restriction",
    "substrate_compatibility",
    "fermentation_source_prior",
    "adaptive_site_enrichment",
    "ancestral_gain",
]


def _present(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, float) and math.isnan(value):
        return False
    return str(value).strip().lower() not in {"", "nan", "none", "null"}


def _float(value: Any) -> float | None:
    if not _present(value):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _truthy(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "reviewed"}


def _clamp(value: float, low: float = 0.0, high: float = 1.0) -> float:
    return max(low, min(high, value))


def has_serine_hydrolase_motif(seq: str) -> bool:
    """Return whether a sequence contains a G-X-S-X-G-like motif."""
    if not _present(seq):
        return False
    clean = re.sub(r"[^A-Za-z]", "", str(seq)).upper()
    return re.search(r"G[A-Z]S[A-Z]G", clean) is not None


def _structure_fit(row: pd.Series) -> tuple[float, float, list[str]]:
    score = 0.0
    coverage = 0.0
    missing: list[str] = []
    qtm = _float(row.get("qtmscore"))
    ttm = _float(row.get("ttmscore"))
    if qtm is not None and ttm is not None:
        score += 25.0 * _clamp(min(qtm, ttm))
        coverage += 25.0
    else:
        missing.append("bidirectional_tm_score")

    qcov = _float(row.get("qcov"))
    tcov = _float(row.get("tcov"))
    if qcov is not None and tcov is not None:
        score += 10.0 * math.sqrt(_clamp(qcov) * _clamp(tcov))
        coverage += 10.0
    else:
        missing.append("alignment_coverage")

    rmsd = _float(row.get("rmsd"))
    if rmsd is not None:
        score += 10.0 * _clamp(1.0 - rmsd / 5.0)
        coverage += 10.0
    else:
        missing.append("alignment_rmsd")
    return score, coverage, missing


def _catalytic_core(row: pd.Series) -> tuple[float, float, list[str]]:
    score = 0.0
    coverage = 0.0
    missing: list[str] = []
    sequence = row.get("final_protein_sequence", "")
    if _present(sequence):
        score += 5.0 if has_serine_hydrolase_motif(str(sequence)) else 0.0
        coverage += 5.0
    else:
        missing.append("protein_sequence_for_motif")

    mapped_count = _float(row.get("mapped_catalytic_identity_count"))
    if mapped_count is not None:
        score += 15.0 * _clamp(mapped_count / 3.0)
        coverage += 15.0
    else:
        missing.append("mapped_catalytic_triad")

    geometry_contacts = _float(row.get("catalytic_geometry_contacts"))
    if geometry_contacts is not None:
        score += 10.0 * _clamp(geometry_contacts / 2.0)
        coverage += 10.0
    else:
        missing.append("catalytic_geometry")
    return score, coverage, missing


def _annotation_confidence(row: pd.Series) -> tuple[float, float, list[str]]:
    score = 0.0
    coverage = 0.0
    missing: list[str] = []
    names = " ".join(
        str(row.get(column, ""))
        for column in (
            "rcsb_protein_name",
            "uniprot_protein_name",
            "ncbi_protein_name",
            "alphafold_protein_description",
        )
        if _present(row.get(column, ""))
    ).lower()
    if names:
        keywords = ("esterase", "hydrolase", "lipase", "depolymerase", "arylesterase")
        score += 5.0 if any(keyword in names for keyword in keywords) else 0.0
        coverage += 5.0
    else:
        missing.append("functional_annotation")

    experimental = _present(row.get("pdb_id")) or str(row.get("annotation_source", "")).lower() == "pdb"
    reviewed = _truthy(row.get("uniprot_is_reviewed"))
    if experimental or _present(row.get("uniprot_is_reviewed")):
        score += 5.0 if experimental or reviewed else 0.0
        coverage += 5.0
    else:
        missing.append("review_or_experimental_status")

    plddt = _float(row.get("alphafold_global_plddt"))
    if experimental:
        score += 5.0
        coverage += 5.0
    elif plddt is not None:
        score += 5.0 * _clamp((plddt - 50.0) / 40.0)
        coverage += 5.0
    else:
        missing.append("structure_confidence")
    return score, coverage, missing


def _experimental_readiness(row: pd.Series) -> tuple[float, float, list[str]]:
    score = 0.0
    coverage = 0.0
    missing: list[str] = []
    sequence = str(row.get("final_protein_sequence", ""))
    if _present(sequence):
        score += 4.0 if "X" not in sequence.upper() else 0.0
        coverage += 4.0
    else:
        missing.append("complete_protein_sequence")

    length = _float(row.get("final_protein_length"))
    if length is not None:
        score += 2.0 if 180 <= length <= 500 else 0.0
        coverage += 2.0
    else:
        missing.append("protein_length")

    tcov = _float(row.get("tcov"))
    if tcov is not None:
        score += 2.0 if tcov >= 0.80 else 0.0
        coverage += 2.0
    else:
        missing.append("structure_sequence_coverage")

    accession_fields = (
        "target_uniprot_accession",
        "uniprot_accessions",
        "ncbi_protein_accession",
        "pdb_id",
        "target",
    )
    accession_available = any(_present(row.get(column)) for column in accession_fields)
    score += 2.0 if accession_available else 0.0
    coverage += 2.0
    if not accession_available:
        missing.append("stable_accession")
    return score, coverage, missing


def score_wase_screen(row: pd.Series) -> dict[str, Any]:
    """Score currently observable screening evidence without claiming activity prediction."""
    structure, structure_coverage, structure_missing = _structure_fit(row)
    catalytic, catalytic_coverage, catalytic_missing = _catalytic_core(row)
    annotation, annotation_coverage, annotation_missing = _annotation_confidence(row)
    readiness, readiness_coverage, readiness_missing = _experimental_readiness(row)
    score = structure + catalytic + annotation + readiness
    evidence_coverage = structure_coverage + catalytic_coverage + annotation_coverage + readiness_coverage
    uncertainty = 100.0 - evidence_coverage

    qtm = _float(row.get("qtmscore")) or 0.0
    ttm = _float(row.get("ttmscore")) or 0.0
    if evidence_coverage < 60.0:
        tier = "Incomplete evidence"
    elif (
        score >= 70.0
        and evidence_coverage >= 75.0
        and min(qtm, ttm) >= 0.50
        and catalytic >= 20.0
    ):
        tier = "High priority"
    elif score >= 50.0:
        tier = "Medium priority"
    else:
        tier = "Low priority"

    missing = structure_missing + catalytic_missing + annotation_missing + readiness_missing
    missing.extend(FORMAL_WASE_MISSING)
    return {
        "structure_fit_score": round(structure, 3),
        "catalytic_core_score": round(catalytic, 3),
        "annotation_confidence_score": round(annotation, 3),
        "experimental_readiness_score": round(readiness, 3),
        "wase_screen_score": round(score, 3),
        "wase_evidence_coverage": round(evidence_coverage, 3),
        "wase_screen_uncertainty": round(uncertainty, 3),
        "wase_tier": tier,
        "wase_missing_evidence": ";".join(dict.fromkeys(missing)),
        "wase_version": WASE_VERSION,
        "wase_scope": "screening_only_not_formal_wase",
    }


def assign_recommendation(row: pd.Series) -> str:
    """Assign an experiment-facing recommendation after WASE-screen scoring."""
    seq = str(row.get("final_protein_sequence", ""))
    if not _present(seq):
        return "Hold: no protein sequence retrieved"
    if "X" in seq.upper():
        return "Hold: sequence contains X"
    tier = str(row.get("wase_tier", ""))
    if tier == "High priority":
        return "Priority: structurally and catalytically supported"
    if tier == "Medium priority":
        return "Explore: useful screening evidence"
    if tier == "Incomplete evidence":
        return "Review: collect missing structural or annotation evidence"
    return "Low priority: weak current screening evidence"


def add_candidate_scores(df: pd.DataFrame) -> pd.DataFrame:
    """Add motif, WASE-screen, evidence-gap, and recommendation columns."""
    out = df.copy()
    if "final_protein_sequence" not in out.columns:
        out["final_protein_sequence"] = pd.Series(dtype=str)
    out["has_GXSXG_like_motif"] = out["final_protein_sequence"].apply(has_serine_hydrolase_motif)
    if out.empty:
        for column in score_wase_screen(pd.Series(dtype=object)):
            out[column] = pd.Series(dtype=object)
        out["synthesis_recommendation"] = pd.Series(dtype=str)
        return out

    score_rows = out.apply(score_wase_screen, axis=1, result_type="expand")
    for column in score_rows.columns:
        out[column] = score_rows[column]
    out["synthesis_recommendation"] = out.apply(assign_recommendation, axis=1)
    return out
