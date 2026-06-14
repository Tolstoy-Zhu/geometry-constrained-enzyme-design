"""Candidate scoring utilities."""
from __future__ import annotations

import re

import pandas as pd


def has_serine_hydrolase_motif(seq: str) -> bool:
    """Check a simple serine hydrolase motif pattern.

    Many alpha/beta hydrolases contain a nucleophilic Ser in a small-residue
    motif such as GXSXG. This is only a rough filter, not proof of activity.
    """
    if not seq or str(seq).lower() == "nan":
        return False
    return re.search(r"G[A-Z]S[A-Z]G", seq) is not None


def assign_recommendation(row: pd.Series) -> str:
    """Assign a simple synthesis recommendation."""
    seq = str(row.get("final_protein_sequence", ""))
    length = row.get("final_protein_length", None)
    has_x = row.get("has_unknown_X", False)
    name = " ".join([
        str(row.get("rcsb_protein_name", "")),
        str(row.get("uniprot_protein_name", "")),
    ]).lower()

    if not seq or seq.lower() == "nan":
        return "Hold: no protein sequence retrieved"
    if has_x is True or "X" in seq:
        return "Hold: sequence contains X"
    try:
        length_int = int(length)
    except Exception:
        length_int = len(seq)
    if length_int < 180 or length_int > 500:
        return "Check: unusual protein length"

    function_hit = any(k in name for k in ["esterase", "hydrolase", "lipase", "depolymerase", "arylesterase"])
    motif_hit = has_serine_hydrolase_motif(seq)

    if function_hit and motif_hit:
        return "Priority: function annotation and motif support"
    if function_hit:
        return "Priority: function annotation support"
    if motif_hit:
        return "Explore: motif support only"
    return "Explore: structure similarity only"


def add_candidate_scores(df: pd.DataFrame) -> pd.DataFrame:
    """Add motif and recommendation columns."""
    out = df.copy()
    out["has_GXSXG_like_motif"] = out.get("final_protein_sequence", "").apply(has_serine_hydrolase_motif)
    out["synthesis_recommendation"] = out.apply(assign_recommendation, axis=1)
    return out
