from __future__ import annotations

import pandas as pd

from bur01_structmine.benchmark import NEGATIVE_LABELS, POSITIVE_LABELS, validate_benchmark_tables


def _row(label: str, score: float, rmsd: float = 2.0) -> dict[str, object]:
    return {
        "target": f"target_{label}",
        "alntmscore": score,
        "qtmscore": score,
        "ttmscore": score,
        "rmsd": rmsd,
    }


def test_fixed_benchmark_acceptance() -> None:
    scores = {"8HFW": 1.0, "1XKL": 0.8, "6TJ2": 0.75, "Q6RYA0": 0.78, "Q9FVW3": 0.64, "P60174": 0.34, "P00520": 0.1}
    raw = pd.DataFrame([_row(label, score, 0.0 if label == "8HFW" else 2.0) for label, score in scores.items()])
    filtered = raw[raw["alntmscore"] >= 0.45].copy()
    report = validate_benchmark_tables(raw, filtered)
    assert report["passed"] is True
    assert not report["missing_controls"]


def test_benchmark_fails_when_negative_outranks_positive() -> None:
    scores = {label: 0.7 for label in POSITIVE_LABELS}
    scores.update({NEGATIVE_LABELS[0]: 0.9, NEGATIVE_LABELS[1]: 0.1})
    scores["8HFW"] = 1.0
    raw = pd.DataFrame([_row(label, score, 0.0 if label == "8HFW" else 2.0) for label, score in scores.items()])
    filtered = raw[raw["alntmscore"] >= 0.45].copy()
    report = validate_benchmark_tables(raw, filtered)
    assert report["passed"] is False
    ranking_check = next(check for check in report["checks"] if check["name"] == "positive_ranking")
    assert ranking_check["passed"] is False


def test_predicted_query_uses_similarity_not_identity_gate() -> None:
    scores = {"8HFW": 0.8, "1XKL": 0.74, "6TJ2": 0.72, "Q6RYA0": 0.70, "Q9FVW3": 0.62, "P60174": 0.31, "P00520": 0.12}
    raw = pd.DataFrame([_row(label, score, 2.1) for label, score in scores.items()])
    filtered = raw[raw["alntmscore"] >= 0.45].copy()
    report = validate_benchmark_tables(raw, filtered, self_tm_threshold=0.75, require_self_rmsd=False)
    assert report["passed"] is True
    rmsd_check = next(check for check in report["checks"] if check["name"] == "self_rmsd")
    assert rmsd_check["passed"] is True
