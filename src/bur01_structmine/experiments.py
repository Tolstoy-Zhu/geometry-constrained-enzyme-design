"""Condition-aware experimental evaluation shared by CLI and web views."""
from __future__ import annotations

import hashlib
import io
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from pydantic import ValidationError

from .models import EvaluationSummary, ExperimentRecord

METRIC_UNITS = {
    "activity": ("U/mL", "U/mg", "U", "kat", "kat/kg", "umol/min", "mmol/L/h"),
    "relative_activity": ("%", "fraction"), "Tm": ("C", "K"),
    "optimal_temperature": ("C", "K"), "residual_activity": ("%", "fraction"),
    "solubility": ("%", "fraction", "mg/mL", "g/L"), "S/H": ("ratio",),
}
ALIASES = {"活性": "activity", "相对活性": "relative_activity", "最适温度": "optimal_temperature",
           "残余活性": "residual_activity", "可溶性": "solubility", "Topt": "optimal_temperature"}
UNITS = {"°C": "C", "℃": "C", "摄氏度": "C", "μmol/min": "umol/min", "µmol/min": "umol/min"}
COLUMNS = list(ExperimentRecord.model_fields)
REQUIRED = ["sample_id", "metric", "value", "unit", "conditions", "method", "replicate_id"]
GROUP = ["sample_id", "candidate_id", "metric", "unit", "conditions", "method", "data_kind"]
COMPARISON = ["metric", "unit", "conditions", "method", "data_kind"]


@dataclass
class EvaluationResult:
    raw: pd.DataFrame
    valid: pd.DataFrame
    errors: pd.DataFrame
    aggregates: pd.DataFrame
    comparisons: pd.DataFrame
    unmatched: pd.DataFrame
    summary: EvaluationSummary


def read_experiments(content: bytes, filename: str) -> pd.DataFrame:
    suffix = Path(filename).suffix.lower()
    if len(content) > 20 * 1024 * 1024:
        raise ValueError("Input exceeds 20 MiB")
    if suffix == ".csv":
        return pd.read_csv(io.BytesIO(content), dtype=object, keep_default_na=False, encoding="utf-8-sig")
    if suffix == ".xlsx":
        return pd.read_excel(io.BytesIO(content), sheet_name=0, dtype=object, keep_default_na=False)
    raise ValueError("Supported formats: .csv, .xlsx (first sheet = measurements)")


def _text(value) -> str:
    return "" if value is None or pd.isna(value) else str(value).strip()


def _condition(value: str) -> str:
    try:
        parsed = json.loads(value)
        if isinstance(parsed, dict):
            return json.dumps(parsed, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    except (ValueError, TypeError):
        pass
    return " ".join(value.split())


def validate_measurements(raw: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    errors: list[dict] = []
    records: list[dict] = []
    missing = sorted(set(REQUIRED) - set(raw.columns))
    if missing:
        return pd.DataFrame(columns=COLUMNS + ["source_row"]), pd.DataFrame(
            [{"row": 1, "column": col, "reason": "missing_required_column"} for col in missing])
    seen: dict[tuple, int] = {}
    duplicates: set[int] = set()
    for offset, (_, row) in enumerate(raw.iterrows(), 2):
        data = {col: _text(row.get(col)) for col in COLUMNS}
        data["metric"] = ALIASES.get(data["metric"], data["metric"])
        data["unit"] = UNITS.get(data["unit"], data["unit"])
        data["conditions"] = _condition(data["conditions"])
        data["data_kind"] = data["data_kind"] or "experimental"
        for col in ("predicted_value", "wase_screen_score"):
            data[col] = data[col] or None
        data["predicted_metric"] = ALIASES.get(data["predicted_metric"], data["predicted_metric"]) or data["metric"]
        data["predicted_unit"] = UNITS.get(data["predicted_unit"], data["predicted_unit"]) or data["unit"]
        try:
            record = ExperimentRecord.model_validate(data)
            if record.unit not in METRIC_UNITS[record.metric]:
                errors.append({"row": offset, "column": "unit", "reason": f"allowed_for_{record.metric}: {', '.join(METRIC_UNITS[record.metric])}"})
                continue
            if (record.unit == "K" and record.value < 0) or (record.unit == "C" and record.value < -273.15):
                errors.append({"row": offset, "column": "value", "reason": "temperature_below_absolute_zero"})
                continue
            if record.value < 0 and record.metric not in {"Tm", "optimal_temperature"}:
                errors.append({"row": offset, "column": "value", "reason": "nonnegative_value_required"})
                continue
            if (record.unit == "%" and record.metric == "solubility" and record.value > 100) or (
                record.unit == "fraction" and record.metric == "solubility" and record.value > 1
            ):
                errors.append({"row": offset, "column": "value", "reason": "fraction_or_percent_out_of_range"})
                continue
            item = record.model_dump()
            key = tuple(item[col] for col in GROUP + ["replicate_id"])
            if key in seen:
                duplicates.update((seen[key], offset))
                errors.append({"row": offset, "column": "replicate_id", "reason": f"duplicate_of_row_{seen[key]}_both_excluded"})
                errors.append({"row": seen[key], "column": "replicate_id", "reason": f"duplicated_by_row_{offset}_both_excluded"})
            else:
                seen[key] = offset
            records.append({**item, "source_row": offset})
        except ValidationError as exc:
            errors.extend({"row": offset, "column": ".".join(map(str, error["loc"])), "reason": error["msg"]}
                          for error in exc.errors())
    table = pd.DataFrame([item for item in records if item["source_row"] not in duplicates], columns=COLUMNS + ["source_row"])
    ambiguous_rows = set()
    for _, group in table.groupby([column for column in GROUP if column != "candidate_id"], sort=False, dropna=False):
        if group.candidate_id.nunique(dropna=False) > 1:
            ambiguous_rows.update(group.source_row)
            errors.extend({"row": int(row), "column": "candidate_id", "reason": "conflicting_candidate_ids_for_same_sample"} for row in group.source_row)
    return table[~table.source_row.isin(ambiguous_rows)], pd.DataFrame(errors, columns=["row", "column", "reason"])


def paired_statistics(measured, predicted) -> dict:
    y, p = np.asarray(measured, dtype=float), np.asarray(predicted, dtype=float)
    if y.ndim != 1 or p.ndim != 1 or len(y) != len(p) or not np.isfinite(y).all() or not np.isfinite(p).all():
        raise ValueError("Paired statistics require finite, equally sized one-dimensional arrays")
    n = len(y)
    result = {"n": n, "mae": None, "rmse": None, "r2": None, "pearson": None, "spearman": None, "statistics_status": "insufficient_samples"}
    if not n:
        return result
    result.update(mae=float(np.mean(np.abs(y - p))), rmse=float(np.sqrt(np.mean((y - p) ** 2))))
    if n < 3:
        return result
    if np.ptp(y) == 0 or np.ptp(p) == 0:
        result["statistics_status"] = "constant_values"
        return result
    result.update(r2=float(1 - np.sum((y - p) ** 2) / np.sum((y - y.mean()) ** 2)),
                  pearson=float(np.corrcoef(y, p)[0, 1]),
                  spearman=float(np.corrcoef(pd.Series(y).rank(), pd.Series(p).rank())[0, 1]),
                  statistics_status="assessed")
    return result


def top_k_hit_rate(frame: pd.DataFrame, score: str, k: int, threshold: float, higher_is_success: bool = True) -> dict:
    if frame.empty:
        return {"top_k_n": 0, "top_k_hit_rate": None, "top_k_status": "no_pairs"}
    rank_ascending = score != "wase_screen_score" and not higher_is_success
    ranked = frame.sort_values([score, "sample_id"], ascending=[rank_ascending, True])
    count = min(k, len(ranked))
    if count < len(ranked) and ranked.iloc[count - 1][score] == ranked.iloc[count][score]:
        return {"top_k_n": count, "top_k_hit_rate": None, "top_k_status": "tie_at_cutoff"}
    values = ranked.head(count)["mean"]
    success = values >= threshold if higher_is_success else values <= threshold
    return {"top_k_n": count, "top_k_hit_rate": float(success.mean()), "top_k_status": "assessed"}


def evaluate_experiments(raw: pd.DataFrame, *, top_k: int = 5, success_threshold: float | None = None,
                         higher_is_success: bool = True) -> EvaluationResult:
    if top_k < 1 or (success_threshold is not None and not np.isfinite(success_threshold)):
        raise ValueError("top_k must be positive and threshold must be finite")
    valid, errors = validate_measurements(raw)
    if success_threshold is not None and not valid.empty and len(valid[COMPARISON].drop_duplicates()) > 1:
        raise ValueError("A success threshold requires one metric/unit/condition/method/data-kind group; select a group first")
    aggregates, comparisons, unmatched = [], [], []
    for key, group in valid.groupby(GROUP, sort=False, dropna=False):
        record = dict(zip(GROUP, key))
        record.update(mean=float(group.value.mean()), std=float(group.value.std(ddof=1)) if len(group) > 1 else None,
                      n=len(group), source_rows=",".join(map(str, group.source_row)))
        record.update(predicted_value=None, prediction_source="", wase_screen_score=None)
        predictions = group[group.predicted_value.notna()]
        reasons = []
        if predictions.empty:
            reasons.append("no_prediction")
        elif predictions.prediction_source.eq("").any():
            reasons.append("prediction_source_missing")
        elif any("wase" in source.lower() for source in predictions.prediction_source):
            reasons.append("wase_is_ranking_only_use_wase_screen_score")
        elif not predictions.predicted_metric.eq(record["metric"]).all() or not predictions.predicted_unit.eq(record["unit"]).all():
            reasons.append("prediction_metric_or_unit_mismatch")
        elif predictions[["predicted_value", "prediction_source"]].drop_duplicates().shape[0] != 1:
            reasons.append("conflicting_replicate_predictions")
        else:
            record.update(predicted_value=float(predictions.iloc[0].predicted_value), prediction_source=predictions.iloc[0].prediction_source)
        scores = group.wase_screen_score.dropna().unique()
        if len(scores) == 1:
            record["wase_screen_score"] = float(scores[0])
        elif len(scores) > 1:
            reasons.append("conflicting_replicate_wase_scores")
        if not record["candidate_id"]:
            reasons.append("candidate_not_linked")
        if reasons:
            unmatched.append({**dict(zip(GROUP, key)), "reason": ";".join(reasons), "source_rows": record["source_rows"]})
        aggregates.append(record)
    agg = pd.DataFrame(aggregates)
    if not agg.empty:
        paired = agg.dropna(subset=["predicted_value"])
        for key, group in paired.groupby(COMPARISON + ["prediction_source"], sort=False, dropna=False):
            record = dict(zip(COMPARISON + ["prediction_source"], key))
            record.update(comparison_type="prediction", **paired_statistics(group["mean"], group.predicted_value))
            if success_threshold is not None:
                record.update(top_k_hit_rate(group, "predicted_value", top_k, success_threshold, higher_is_success))
            comparisons.append(record)
        activity = agg[agg.metric.isin(["activity", "relative_activity", "residual_activity", "S/H"])].dropna(subset=["wase_screen_score"])
        for key, group in activity.groupby(COMPARISON, sort=False, dropna=False):
            stats = paired_statistics(group["mean"], group.wase_screen_score)
            record = dict(zip(COMPARISON, key))
            record.update(comparison_type="wase_ranking_only", prediction_source="WASE-screen", n=len(group),
                          mae=None, rmse=None, r2=None, pearson=None, spearman=stats["spearman"], statistics_status=stats["statistics_status"])
            if success_threshold is not None:
                record.update(top_k_hit_rate(group, "wase_screen_score", top_k, success_threshold, higher_is_success))
            comparisons.append(record)
    summary = EvaluationSummary(raw_rows=len(raw), valid_rows=len(valid), invalid_rows=len(raw) - len(valid),
                                sample_groups=len(agg), comparison_groups=len(comparisons),
                                simulated_rows=int(valid.data_kind.eq("simulated").sum()), top_k=top_k, success_threshold=success_threshold,
                                success_direction="at_least" if higher_is_success else "at_most",
                                notes=["No automatic unit conversion; exact condition and method groups.",
                                       "R2 = 1 - SSE/SST; correlations withheld for n<3 or constant values.",
                                       "WASE-screen is evaluated for ranking only; no success probability.",
                                       "Replicates are not independent samples. Top-k cutoff ties are withheld."])
    return EvaluationResult(raw.copy(), valid, errors, agg, pd.DataFrame(comparisons), pd.DataFrame(unmatched), summary)


def safe_table(frame: pd.DataFrame) -> pd.DataFrame:
    # Prevent spreadsheet formula execution when importing user text into Excel.
    return frame.map(lambda value: "'" + value if isinstance(value, str) and value.startswith(("=", "+", "-", "@")) else value)


def evaluation_excel(result: EvaluationResult) -> bytes:
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        for name in ("raw", "valid", "errors", "aggregates", "comparisons", "unmatched"):
            safe_table(getattr(result, name)).to_excel(writer, sheet_name=name, index=False)
        pd.DataFrame([result.summary.model_dump()]).astype(str).to_excel(writer, sheet_name="summary", index=False)
    return buffer.getvalue()


def template_bytes() -> bytes:
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        pd.DataFrame(columns=COLUMNS).to_excel(writer, sheet_name="measurements", index=False)
        pd.DataFrame([{"metric": name, "units": ", ".join(units)} for name, units in METRIC_UNITS.items()]).to_excel(writer, sheet_name="units", index=False)
    return buffer.getvalue()


def simulated_experiments() -> pd.DataFrame:
    rows = []
    for i, (value, predicted, score) in enumerate(zip([10, 20, 30, 40, 50], [12, 18, 33, 38, 49], [40, 65, 55, 80, 90]), 1):
        for replicate, deviation in enumerate([-1, 1], 1):
            rows.append(ExperimentRecord(sample_id=f"SIM-{i}", candidate_id=f"SIM-C{i}", metric="activity",
                                         value=value + deviation, unit="U/mg", conditions="pH=7;T=30C;substrate=SIM",
                                         method="simulated", replicate_id=str(replicate), predicted_value=predicted,
                                         prediction_source="simulated_predictor", wase_screen_score=score, data_kind="simulated").model_dump())
    return pd.DataFrame(rows)


def save_evaluation(result: EvaluationResult, output_dir: str | Path, source: bytes = b"") -> Path:
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=False)
    for name in ("raw", "valid", "errors", "aggregates", "comparisons", "unmatched"):
        safe_table(getattr(result, name)).to_csv(directory / f"{name}.csv", index=False, encoding="utf-8-sig")
    (directory / "evaluation.xlsx").write_bytes(evaluation_excel(result))
    (directory / "summary.json").write_text(result.summary.model_dump_json(indent=2), encoding="utf-8")
    if source:
        (directory / "source.sha256").write_text(hashlib.sha256(source).hexdigest() + "\n", encoding="ascii")
    return directory
