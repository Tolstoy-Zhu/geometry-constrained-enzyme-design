from __future__ import annotations

import io
import math

import pandas as pd
import pytest

from bur01_structmine.experiments import (evaluate_experiments, evaluation_excel, paired_statistics, read_experiments,
                                         save_evaluation, simulated_experiments, template_bytes)


def test_hand_calculated_simulated_statistics_and_replicates():
    result = evaluate_experiments(simulated_experiments(), top_k=3, success_threshold=30)
    assert result.summary.raw_rows == result.summary.valid_rows == 10
    assert result.summary.sample_groups == 5
    assert result.aggregates.n.tolist() == [2] * 5
    assert result.aggregates["mean"].tolist() == [10, 20, 30, 40, 50]
    assert result.aggregates.iloc[0]["std"] == pytest.approx(math.sqrt(2))
    quantitative = result.comparisons.iloc[0]
    assert quantitative.mae == 2
    assert quantitative.rmse == pytest.approx(math.sqrt(4.4))
    assert quantitative.r2 == pytest.approx(0.978)
    ranking = result.comparisons.iloc[1]
    assert ranking.spearman == pytest.approx(0.9)
    assert pd.isna(ranking.r2) and pd.isna(ranking.mae) and pd.isna(ranking.pearson)
    assert ranking.top_k_hit_rate == pytest.approx(2 / 3)
    assert len(result.raw) == 10


def test_r2_is_not_correlation_squared():
    result = paired_statistics([1, 2, 3], [11, 12, 13])
    assert result["pearson"] == pytest.approx(1)
    assert result["r2"] == pytest.approx(-149)


@pytest.mark.parametrize("y,p,status", [([1, 2], [2, 3], "insufficient_samples"), ([1, 1, 1], [1, 2, 3], "constant_values"), ([1, 2, 3], [2, 2, 2], "constant_values")])
def test_uninterpretable_correlations_withheld(y, p, status):
    result = paired_statistics(y, p)
    assert result["r2"] is result["pearson"] is result["spearman"] is None
    assert result["statistics_status"] == status


def test_conditions_units_metrics_methods_and_simulation_never_pool():
    first = simulated_experiments()
    other = first.copy()
    for column, value in [("conditions", "pH=8"), ("unit", "U/mL"), ("method", "other"), ("data_kind", "experimental")]:
        other = first.copy()
        other[column] = value
        result = evaluate_experiments(pd.concat([first, other], ignore_index=True))
        assert result.summary.sample_groups == 10
    with pytest.raises(ValueError, match="one metric"):
        evaluate_experiments(pd.concat([first, other]), success_threshold=30)


def test_import_error_rows_and_raw_are_preserved():
    raw = simulated_experiments().astype(object)
    raw.loc[0, "value"] = "bad"
    raw.loc[1, "value"] = float("inf")
    raw.loc[2, "unit"] = "C"
    result = evaluate_experiments(raw)
    assert set(result.errors.row) == {2, 3, 4}
    assert result.summary.invalid_rows == 3
    assert result.raw.loc[0, "value"] == "bad"
    assert evaluate_experiments(raw.drop(columns="sample_id")).errors.iloc[0].row == 1


def test_duplicate_replicates_both_excluded():
    raw = simulated_experiments()
    result = evaluate_experiments(pd.concat([raw, raw.iloc[:1]], ignore_index=True))
    assert result.summary.valid_rows == 9
    assert result.summary.invalid_rows == 2
    assert {2, 12} == set(result.errors.row)


@pytest.mark.parametrize("field,value,reason", [("predicted_unit", "U/mL", "mismatch"), ("predicted_metric", "Tm", "mismatch"), ("prediction_source", "WASE-screen v0.2", "ranking_only"), ("prediction_source", "", "source_missing")])
def test_prediction_matching_failures_visible(field, value, reason):
    raw = simulated_experiments()
    raw[field] = value
    result = evaluate_experiments(raw)
    assert result.aggregates.predicted_value.isna().all()
    assert result.unmatched.reason.str.contains(reason).all()


def test_conflicting_prediction_and_top_k_ties():
    raw = simulated_experiments()
    raw.loc[0, "predicted_value"] = 999
    raw["wase_screen_score"] = 50
    result = evaluate_experiments(raw, top_k=2, success_threshold=30)
    assert "conflicting_replicate_predictions" in result.unmatched.iloc[0].reason
    assert result.comparisons.iloc[-1].top_k_status == "tie_at_cutoff"
    assert pd.isna(result.comparisons.iloc[-1].top_k_hit_rate)


def test_excel_csv_roundtrip_and_exclusive_export(tmp_path):
    raw = simulated_experiments()
    raw["sample_id"] = ["001"] * 2 + raw.sample_id.tolist()[2:]
    data = raw.to_csv(index=False).encode("utf-8-sig")
    assert read_experiments(data, "test.csv").sample_id.iloc[0] == "001"
    result = evaluate_experiments(raw)
    content = evaluation_excel(result)
    assert read_experiments(content, "test.xlsx").sample_id.iloc[0] == "001"
    assert list(pd.read_excel(io.BytesIO(template_bytes())).columns) == list(raw.columns)
    path = save_evaluation(result, tmp_path / "output", data)
    assert (path / "source.sha256").is_file()
    with pytest.raises(FileExistsError):
        save_evaluation(result, path)


def test_spreadsheet_formula_neutralized():
    raw = simulated_experiments()
    raw["method"] = '=HYPERLINK("https://example.invalid")'
    result = evaluate_experiments(raw)
    exported = pd.read_excel(io.BytesIO(evaluation_excel(result)), sheet_name="raw")
    assert exported.method.iloc[0].startswith("'=")


def test_top_k_lower_success_direction_is_recorded():
    result = evaluate_experiments(simulated_experiments(), top_k=2, success_threshold=20, higher_is_success=False)
    assert result.comparisons.iloc[0].top_k_hit_rate == 1
    assert result.summary.success_direction == "at_most"


def test_ambiguous_candidate_links_do_not_inflate_samples():
    raw = simulated_experiments()
    raw.loc[0, "candidate_id"] = "different_candidate"
    result = evaluate_experiments(raw)
    assert result.summary.sample_groups == 4
    assert result.summary.invalid_rows == 2
    assert result.errors.reason.eq("conflicting_candidate_ids_for_same_sample").all()


def test_residual_activity_above_100_can_be_measured():
    raw = simulated_experiments().iloc[:1].copy()
    raw["metric"], raw["unit"], raw["value"] = "residual_activity", "%", 110
    assert evaluate_experiments(raw).summary.valid_rows == 1
