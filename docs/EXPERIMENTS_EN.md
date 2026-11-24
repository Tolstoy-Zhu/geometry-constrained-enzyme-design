# Experimental Template and Assessment

Download the template in **Experimental assessment**, or use `data/templates/experiment_template.xlsx`. The first sheet (`measurements`) contains records; `units` lists supported units. Inputs: `.xlsx` or UTF-8 CSV, maximum 20 MiB. Macro-enabled workbooks are not supported.

## Fields

Each row is one measurement. Required: `sample_id`, `metric`, finite `value`, `unit`, `conditions`, `method`, `replicate_id`. Leading zeros in identifiers are preserved.

Optional: `candidate_id`, `predicted_value`, `prediction_source`, `predicted_metric`, `predicted_unit`, `wase_screen_score`, `wase_source_run_id`, `wase_version`, `data_kind`. The WASE source run/version are recorded when linking from the web candidate table. Prediction endpoint/unit default to the measured endpoint/unit. Predictions lacking a source are not paired. `data_kind` is `experimental` (default) or `simulated`; supplied examples are explicitly simulated.

Candidate IDs link to exact accession/target matches. The workbench can supplement missing WASE-screen scores from the selected candidate table; it does not overwrite supplied scores or resolve ambiguous matches automatically. Never enter WASE-screen as an activity prediction.

| Metric | Units |
|---|---|
| activity | U/mL, U/mg, U, kat, kat/kg, umol/min, mmol/L/h |
| relative_activity | %, fraction |
| Tm, optimal_temperature | C, K |
| residual_activity | %, fraction |
| solubility | %, fraction, mg/mL, g/L |
| S/H | ratio |

`fraction=1` means 100%. Celsius glyphs are normalized to C, but units are not numerically converted. Tm and optimum temperature remain different endpoints. Describe substrate, temperature, pH and other assay conditions explicitly; a JSON object is accepted. Conditions are whitespace-normalized, JSON keys sorted, and otherwise compared exactly. The system does not infer equivalence between differently written conditions.

## Analysis

Replicates group by sample, candidate, endpoint, unit, conditions, method and data kind. Report mean, sample SD (`ddof=1`), n and source rows. Single measurements have null SD. Replicates are not independent samples. Duplicate replicate IDs exclude both records.

Prediction comparisons further separate prediction sources. Conflicting predictions across replicates are not averaged. Matching endpoint and unit are required.

- MAE = mean(abs(y−p)); RMSE = sqrt(mean((y−p)²)).
- R² = 1−SSE/SST, not Pearson²; it can be negative.
- Pearson, Spearman and R² are withheld for fewer than three samples or constant values on either side. Spearman uses average ranks for ties. MAE/RMSE for small groups are descriptive, not generalization evidence.
- WASE-screen versus activity, relative/residual activity or S/H reports Spearman and Top-k only. No MAE, RMSE, R² or success probability is calculated for that comparison.
- Set Top-k and an experimental success threshold within one endpoint/unit/condition/method group. Predicted values rank in the success direction; WASE ranks descending. Effective k is capped at sample count. A tie across the kth cutoff withholds the hit rate rather than selecting arbitrary winners.

## Errors and Exports

Diagnostics identify the worksheet/CSV record number (including header), column and reason. Invalid rows remain in `raw`; unmatched predictions have explicit reasons. Excel exports include raw, valid, errors, aggregates, comparisons, unmatched and summary. A separate selected-group export retains threshold settings. Spreadsheet formulas in input text are neutralized during export. Production candidate files are not modified.

```powershell
python -m bur01_structmine.cli evaluate-experiments --input data/examples/simulated/experimental_measurements.csv --output results/evaluations/example_01 --top-k 3 --success-threshold 30
```

Existing output directories are rejected. Invalid input rows still produce review files and exit code 2. A shared threshold is rejected for multiple condition groups; split/select the group first.

`data/examples/simulated/` is separate from real experiments. No historical research is claimed to have been generated or validated by this platform.
