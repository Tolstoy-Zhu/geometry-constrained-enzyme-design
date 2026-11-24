# WASE-screen 0.2.0

Application version: 0.5.0.dev0. This transparent triage score does not predict activity, S/H or experimental success probability.

## Components

| Component | Points | Rule |
|---|---:|---|
| Structure fit | 45 | 25×min(qTM,tTM) + 10×sqrt(qcov×tcov) + 10×clamp(1−RMSD/5,0,1) |
| Catalytic core | 30 | Mapped Ser-centered GXSXG:5; mapped S/D/H identities:5 each; two geometry distances:5 each |
| Annotation confidence | 15 | Relevant function:5; reviewed/experimental evidence:5; structure confidence:5 |
| Experimental readiness | 10 | Available sequence without X:4; 180–500 aa:2; tcov≥0.8:2; traceable identifier:2 |

TM and coverage inputs are clamped to [0,1] for scoring, without changing the original search metrics. Approximate Foldseek `alntmscore` may slightly exceed 1; query- and target-normalized TM values retain their separate definitions.

High priority requires score≥70, assessable field weight≥75, min(qTM,tTM)≥0.5 and catalytic score≥20. Score≥50 is medium unless high-priority criteria pass; assessable weight<60 is first classified as incomplete evidence. Structural filtering remains independent and uses each run's thresholds. High priority does not automatically imply passing every structural filter.

## Motif and Geometry

TM-align maps reference sites to the target chain. Only the five-residue continuous coordinate-sequence window centered on mapped Ser can earn motif points. Missing flanks or unresolved gaps are missing evidence. Whole-sequence `has_GXSXG_like_motif` is auxiliary, with no catalytic credit.

Relevant Ser OG, Asp OD1/OD2 and His ND1/NE2 atoms are checked before taking the minimum Ser-His or Asp-His distance. Identity mismatches remain explicit. Each fully assessable distance earns 5 points within 2–5 Angstrom. These are proximity criteria, not confirmed hydrogen bonds or proton-transfer evidence.

- `compatible`: both distances assessed and in range.
- `out_of_range`: both assessed, at least one outside range.
- `not_assessable`: not fully assessed; total contact count is null. An individually assessable pair can still contribute 5 points and its field weight; the other remains missing.

Supplemental mapping preserves Foldseek metrics and records a TM-align mapping source and structure SHA-256. Rescoring uses local coordinates only, without downloads or Foldseek execution.

## Missingness and Error

`wase_information_missingness = 100 − wase_evidence_coverage`. Coverage measures the assessable weight of the 100-point rubric, not completeness of all biological evidence. Tunnel, hydration, alcohol residence, substrate, ecological and evolutionary gaps remain listed even at 100% rubric coverage.

New `wase_screen_uncertainty` and `wase_prediction_error` are null; `wase_prediction_error_status=not_estimated`. Legacy uncertainty is displayed only as information missingness. There is no probability calibration or predictive-error interval.

## Independent Rescoring

Use **Save independent rescore** or:

```powershell
python -m bur01_structmine.cli rescore --run-dir results/runs/<run_id>
```

New TSV, Excel and `rescore_manifest.json` are stored in `results/rescores/<id>/`. They retain source run ID, source hashes, previous version/score and score delta. Original tables and manifests are untouched. Reannotation also creates a new run. Historical predicted queries without reference-site mappings remain unmapped; 8HFW residue numbering is never assumed for them.
