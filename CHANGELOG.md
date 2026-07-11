# Changelog

## 0.3.0 - 2026-07-11

### Local AlphaFold2 Prediction

- Added monomer-only LocalColabFold 1.6.1 / AlphaFold2-PTM prediction in the existing `Ubuntu-Foldseek` WSL2 distribution.
- Added strict 20-1000 aa FASTA validation, no more than 5% `X`, and candidate-run/accession input.
- Split public MMseqs2 MSA retrieval from local RTX GPU inference, with sequence-hash caching and three retries.
- Added 15/30-second exponential MSA retry delays and preserved upstream service errors when ColabFold exits successfully without an A3M.
- Added per-job third-party transfer consent, serialized MSA access, one-job GPU locking, stale-lock recovery, cancellation, and explicit CUDA/OOM failures without CPU fallback.
- Added a hash-recorded ColabFold 1.6.1 result-download timeout patch after the upstream 6.02-second read timeout caused completed public MSA jobs to loop during archive retrieval.
- Added quick, standard, and high presets; standard uses five models, three recycles, one seed, and rank_001 Amber relaxation.
- Disabled ColabFold unified memory for 16 GB GPU stability and kept rank_001 Amber relaxation on CPU OpenMM; scientific model parameters are unchanged.
- Added verified AlphaFold parameter download, archive and per-model SHA-256 provenance, and explicit CC BY-NC-SA 4.0 notice.

### Scientific Integration

- Added `PredictionRequest`, `PredictionManifest`, `PredictionModelRecord`, and `PredictionSummary` contracts.
- Added `predict` and `predict-and-mine` CLI actions and candidate re-prediction from an existing run.
- Added global 8HFW alignment for dynamic `Ser111/Asp241/His274` mapping plus explicit manual correction.
- Made `query_structure`, `query_chain`, `parent_prediction_id`, and catalytic reference sites effective in the shared mining pipeline.
- Added bidirectional prediction/run lineage and separate pLDDT, pTM, PAE, model-parameter checksum, event, warning, and failure evidence.

### Interface And Delivery

- Added the bilingual AF Prediction workbench with privacy confirmation, mapping preview, progress, ranking, pLDDT, PAE, and confidence-colored 3D structure.
- Extended `doctor` with WSL GPU, JAX/CUDA, ColabFold, model parameters, MSA endpoint, storage, and GPU-lock checks.
- Added a deterministic fake ColabFold executable for non-GPU CI end-to-end tests.
- Added Python 3.11 GitHub Actions, prediction runtime locking, sanitized publication rules, and updated bilingual documentation.
- Added repeatable Playwright checks at 1440x900 and 390x844 for layout, focus, mapping, pLDDT, PAE, and nonblank py3Dmol pixels.

## 0.2.0 - 2026-07-11

### Scientific Corrections

- Replaced the P00520 false-positive target workflow with a fixed positive/negative benchmark.
- Removed zero-threshold acceptance as evidence of a successful scientific workflow.
- Marked the query as PDB 8HFW chain A, residues 12-294, with input SHA-256.
- Added chain-independent parsing and reporting.
- Replaced residue-number truncation with real TM-align through `tmtools==0.3.0`.
- Stopped fabricating E-value, bit score, and probability for local alignments.

### WASE-screen v0.1

- Added transparent structure, catalytic-core, annotation-confidence, and experimental-readiness scores.
- Added evidence coverage, uncertainty, tier, version, scope, and missing-evidence fields.
- Enforced high-priority gates for score, evidence coverage, bidirectional TM-score, and catalytic core.
- Explicitly separated WASE-screen from future formal WASE and S/H-ratio prediction.

### Reproducibility

- Added Pydantic public contracts and pre-run configuration validation.
- Added versioned run directories, input hashes, configuration snapshots, progress events, artifacts, warnings, failures, and cancellation state.
- Added retrying, cached, batch annotation with explicit error columns.
- Added a WSL-aware Foldseek runtime and environment doctor.

### Interface

- Added a bilingual Streamlit workbench with run, candidate, structure, history, and doctor views.
- Added Plotly candidate evidence plots and py3Dmol TM-align overlays.
- Added double-click `start_ui.bat` launcher with automatic port selection.
- Unified CLI and web execution through the same service layer.

### Validation And Delivery

- Added isolated pytest coverage and a fixed benchmark acceptance report.
- Added WSL/Foldseek setup and full production acceptance scripts.
- Added a pinned official-S3 fallback with resumable transfer, length/SHA-256 checks, safe extraction, and provenance output.
- Added bilingual documentation, dependency locking, third-party notices, and validation screenshots.
- Passed the real Foldseek benchmark and the complete AlphaFold/Swiss-Prot v6 2025_03 production search (999 raw, 9 filtered, 9 annotated).

## 0.1.0

- Initial command-line prototype for structural hit parsing, annotation, and export.
