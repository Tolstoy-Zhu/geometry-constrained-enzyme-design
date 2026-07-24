# Changelog

## Unreleased

### v0.4 NR Integration

- Added a separate NCBI NR service layer for Foldseek target-sequence seed export, two-iteration MMseqs2 search, accession and sequence-hash deduplication, 50% identity clustering, direct 8HFW catalytic-site mapping, and AF2-ready queues.
- Added real NCBI remote hit-snapshot construction with contact email, explicit per-request transfer consent, serial polling, bounded retries, disk cache, batched sequence retrieval, and `nr_remote_snapshot` provenance.
- Added credential-free complete-NR Bash, Slurm, and PBS job packages, explicit `mmseqs databases NR` management, resource preflight, result manifests, SHA-256 files, safe ZIP import, and `full_nr` scope validation.
- Added independent NR nearest-neighbor novelty audit. A natural NR sequence's own 100% hit is explicitly not novelty evidence.
- Added `nr-demo-build`, `nr-expand`, `nr-export-job`, `nr-import-job`, and `nr-novelty`; `predict` and `predict-and-mine` now accept `--nr-run`.
- Added the bilingual **NR Mining** web view, background local/server tasks, evidence filters, exports, and NR-to-AF2 handoff.
- Added a stage-by-stage NR processing trace and reproducible 2x-density desktop/mobile UI capture validation with automatic hidden-service shutdown.
- Added MMseqs2 `18.8cc5c` and NCBI BLAST+ `2.17.0` diagnostics and runtime provenance.
- Kept all sequence-only candidates at `WASE not evaluated` until AF2 and Foldseek/TM-align validation.
- Hardened complete-NR imports against Windows/POSIX path traversal, symlinks, duplicate members, archive expansion, unsafe run IDs, checksum mismatch, and inconsistent database manifests.
- Ensured direct 8HFW mapping covers every full-NR cluster representative and that local expansion, novelty, export, and import runs write SHA-256 inventories.
- Added fake-MMseqs2/fake-BLAST CI coverage, full-NR package integrity tests, novelty end-to-end tests, old-manifest compatibility, and bidirectional NR-to-AF2 lineage tests.
- Formal `v0.4.0` merge, tag, and Release remain blocked until a complete-NR server result is imported and passes acceptance.

### Persistent Local Web Service

- Replaced the terminal-bound Windows launcher with a fixed-port background watchdog at `127.0.0.1:8501`.
- Added per-user Windows logon startup, five-second health checks, automatic Streamlit recovery, and a one-minute Task Scheduler recovery trigger for the watchdog itself.
- Added double-click install, start, stop, status, and auto-start removal entry points.
- Kept the service local-only and stored machine-specific state and logs under the ignored `.cache/ui_service/` directory.
- Added bilingual operating documentation and source-level regression checks for the persistent UI assets.

### Authenticated Online Demonstration

- Added a Cloudflare HTTPS tunnel controller so invited experts can operate the local workbench from other computers without opening an inbound router port.
- Added quick demonstration mode plus token-based named-tunnel support for a future fixed hostname.
- Added 600,000-iteration PBKDF2-SHA256 application authentication, password rotation, current-user secret ACLs, and ignored share files.
- Added logon startup, one-minute watchdog recovery, tunnel-process recovery, current-URL tracking, and direct no-proxy health checks.
- Disabled server-side "open results folder" actions while online access is enabled.
- Added bilingual online demonstration guides and automated authentication/asset tests.

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
