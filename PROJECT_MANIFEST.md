# AquaEsterNet-StructMine v0.4 File Manifest

## User Entry Points

- `start_ui.bat`: starts or reconnects to the persistent local UI at `127.0.0.1:8501` and opens the browser.
- `install_ui_autostart.bat`: installs the current-user Windows logon task and starts the UI.
- `stop_ui.bat`: stops the watchdog and Streamlit process without removing data.
- `ui_status.bat`: displays health, task state, process IDs, restart count, and log path.
- `uninstall_ui_autostart.bat`: stops the UI and removes only the scheduled task.
- `install_online_ui.bat` / `start_online_ui.bat`: install or start authenticated HTTPS expert access.
- `online_ui_status.bat` / `stop_online_ui.bat`: inspect or stop public access.
- `reset_online_password.bat`: rotate the shared expert password and rebuild the tunnel.
- `uninstall_online_ui.bat`: remove the remote task and local access secrets without deleting project data.
- `start_here.py`: compatibility launcher for doctor, benchmark, production, validation, and file checks.
- `run_pipeline.py`: direct shared-pipeline runner.
- `prepare_databases.py`: Foldseek database download and preparation.

## Documentation

- `README_FIRST.txt`
- `README_CN.md`
- `README.md`
- `MODEL_GUIDE_CN.md`
- `MODEL_GUIDE_EN.md`
- `ONLINE_DEMO_GUIDE_CN.md`
- `ONLINE_DEMO_GUIDE_EN.md`
- `PROJECT_OVERVIEW.md`
- `PROJECT_MANIFEST.md`
- `VALIDATION_REPORT.md`
- `COMPLETION_AUDIT.md`
- `CHANGELOG.md`
- `THIRD_PARTY_LICENSES.md`
- `docs/RELEASE_NOTES_CN.md` / `docs/RELEASE_NOTES_EN.md`: bilingual v0.3 release notes.
- `docs/NR_WORKFLOW_CN.md` / `docs/NR_WORKFLOW_EN.md`: bilingual v0.4 NR local/server workflow and evidence boundaries.
- `docs/validation/v0.4_development_validation.json`: sanitized v0.4 development-gate evidence.
- `docs/validation/v0.4_nr_ui_validation.json`: NR desktop/mobile screenshot, layout, focus, checksum, and temporary-service evidence.
- `docs/validation/v0.3_validation_summary.json`: sanitized acceptance summary.
- `docs/validation/persistent_online_ui_validation.json`: sanitized local/online service recovery and security evidence.
- `docs/screenshots/af_prediction_desktop.png` / `af_prediction_mobile.png`: sanitized UI evidence.
- `docs/screenshots/nr_workflow_01_workbench.png` through `nr_workflow_05_mobile.png`: 2x-density sanitized NR workbench evidence.
- `requirements.lock.txt`: exact verification-environment dependency lock.
- `configs/prediction_runtime.lock.yaml`: isolated WSL prediction-runtime lock and source inventory.
- `configs/colabfold-pip-lock.txt`: exact installed Python/CUDA transitive package versions.
- `configs/cloudflared_runtime.lock.json`: verified Windows tunnel binary source, checksum, license, and deployment scope.
- `.github/workflows/ci.yml`: Python 3.11 non-GPU validation.

## Configuration

- `config.yaml`: default PDB Foldseek route.
- `configs/benchmark_foldseek.yaml`: formal fixed benchmark.
- `configs/benchmark_tmtools.yaml`: fixed small-set TM-align preview.
- `configs/alphafold_swissprot.yaml`: production AlphaFold/Swiss-Prot route.
  It pins the official v6 2025_03 archive URL, byte length, and SHA-256.
- `configs/prediction_runtime.lock.yaml`: LocalColabFold, JAX, OpenMM, MMseqs2, HHsuite, Kalign, and preset pins.
- `configs/alphafold_local_structures.yaml`: local AlphaFold folder through Foldseek.
- `configs/alphafold_local_python.yaml`: local AlphaFold folder through TM-align.
- `configs/swissprot_local_python.yaml`: reviewed local structures through TM-align.
- `configs/ncbi_prostt5.yaml`: NCBI FASTA through ProstT5/Foldseek.
- `configs/nr_local_demo.yaml`: query-derived NCBI NR hit-snapshot preset.
- `configs/nr_full_server.yaml`: complete-NR server/HPC preset.
- `configs/nr_runtime.lock.yaml`: pinned MMseqs2 and NCBI BLAST+ sources and checksums.
- `configs/ncbi_predicted_structures_python.yaml`: NCBI predicted structures through TM-align.
- `configs/multidb_template.yaml`: multi-source template.

## Included Data

- `data/query/Bur01m.pdb`: PDB 8HFW chain A query, residues 12-294.
- `data/query/Bur01m_seed.fasta`: query sequence.
- `data/example_foldseek_result.tsv`: parser and pipeline fixture.
- `data/benchmark/benchmark_manifest.yaml`: fixed benchmark definition.
- `data/benchmark/download_manifest.json`: benchmark source URLs and hashes.
- `data/benchmark/structures/`: seven fixed benchmark controls.
- `data/benchmark/README.md`: benchmark documentation.
- `data/alphafold/README.md`: AlphaFold data notes.
- `data/ncbi/README.md`: NCBI data notes.

## Web Interface

- `ui/app.py`: bilingual Streamlit workbench.
- `.streamlit/config.toml`: fixed color theme and local server defaults.
- `src/bur01_structmine/visualization.py`: TM-align superposition for py3Dmol.
- `scripts/windows/ui_service.ps1`: local-only watchdog, health check, automatic recovery, state/log handling, and Windows logon-task controller.
- `scripts/windows/remote_ui_service.ps1`: authenticated Cloudflare HTTPS tunnel, quick/named modes, password rotation, watchdog, and online state controller.
- `src/bur01_structmine/remote_auth.py`: bounded PBKDF2-SHA256 access-record loading and constant-time password verification.

## Core Package

- `src/bur01_structmine/models.py`: Pydantic public contracts.
- `src/bur01_structmine/prediction.py`: shared prediction and prediction-to-mining service.
- `src/bur01_structmine/prediction_sequence.py`: monomer FASTA validation and 8HFW catalytic mapping.
- `src/bur01_structmine/colabfold_runtime.py`: command construction, MSA/GPU locks, cancellation, and model parsing.
- `src/bur01_structmine/config.py`: validated configuration loading.
- `src/bur01_structmine/pipeline.py`: shared search, filter, annotation, scoring, and export service.
- `src/bur01_structmine/run_context.py`: run manifests and progress events.
- `src/bur01_structmine/foldseek_runner.py`: native and WSL Foldseek runtime.
  It also implements resumable pinned downloads and safe archive extraction.
- `src/bur01_structmine/structure_search.py`: chain-aware TM-align preview.
- `src/bur01_structmine/results_parser.py`: Foldseek parsing and structural filters.
- `src/bur01_structmine/scoring.py`: WASE-screen v0.1.
- `src/bur01_structmine/mmseqs_runner.py`: native/WSL MMseqs2 commands and execution.
- `src/bur01_structmine/ncbi_blast.py`: consent-gated, serial, cached NCBI remote BLAST and sequence retrieval.
- `src/bur01_structmine/nr_pipeline.py`: seed export, search, accession/sequence deduplication, clustering, catalytic mapping, sequence screen, novelty audit, and manifests.
- `src/bur01_structmine/nr_job_bundle.py`: credential-free complete-NR export and verified import.
- `src/bur01_structmine/benchmark.py`: fixed benchmark preparation and acceptance.
- `src/bur01_structmine/doctor.py`: environment diagnostics.
- `src/bur01_structmine/annotator.py`: cross-database annotation.
- `src/bur01_structmine/http_utils.py`: retries, backoff, and disk cache.
- `src/bur01_structmine/*_client.py`: RCSB, UniProt, AlphaFold, and NCBI clients.
- `src/bur01_structmine/cli.py`: unified command-line interface.

## Setup And Acceptance Scripts

- `scripts/setup_wsl_foldseek.ps1`: WSL2 and Foldseek installation under D:.
- `scripts/install_colabfold_wsl.sh`: pinned LocalColabFold 1.6.1 / CUDA 12 environment under `/opt/colabfold`.
- `scripts/install_alphafold_params.sh`: resumable verified AlphaFold2 parameter installer.
- `scripts/run_af2_acceptance.py`: real 8HFW/Q9FVW3 and linked-mining acceptance report.
- `scripts/validate_ui_playwright.py`: desktop/mobile interaction, focus, chart, and WebGL pixel acceptance.
- `scripts/validate_nr_ui_playwright.py`: isolated NR fixture, hidden temporary UI, desktop/mobile capture, overflow/focus checks, and verified shutdown.
- `scripts/run_production_acceptance.ps1`: tests, Foldseek benchmark, Swiss-Prot database, and production run.
- `scripts/prepare_benchmark.py`
- `scripts/run_benchmark.py`
- `scripts/validate_real_search_outputs.py`
- `scripts/write_database_manifest.py`
- `scripts/build_verified_archive.py`
- `scripts/build_public_release.py`: sanitized v0.3 ZIP and SHA-256 builder.
- `scripts/install_ncbi_blast_wsl.sh`: pinned and checksum-verified NCBI BLAST+ installer.
- `scripts/nr_server_postprocess.py`: dependency-free server deduplication, manifest, checksum, and archive helper.
- `scripts/inspect_enabled_sources.py`
- Existing download and database-preparation helpers under `scripts/` remain available.

## Tests

`tests/` covers chain handling, TM-align, overlay generation, filtering, WASE-screen, configuration validation, cache/retry behavior, annotation state, benchmark acceptance, run manifests, cancellation, pinned database validation/extraction, pipeline exports, Streamlit loading, FASTA validation, presets, path conversion, model parsing, MSA consent/cache/retry, GPU locks, MMseqs2 command construction, NCBI consent and rate limiting, accession/sequence deduplication, clustering, catalytic mapping, complete-NR bundle integrity, safe import, and NR-to-AF2 lineage. `tests/fixtures/fake_colabfold.py` and `tests/fixtures/fake_mmseqs.py` provide deterministic non-GPU/non-NR CI execution.

## Runtime Data

The following are generated or downloaded and are not required in a source-only checkout:

- `results/runs/`
- `results/predictions/`
- `results/nr_runs/`
- `results/nr_novelty/`
- `results/ui_prediction_inputs/`
- `results/ui_configs/`
- `tmp/`
- `.cache/`
- `databases/`

Validated installations write `databases/environment/foldseek_install_manifest.json`, `afdb_swissprot_download_provenance.json`, and `alphafold_swissprot_database_manifest.json`.

The verified delivery archive may include selected benchmark structures and validation screenshots, but existing delivery archives are preserved.
