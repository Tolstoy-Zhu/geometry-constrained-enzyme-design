# AquaEsterNet-StructMine v0.3 File Manifest

## User Entry Points

- `start_ui.bat`: double-click bilingual web launcher with automatic free-port selection.
- `start_here.py`: compatibility launcher for doctor, benchmark, production, validation, and file checks.
- `run_pipeline.py`: direct shared-pipeline runner.
- `prepare_databases.py`: Foldseek database download and preparation.

## Documentation

- `README_FIRST.txt`
- `README_CN.md`
- `README.md`
- `MODEL_GUIDE_CN.md`
- `MODEL_GUIDE_EN.md`
- `PROJECT_OVERVIEW.md`
- `PROJECT_MANIFEST.md`
- `VALIDATION_REPORT.md`
- `COMPLETION_AUDIT.md`
- `CHANGELOG.md`
- `THIRD_PARTY_LICENSES.md`
- `docs/RELEASE_NOTES_CN.md` / `docs/RELEASE_NOTES_EN.md`: bilingual v0.3 release notes.
- `docs/validation/v0.3_validation_summary.json`: sanitized acceptance summary.
- `docs/screenshots/af_prediction_desktop.png` / `af_prediction_mobile.png`: sanitized UI evidence.
- `requirements.lock.txt`: exact verification-environment dependency lock.
- `configs/prediction_runtime.lock.yaml`: isolated WSL prediction-runtime lock and source inventory.
- `configs/colabfold-pip-lock.txt`: exact installed Python/CUDA transitive package versions.
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
- `scripts/run_production_acceptance.ps1`: tests, Foldseek benchmark, Swiss-Prot database, and production run.
- `scripts/prepare_benchmark.py`
- `scripts/run_benchmark.py`
- `scripts/validate_real_search_outputs.py`
- `scripts/write_database_manifest.py`
- `scripts/build_verified_archive.py`
- `scripts/build_public_release.py`: sanitized v0.3 ZIP and SHA-256 builder.
- `scripts/inspect_enabled_sources.py`
- Existing download and database-preparation helpers under `scripts/` remain available.

## Tests

`tests/` covers chain handling, TM-align, overlay generation, filtering, WASE-screen, configuration validation, cache/retry behavior, annotation state, benchmark acceptance, run manifests, cancellation, pinned database validation/extraction, pipeline exports, Streamlit loading, FASTA validation, presets, path conversion, model parsing, MSA consent/cache/retry, GPU locks, and prediction-to-mining lineage. `tests/fixtures/fake_colabfold.py` provides deterministic non-GPU CI execution.

## Runtime Data

The following are generated or downloaded and are not required in a source-only checkout:

- `results/runs/`
- `results/predictions/`
- `results/ui_prediction_inputs/`
- `results/ui_configs/`
- `tmp/`
- `.cache/`
- `databases/`

Validated installations write `databases/environment/foldseek_install_manifest.json`, `afdb_swissprot_download_provenance.json`, and `alphafold_swissprot_database_manifest.json`.

The verified delivery archive may include selected benchmark structures and validation screenshots, but existing delivery archives are preserved.
