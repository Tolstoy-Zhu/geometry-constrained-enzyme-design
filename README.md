# AquaEsterNet-StructMine v0.4

**AquaEsterNet-StructMine** is a joint structure-and-sequence mining tool for Bur01-like enzymes. It uses [PDB 8HFW](https://www.rcsb.org/structure/8HFW) chain A, Foldseek/TM-align structure evidence, NCBI NR and MMseqs2 sequence expansion, local AlphaFold2 validation, cross-database annotation, transparent WASE-screen triage, and traceable manifests.

This project is not a black-box predictor of esterification selectivity or S/H ratio. It addresses an earlier decision: which Bur01-like proteins are sufficiently supported to justify expression and enzymology experiments.

## Quick Start

Run **aquaesternet-structmine properties-web** (or start_properties_web.bat on Windows) to score design cohorts and their parents with local SWI, TmProt and EnzRank models. Each completed run includes a portable interactive result page. [Property models and setup](docs/LOCAL_PROPERTY_MODELS_CN.md)

The **Structure** view includes native Bur01-family sequence design. AquaEsterNet selects mutable sites, runs ProteinMPNN or LigandMPNN locally, and attributes reviewed coordinates to an exact design sequence. See the [design and validation guide](docs/NATIVE_ENZYME_DESIGN_CN.md).

Double-click:

```text
start_ui.bat
```

The launcher starts a background watchdog and opens the fixed local address `http://127.0.0.1:8501`. The web app no longer depends on an open terminal, and the watchdog restarts Streamlit after an unexpected exit. Chinese is the default language.

To enable automatic startup at Windows logon, double-click:

```text
install_ui_autostart.bat
```

After installation, the workbench starts in the background when the current Windows user signs in. It is not exposed to the LAN or public internet. Use `stop_ui.bat` to stop it, `ui_status.bat` to inspect it, and `uninstall_ui_autostart.bat` to remove only the scheduled task without deleting project data.

Command-line entry points:

```powershell
python start_here.py --action doctor
python start_here.py --action benchmark
python start_here.py --action production-run
python start_here.py --action validate
python start_here.py --action check-files
python -m bur01_structmine.cli predict --fasta query.fasta --preset standard --confirm-public-msa
python -m bur01_structmine.cli predict-and-mine --fasta query.fasta --preset standard --confirm-public-msa --mining-config configs/alphafold_swissprot.yaml
python -m bur01_structmine.cli nr-demo-build --query-fasta data/query/Bur01m_seed.fasta --email name@example.org --confirm-ncbi
python -m bur01_structmine.cli nr-expand --source-run <run_id> --config configs/nr_local_demo.yaml
python -m bur01_structmine.cli nr-export-job --source-run <run_id> --config configs/nr_full_server.yaml
python -m bur01_structmine.cli nr-import-job --bundle <nr_run_id>_server_result.zip
python -m bur01_structmine.cli nr-novelty --fasta designed.fasta --backend ncbi_remote --email name@example.org --confirm-ncbi
```

`target-run` remains as a compatibility alias for the real Foldseek benchmark. It no longer treats P00520 as a positive hit.

## Purpose

1. Discover proteins structurally related to Bur01 in AlphaFold/Swiss-Prot or PDB.
2. Retain remote candidates when sequence identity is modest but fold, coverage, and catalytic-core evidence are conserved.
3. Normalize Foldseek, TM-align, UniProt, AlphaFold DB, RCSB, and NCBI evidence into one table.
4. Export Excel, TSV, and FASTA files for synthesis, expression, and enzymology planning.
5. Provide transparent structural evidence for future formal WASE model development.

## NR Sequence Expansion And Novelty Audit

v0.4 follows the public VenusMine MMseqs2 search and clustering steps: two iterations, `max-seqs=2500`, `search-type=1`, and 50% sequence-identity clustering. Foldseek now exports target sequences; seeds and hits are deduplicated by accession and sequence SHA-256. Every representative is directly realigned to 8HFW chain A so Ser111, Asp241, and His274 are mapped independently.

The workstation stores only a real query-derived NCBI remote hit subset, explicitly scoped as `nr_remote_snapshot`. Complete NR runs through a credential-free Bash/Slurm/PBS package on a server or HPC. Normal searches never download NR, and imported results must pass scope, version, path-safety, and SHA-256 checks.

NR sequence evidence and structure evidence are separate. Before AF2 plus Foldseek/TM-align validation, every candidate remains:

```text
wase_tier = Not evaluated - structure required
wase_scope = not_evaluated_sequence_only
```

The **NR Mining** view supports local-subset execution, server export/import, evidence filtering, and handoff to AF2. See [`docs/NR_WORKFLOW_EN.md`](docs/NR_WORKFLOW_EN.md). `configs/ncbi_prostt5.yaml` is retained only as a legacy small-custom-FASTA route and is not complete NR.

## AlphaFold2 Prediction

v0.3 added LocalColabFold 1.6.1 and AlphaFold2-PTM to `Ubuntu-Foldseek`. v0.4 also accepts an NR run ID and accession, creating bidirectional NR → AF2 → Foldseek/TM-align lineage.

### Input Contract

- Exactly one FASTA record, 20-1000 aa.
- The canonical 20 amino acids plus at most 5% `X`.
- An annotated candidate can be selected by run ID and accession.
- The current AF2 workflow is monomer-only; complexes, ligands, nucleic acids, and covalent modifications are out of scope.

### Presets

| Preset | Models | Recycles | Seeds | Amber relaxation |
|---|---:|---:|---:|---:|
| Quick | 1 | 3 | 1 | None |
| Standard, default | 5 | 3 | 1 | rank_001 |
| High sampling | 5 | 12 | 3 | rank_001 |

The model type is fixed to `alphafold2_ptm`. Templates are disabled unless explicitly enabled. GPU inference is mandatory: JAX/CUDA or memory failures are reported and never silently fall back to CPU.

### MSA Privacy And Caching

Prediction is split into a public ColabFold MMseqs2 MSA request and local RTX 4080 inference. Model ranking is followed by CPU OpenMM relaxation of rank_001 so relaxation does not compete with resident model weights for GPU memory. MSA requests are serialized and cached by normalized-sequence SHA-256. Every cache miss requires explicit confirmation that the sequence will be sent to the third-party service. Network or service failures are retried up to three times with 15- and 30-second exponential delays, and the original service error is written to events and the failure manifest. The approximately 940 GB local ColabFold database is intentionally not installed.

Inside `Ubuntu-Foldseek`, from the mounted project root, install the pinned runtime and verified parameters with:

```bash
sudo bash scripts/install_colabfold_wsl.sh
sudo bash scripts/install_alphafold_params.sh
```

The archive and all five PTM parameter files are SHA-256 recorded. AlphaFold parameters remain under their bundled CC BY-NC-SA 4.0 license and are excluded from Git and release archives.

### Prediction Contract

Every task writes `results/predictions/<prediction_id>/` with the normalized FASTA, hash, A3M, ranked PDBs, canonical `rank_001.pdb`, pLDDT, pTM, PAE, logs, events, runtime versions, parameter checksums, cache status, and failures. Global alignment maps 8HFW `Ser111/Asp241/His274` to the new sequence and accepts explicit corrections. Mapping or identity failures do not stop structure search; catalytic evidence remains missing and is reported as such.

Prediction and mining manifests link to each other through `prediction_id` and `run_id`.

## Query Structure

| Item | Value |
|---|---|
| PDB | 8HFW |
| Query chain | A |
| File | `data/query/Bur01m.pdb` |
| Residue range | 12-294 |
| Chain length | 283 aa |
| Catalytic residues | Ser111 / Asp241 / His274 |
| SHA-256 | `9f51637496713b18d3ec8f0586a1aef19d3c8f96bd3e48fda9e9af6e88da24fe` |

Every run records this provenance in `run_manifest.json`. Multi-chain structures are never concatenated implicitly. Each chain is searched independently and reported through `query_chain` and `target_chain`.

## Fixed Benchmark

| Role | Controls | Purpose |
|---|---|---|
| Self positive | 8HFW chain A | TM-score and RMSD ceiling |
| Related positives | 1XKL, 6TJ2, Q6RYA0, Q9FVW3 | Bur01-like ranking |
| Length-matched negative | P60174 | Reject length-driven similarity |
| Remote negative | P00520 | Reject the former false positive |

Acceptance criteria:

- 8HFW self-alignment TM-score is greater than 0.98 and RMSD is near zero.
- Every positive ranks above both negatives.
- P00520 does not pass the default structural filter.
- Final acceptance must use Foldseek; the TM-align preview cannot replace it.

Benchmark provenance and checksums are stored in:

```text
data/benchmark/benchmark_manifest.yaml
data/benchmark/download_manifest.json
```

## Search Engines

### Foldseek

Foldseek is used for formal benchmark acceptance and full AlphaFold/Swiss-Prot search. The default WSL layout is:

```text
Distribution: Ubuntu-Foldseek
Windows location: configurable WSL import directory
Foldseek: /opt/foldseek/bin/foldseek
Database root: <project>/databases
```

Install and validate the environment:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\setup_wsl_foldseek.ps1
powershell -ExecutionPolicy Bypass -File scripts\run_production_acceptance.ps1
```

### TM-align Small-Set Preview

`tmtools==0.3.0` runs the real TM-align algorithm for small local sets:

```powershell
python -m bur01_structmine.cli benchmark --engine tmtools
```

It reports real qTM, tTM, RMSD, aligned sequences, and residue mappings. Local TM-align rows leave `evalue`, `bits`, and `prob` empty rather than fabricating Foldseek metrics.

## Web Workbench

The workbench has seven views:

- **Run**: query, database, thresholds, preflight checks, run, cancel, progress events, and errors.
- **NR Mining**: local hit-snapshot execution, complete-NR job export/import, sequence evidence, AF2-ready filtering, and AF2 handoff.
- **AF Prediction**: FASTA upload/paste or candidate input, privacy confirmation, catalytic mapping, prediction or predict-and-mine, model ranking, pLDDT, PAE, and confidence-colored 3D structure.
- **Candidates**: filters for source, TM-score, coverage, WASE tier, reviewed state, and evidence coverage.
- **Structure**: TM-align query-candidate overlay with Ser111/Asp241/His274 and mapped target residues.
- **History**: run status, counts, warnings, errors, and result folders.
- **Doctor**: WSL, Foldseek, AVX2, databases, network, storage, Python packages, and output permissions.

Candidates can be exported as Excel, TSV, or FASTA and can be re-annotated from the same service layer.

The NR view also shows a stage-by-stage trace for seed export, iterative MMseqs2 search, deduplication, 50% clustering, catalytic-site mapping, sequence prescreening, and AF2 queue creation. Public v0.4 UI captures are stored under `docs/screenshots/nr_workflow_*.png`; their validation report explicitly marks the data as a sanitized interface fixture rather than a real NCBI result.

### Persistent local operation

- Stable address: `http://127.0.0.1:8501`.
- Automatic startup after Windows logon, without keeping a browser or console window open.
- Health checks every five seconds and automatic Streamlit restart after failure.
- A one-minute Task Scheduler recovery trigger restarts the watchdog itself when needed.
- State and logs stay under `.cache/ui_service/`; they are excluded from releases and do not overwrite scientific results.
- Access remains local-only. The page is unavailable while Windows is asleep, powered off, or awaiting user logon, and returns after resume/sign-in.

### Online expert demonstration

Double-click `install_online_ui.bat` to create a password-protected HTTPS demonstration tunnel. Experts can operate the workbench from another computer without installing the project. The launcher opens the current online URL and `.cache/remote_ui/SHARE_WITH_EXPERTS.txt`.

The default Quick Tunnel URL may change after restart and is intended for temporary demonstrations rather than formal hosting. A durable hostname requires a Cloudflare named Tunnel, a domain, and an Access allowlist. See `ONLINE_DEMO_GUIDE_EN.md` for controls, security boundaries, and fixed-hostname setup.

## Production Search

Prepare the complete AlphaFold/Swiss-Prot Foldseek database:

```powershell
python prepare_databases.py --config configs/alphafold_swissprot.yaml --only alphafold_swissprot
```

Run production search:

```powershell
python -m bur01_structmine.cli production-run --config configs/alphafold_swissprot.yaml
```

This 32 GB workstation intentionally does not download UniProt50. The full AlphaFold/Swiss-Prot database is the supported production target.

The production preset pins the official AlphaFold/Swiss-Prot `v6 2025_03` archive by byte length and SHA-256. If WSL cannot reach Foldseek's Cloudflare redirect, the preparer resumes from the configured official S3 object, validates it, extracts only safe members, and writes a provenance manifest.

## Versioned Run Outputs

Every run writes to:

```text
results/runs/<run_id>/
```

Main files:

```text
run_manifest.json                 input hashes, software, database versions, status
config.snapshot.yaml              validated configuration snapshot
events.jsonl                      stage-by-stage progress events
*.search.raw.tsv/.xlsx            raw structural hits
*.search.filtered.tsv/.xlsx       structural filter output
*.annotated.tsv/.xlsx             annotations and WASE-screen output
*.candidate_proteins.fasta        candidate protein sequences
*.download_ncbi_cds.sh            NCBI CDS helper script
benchmark_report.json             fixed benchmark acceptance report
```

`results/runs/latest.json` points to the latest run. Failures and cancellation states are explicit.

## Important Fields

Structure:

- `search_engine`, `metric_source`, `alignment_version`
- `query_chain`, `target_chain`
- `qtmscore`, `ttmscore`, `alntmscore`
- `rmsd`, `qcov`, `tcov`, `alnlen`, `fident`
- `aligned_query`, `alignment_match_line`, `aligned_target`

Catalytic core:

- `mapped_catalytic_residues`
- `mapped_catalytic_identity_count`
- `catalytic_ser_his_distance`
- `catalytic_asp_his_distance`
- `catalytic_geometry_contacts`

Annotation:

- `target_uniprot_accession`, `pdb_id`, `ncbi_protein_accession`
- `uniprot_is_reviewed`, `alphafold_global_plddt`
- `final_protein_sequence`, `final_sequence_source`
- `annotation_status`, `annotation_errors`

## WASE-screen v0.1

WASE-screen is transparent candidate triage. It is not the future formal WASE model and does not predict S/H ratio.

```text
structure_fit_score           45 points
catalytic_core_score          30 points
annotation_confidence_score   15 points
experimental_readiness_score  10 points
```

Structural score:

```text
25 * min(qTM, tTM)
+ 10 * sqrt(qcov * tcov)
+ 10 * clamp(1 - RMSD / 5, 0, 1)
```

High priority requires all of the following:

- Total score at least 70.
- Evidence coverage at least 75%.
- `min(qTM,tTM) >= 0.5`.
- Catalytic-core score at least 20.

Outputs include `wase_screen_score`, `wase_evidence_coverage`, `wase_screen_uncertainty`, `wase_tier`, `wase_missing_evidence`, and `wase_version=0.1.0`.

Tunnel specialization, alcohol retention, water exclusion, product escape, fermentation-source priors, and evolutionary evidence remain explicit missing evidence. WASE-screen therefore cannot establish a water-phase esterification mechanism.

## Annotation Reliability

External requests use batch retrieval, retries for 429 and 5xx responses, exponential backoff, a 30-day disk cache, stale-cache fallback, and explicit `annotation_status` / `annotation_errors` columns. Network failures are not silently converted to blanks.

## Validation

```powershell
python -m pytest -q
python -m bur01_structmine.cli doctor --config configs/benchmark_foldseek.yaml
python -m bur01_structmine.cli benchmark --engine foldseek
```

Tests use isolated temporary directories and do not write to formal `results`. See `VALIDATION_REPORT.md` for current evidence.

Retained v0.3 acceptance on 2026-07-11:

- `43 passed`; doctor reported `31 passed, 0 warnings, 0 failed`.
- The standard 8HFW prediction produced five models; rank_001 has mean pLDDT `98.07`, pTM `0.950`, and `min(qTM,tTM)=0.9925` against the experimental structure.
- The traceable Q9FVW3 candidate prediction produced five models; rank_001 has mean pLDDT `83.54`, pTM `0.740`, and `min(qTM,tTM)=0.7815` against AlphaFold DB v6.
- The predicted-8HFW Foldseek benchmark placed the lowest positive at `0.8091`, above the highest negative at `0.2493`, and rejected P00520.
- The predicted query searched all 590,183 Swiss-Prot structures and produced 999 raw hits, 12 default-filter candidates, and 12 valid annotations with zero warnings.
- The AF workbench passed screenshot, horizontal-overflow, keyboard-focus, pLDDT, PAE, and nonblank WebGL checks at `1440x900` and `390x844`.

v0.4 development validation on 2026-07-24:

- `65 passed`; the no-network doctor profile reported `29 passed, 3 expected warnings, 0 failed`.
- MMseqs2 `18.8cc5c` and NCBI BLAST+ `2.17.0` pass real runtime checks; BLAST+ MD5 and SHA-256 are recorded.
- The NR UI passed desktop/mobile screenshot, horizontal-overflow, keyboard-focus, fixture-disclosure, and hidden-service-shutdown checks.
- The real NCBI hit snapshot awaits explicit transfer consent, and complete-NR server acceptance is still pending.
- The development branch may be pushed, but merge, `v0.4.0` tag, and Release remain blocked.

## Scientific Boundaries

- Structural similarity does not prove identical catalytic activity.
- WASE-screen prioritizes candidates; it is not an activity or selectivity predictor.
- AlphaFold hits require interpretation with pLDDT, PAE, annotation quality, and experiments.
- Catalytic-triad mapping cannot replace tunnel, hydration, and dynamics evidence.
- Final candidates require synthesis, expression, purification, substrate profiling, and water-phase esterification assays.

## Main Documents

- `MODEL_GUIDE_CN.md` and `MODEL_GUIDE_EN.md`: purpose, principles, innovation, and scope.
- `PROJECT_OVERVIEW.md`: concise route map.
- `PROJECT_MANIFEST.md`: project file manifest.
- `VALIDATION_REPORT.md`: validation evidence.
- `COMPLETION_AUDIT.md`: requirement audit.
- `CHANGELOG.md`: version history.
- `THIRD_PARTY_LICENSES.md`: third-party software notices.
- `docs/NR_WORKFLOW_CN.md` and `docs/NR_WORKFLOW_EN.md`: NR operation and evidence boundary.
- `docs/V0.4_DEVELOPMENT_NOTES_CN.md` and `docs/V0.4_DEVELOPMENT_NOTES_EN.md`: current development status and open release gates.
- `docs/RELEASE_NOTES_CN.md` and `docs/RELEASE_NOTES_EN.md`: bilingual v0.3 release notes.
