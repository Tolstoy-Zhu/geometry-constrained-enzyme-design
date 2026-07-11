# AquaEsterNet-StructMine v0.3 Completion Audit

Audit date: 2026-07-11

Local implementation status: **COMPLETE AND PASSING**

Publication status: **ready for non-draft PR; merge and Release remain user-controlled checkpoints**

## Preservation

| Requirement | Status | Evidence |
|---|---|---|
| Preserve v0.2 results | PASS | Existing `results/runs/` content retained |
| Preserve complete Foldseek database | PASS | AlphaFold/Swiss-Prot v6 2025_03 remains installed and validated |
| Preserve Obsidian vault | PASS | `Github-Bur01` excluded from code and package operations |
| Preserve earlier delivery archives | PASS | Existing archives are not overwritten |
| Create pre-v0.3 checkpoint | PASS | Separate snapshot created before implementation |

## Prediction Runtime

| Requirement | Status | Evidence |
|---|---|---|
| LocalColabFold 1.6.1 in existing WSL | PASS | Isolated `/opt/colabfold` runtime |
| Python 3.13 and CUDA/JAX GPU | PASS | Python 3.13.14, JAX 0.5.3, RTX 4080 detected |
| Fixed dependency versions | PASS | WSL package lock, runtime lock, install provenance |
| AlphaFold2-PTM parameters | PASS | Archive plus five model SHA-256 values recorded |
| Parameter licensing | PASS | Bundled CC BY-NC-SA 4.0 identified; weights excluded from delivery |
| No Docker or local 940 GB MSA database | PASS | Public MSA plus local inference only |
| No CPU inference fallback | PASS | Missing GPU and OOM are explicit failures |
| Stable 16 GB allocation | PASS | Unified memory disabled; model parameters unchanged |
| rank_001 Amber relaxation | PASS | CPU OpenMM after pLDDT ranking |

## Input, Privacy, And Resources

| Requirement | Status | Evidence |
|---|---|---|
| Monomer-only 20-1000 aa FASTA | PASS | Aggregated Pydantic/sequence validation |
| Canonical residues and <=5% X | PASS | Unit and end-to-end tests |
| FASTA or candidate input | PASS | CLI and web use shared resolver |
| Per-job public-MSA consent | PASS | Cache miss blocked without confirmation |
| Sequence-hash MSA cache | PASS | Real 8HFW and Q9FVW3 cache records |
| Serialized public MSA | PASS | Process lock and stale-owner recovery |
| Three retries with backoff | PASS | 15/30-second delays and upstream error retention |
| One prediction per GPU | PASS | GPU process lock |
| Cancel and stale lock recovery | PASS | Unit tests and explicit cancel marker |

## Scientific Integration

| Requirement | Status | Evidence |
|---|---|---|
| Prediction Pydantic contracts | PASS | Request, manifest, model record, summary |
| `predict` CLI | PASS | FASTA and candidate routes |
| `predict-and-mine` CLI | PASS | Shared prediction and mining service |
| Dynamic query structure and chain | PASS | Effective in RunRequest and pipeline |
| Parent prediction linkage | PASS | Bidirectional prediction/run manifests |
| Dynamic catalytic sites | PASS | 8HFW global mapping and manual override |
| Continue when mapping is incomplete | PASS | Search continues with missing catalytic evidence |
| Five ranked model outputs | PASS | Real 8HFW and Q9FVW3 runs |
| pLDDT, pTM, PAE, logs, events | PASS | Prediction artifact contract |

## Real Acceptance

| Check | Status | Observed |
|---|---|---|
| Automated tests | PASS | 34 passed, 0 failed |
| Doctor | PASS | 31 passed, 0 warnings, 0 failed |
| 8HFW prediction | PASS | pLDDT 98.07, pTM 0.950, min TM 0.9925 |
| Q9FVW3 prediction | PASS | pLDDT 83.54, pTM 0.740, min TM 0.7815 |
| Q9FVW3 lineage | PASS | Source run and accession recorded |
| Predicted-query benchmark | PASS | All positives above negatives; P00520 rejected |
| Complete Swiss-Prot search | PASS | 999 raw, 12 filtered, 12 annotated, no warnings |
| Silent-error rejection | PASS | Semantic MSA failure and postprocessing failure did not become success |

## Web Workbench

| Requirement | Status | Evidence |
|---|---|---|
| Six bilingual views | PASS | Run, AF Prediction, Candidates, Structure, History, Doctor |
| FASTA upload/paste/candidate | PASS | Shared validation and candidate resolver |
| Presets, template switch, consent | PASS | Feature-complete controls |
| Catalytic preview and correction | PASS | Mapping table and manual inputs |
| Predict, predict-and-mine, cancel | PASS | Shared subprocess service |
| Ranking, pLDDT, PAE, 3D structure | PASS | Completed prediction display |
| Candidate re-prediction | PASS | Candidate-table handoff |
| 1440x900 acceptance | PASS | No overflow; focus, charts, and WebGL nonblank |
| 390x844 acceptance | PASS | No overflow; focus, charts, and WebGL nonblank |
| Specified visual palette | PASS | Structure-green theme, amber/red evidence colors, IBM Plex |

## CI, Documentation, And Packaging

| Requirement | Status | Evidence |
|---|---|---|
| Non-GPU GitHub Actions | PASS | Python 3.11, fake ColabFold, tests, file and dependency checks |
| Bilingual README | PASS | `README_CN.md` and `README.md` synchronized |
| Bilingual model guide | PASS | Purpose, principle, innovation, use, and boundaries |
| Bilingual release notes | PASS | `docs/RELEASE_NOTES_CN.md` and `docs/RELEASE_NOTES_EN.md` |
| Runtime and third-party notices | PASS | Exact pins, parameter provenance, licenses |
| Sanitized validation summary | PASS | No machine paths, logs, or candidate sequences |
| Sanitized screenshots | PASS | Public 8HFW fixture only |
| Public archive guard | PASS | Excludes results, databases, weights, caches, vault, logs, and local paths |
| Preserve prior archives | PASS | New v0.3 archive receives a new name |

## GitHub Publication Gate

The audited source tree is intended for branch `agent/af-prediction-v0.3` with commit message `Add local AlphaFold2 prediction workflow` and a ready, non-draft PR titled `Add local AlphaFold2 prediction and AquaEsterNet-StructMine v0.3`.

The PR is not auto-merged. Tag `v0.3.0` and the GitHub Release are created only after the user merges the PR. Release assets must be the sanitized ZIP and its independent SHA-256 file; databases, candidate sequences, model parameters, WSL state, caches, and logs remain excluded.

## Scientific Boundary

This audit validates software behavior, reproducibility, structure prediction, and structure-mining acceptance. It does not validate enzyme activity, S/H ratio, substrate-channel behavior, water-phase esterification, or a formal WASE model. Those claims require the missing mechanistic evidence and wet-lab experiments listed by WASE-screen.
