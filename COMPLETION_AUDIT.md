# AquaEsterNet-StructMine v0.4 Completion Audit

Audit date: 2026-07-24

Local implementation status: **IMPLEMENTED AND PASSING**

Publication status: **development branch only; do not merge, tag, or release before complete-NR server acceptance**

## v0.4 NR Integration

| Requirement | Status | Evidence |
|---|---|---|
| Preserve v0.3 data, databases, Obsidian, and archives | PASS | New code is additive; generated NR paths remain ignored |
| Do not restore background startup | PASS | No service installation or scheduled-task change |
| Foldseek target sequences | PASS | `tseq` and `theader` added to formal output |
| Independent NR contracts | PASS | Expansion, candidate, run, summary, novelty request/manifest models |
| Two-iteration VenusMine preset | PASS | Fixed in config, command builder, server scripts, and tests |
| Accession plus sequence-hash deduplication | PASS | Local and server post-processing |
| 50% clustering and direct 8HFW mapping | PASS | Representative membership and three-site status fields |
| Transparent sequence-only ranking | PASS | No new PLM composite score |
| WASE evidence boundary | PASS | Sequence-only candidates cannot receive WASE priority |
| Real local NR snapshot implementation | PASS | Consent, email, serial polling, retry, cache, batched efetch |
| Real local NR snapshot execution | PENDING | Awaiting explicit third-party sequence-transfer confirmation |
| Credential-free server package | PASS | Generic Bash/Slurm/PBS plus checksums and bilingual instructions |
| Safe server result import | PASS | Path traversal, scope, version, missing file, and SHA-256 validation |
| Full NR server execution | PENDING | External server/HPC is required |
| NR novelty audit | PASS | Nearest-neighbor evidence and natural self-hit warning |
| NR-to-AF2 linkage | PASS | Bidirectional manifests and CLI/web handoff |
| Seven-view local web workbench | PASS | NR Mining view added; isolated Streamlit tests pass |
| Automated tests | PASS | 65 isolated tests |
| MMseqs2 runtime | PASS | 18.8cc5c |
| NCBI BLAST+ runtime | PASS | 2.17.0, MD5 and SHA-256 recorded |
| Desktop/mobile NR screenshots | PASS | Five 2x-density captures, overflow/focus checks, fixture disclosure, and verified temporary-service shutdown |
| Formal v0.4.0 release | BLOCKED | Complete-NR imported acceptance is mandatory |

The remainder of this document retains the prior v0.3 acceptance evidence. It must not be interpreted as satisfying the new complete-NR release gate.

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
| Automated tests | PASS | 43 passed, 0 failed |
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

## Persistent And Online UI Extension

| Requirement | Status | Evidence |
|---|---|---|
| Background local service | PASS | Fixed localhost port, hidden watchdog, five-second health loop |
| Windows logon startup | PASS | Per-user scheduled task with logon and one-minute recovery triggers |
| Local process recovery | PASS | Streamlit and watchdog fault injection both recovered automatically |
| Authenticated online access | PASS | HTTPS tunnel plus PBKDF2 application password |
| Other-computer interaction | PASS | Public HTTPS health and Streamlit WebSocket negotiation succeeded |
| Tunnel recovery | PASS | New URL propagated to the ignored share file without changing the password |
| Safe shutdown | PASS | Public access can be stopped without deleting local UI, results, databases, or vault files |
| Optional fixed-hostname upgrade | READY | Token mode implemented; owner can add a Cloudflare domain, hostname, and Access policy when a permanent branded URL is desired |

The immediate expert-demonstration route is operational. A Quick Tunnel is not treated as permanent hosting or an SLA-backed endpoint; a stable hostname remains an external-account configuration step.

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

The v0.4 source tree is intended for development branch `agent/nr-integration-v0.4`. It may be pushed for review only after preserving the existing user-modified benchmark manifest. It must not target a merge before the current v0.3 PR is merged and the complete-NR external acceptance gate passes.

No v0.4 tag or GitHub Release is created in this development state. The prior v0.3 branch, PR, tag, and release instructions remain historical evidence. Any future public asset must exclude databases, real NR candidates, model parameters, WSL state, caches, server configuration, credentials, and logs.

## Scientific Boundary

This audit validates software behavior, reproducibility, structure prediction, and structure-mining acceptance. It does not validate enzyme activity, S/H ratio, substrate-channel behavior, water-phase esterification, or a formal WASE model. Those claims require the missing mechanistic evidence and wet-lab experiments listed by WASE-screen.
