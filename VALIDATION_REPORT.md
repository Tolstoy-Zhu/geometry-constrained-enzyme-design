# AquaEsterNet-StructMine v0.4 Development Validation Report

Validation date: 2026-07-24

Project version: 0.4.0

Development result: **PASS**

Formal release result: **BLOCKED PENDING COMPLETE-NR SERVER ACCEPTANCE**

## v0.4 NR Development Gate

| Requirement | Status | Evidence |
|---|---|---|
| Existing 43-test baseline remains stable | PASS | Full isolated suite completes without WSL/doctor hang |
| Expanded automated suite | PASS | `65 passed in 15.97s` |
| MMseqs2 runtime | PASS | `18.8cc5c` at `/opt/colabfold/env/bin/mmseqs` |
| NCBI BLAST+ runtime | PASS | `2.17.0` at `/opt/ncbi-blast/bin/blastp` |
| BLAST+ official archive MD5 | PASS | `bdec166721de3b55f90a3badc83538e8` |
| BLAST+ locally calculated SHA-256 | PASS | `3888112d8207831aa47371d93583c601f058f88b5db22dc782438b039a3a411b` |
| Sequence/structure evidence separation | PASS | NR candidates force `wase_scope=not_evaluated_sequence_only` |
| MMseqs2 paper preset | PASS | 2 iterations, max-seqs 2500, search-type 1, 50% clustering |
| Accession and sequence deduplication | PASS | Best-hit selection plus complete-sequence SHA-256 |
| Catalytic mapping | PASS | Separate Ser111/Asp241/His274 positions, identities, and missing reasons |
| P00520 default sequence screen | PASS in fixed CI fixture | Negative remains `hold` and is not AF2-ready |
| Complete-NR job export/import | PASS in deterministic fixture | Bash/Slurm/PBS, resource check, safe ZIP, manifest, SHA-256 |
| NR to AF2 lineage | PASS in deterministic fixture | `source_nr_run_id` and `child_prediction_ids` are bidirectional |
| Streamlit isolation | PASS | Seven-view UI and remote-auth tests finish without real WSL calls |
| NR desktop/mobile UI | PASS | Four `2880x1800` desktop captures and one `780x1688` mobile capture; no horizontal overflow; keyboard focus visible |
| Screenshot service isolation | PASS | Hidden localhost-only Streamlit process stopped after capture; no devtunnel started |
| Real local NCBI NR hit snapshot | PENDING CONSENT | Contact email and explicit sequence-transfer confirmation required |
| At least one real NR candidate through AF2 and structure validation | PENDING | Requires the real local snapshot |
| Complete-NR server/HPC import | PENDING | External server resources and returned result required |
| v0.4.0 tag and Release | BLOCKED | Allowed only after complete-NR acceptance |

The no-network doctor run against the TM-align benchmark profile reports 29 passes, three expected warnings, and zero failures. The warnings are the missing real NR snapshot, no exported full-NR job yet, and no `AQUA_NCBI_EMAIL`; they are explicit gates rather than silent success.

The NR screenshots use a sanitized deterministic interface fixture and say so visibly. They demonstrate the workbench, processing trace, sequence evidence, complete-NR server package controls, and responsive mobile layout without claiming that a real NCBI query has already run. Image hashes and layout checks are recorded in `docs/validation/v0.4_nr_ui_validation.json`.

## Retained v0.3 Acceptance

## Scope

The retained section below covers the previously accepted v0.3 monomer AlphaFold2 workflow, dynamic prediction-to-mining integration, real Foldseek acceptance, complete AlphaFold/Swiss-Prot search, bilingual UI, non-GPU CI path, reproducibility records, and public-package safeguards.

Existing v0.2 results, databases, the Obsidian vault, and prior delivery archives were preserved. A separate pre-v0.3 snapshot was created before implementation.

## Environment

| Component | Verified value |
|---|---|
| Windows Python | 3.11.9 |
| WSL distribution | Ubuntu-Foldseek |
| GPU | NVIDIA GeForce RTX 4080, 16 GB |
| LocalColabFold | 1.6.1 |
| AlphaFold-ColabFold | 2.3.13 |
| JAX | 0.5.3, CUDA device detected |
| OpenMM | 8.5.2, CUDA platform available |
| MMseqs2 | 18.8cc5c |
| NCBI BLAST+ | 2.17.0 |
| HHsuite | 3.3.0 |
| Kalign | 2.04 |
| Foldseek | `7f4b94644d6b210af4530e04a9fdc9e4bb64ea66` |
| AlphaFold parameter archive SHA-256 | `ad158bc5e80eb6fa339928c2ad9a052e5386bb1a650acf3c39b4fb57da4e369d` |
| PTM parameter checksums | Five recorded and present |

The ColabFold result-archive read timeout patch is provenance-recorded with original and patched source hashes and declares no algorithm change. JAX unified memory is disabled for stable 16 GB allocation; model type, weights, recycles, seeds, MSA, and ranking are unchanged.

## Automated Checks

```text
65 passed in 15.97s
v0.4 no-network doctor: 29 passed, 3 expected warnings, 0 failed
```

Tests use isolated temporary directories. The suite covers sequence validation, presets, command construction, Windows/WSL path conversion, cache and retry behavior, semantic MSA failures, process locks, cancellation, parameter provenance, result parsing, catalytic mapping, manifests, prediction-to-mining lineage, chain parsing, TM-align, Foldseek filtering, annotation, WASE-screen, configuration, output contracts, MMseqs2, NCBI consent/rate limits/cache, NR deduplication and clustering, complete-NR package integrity, safe import, and NR-to-AF2 lineage.

The current suite additionally covers persistent Windows launch assets, local-only origin binding, remote access-record parsing, bounded PBKDF2 password verification, authenticated Streamlit gate behavior, online entry points, and secret-file exclusion.

## Persistent And Online UI Extension

| Requirement | Status | Evidence |
|---|---|---|
| Stable local address | PASS | `127.0.0.1:8501`, health response `ok` |
| Logon auto-start | PASS | Current-user Task Scheduler logon trigger |
| Streamlit recovery | PASS | Forced process termination recovered in 10.9 s with a new PID |
| Local watchdog recovery | PASS | One-minute recovery trigger restored the watchdog in 51.3 s |
| HTTPS expert access | PASS | Cloudflare Tunnel health endpoint reachable without an inbound port |
| Real-time protocol | PASS | Public `wss://.../_stcore/stream` negotiated the `streamlit` subprotocol |
| Application authentication | PASS | Random-salt PBKDF2-SHA256, 600,000 iterations, constant-time comparison |
| Tunnel recovery | PASS | Forced termination recovered in 26.1 s; share URL updated and password preserved |
| Remote watchdog recovery | PASS | One-minute recovery trigger restored the watchdog in 53.4 s without dropping the tunnel |
| Secret exclusion | PASS | Access/share/state files matched the ignored `.cache/` rule |
| Remote host-side action safety | PASS | "Open results folder" disabled while online access is enabled |

Machine-specific URLs, passwords, process IDs, logs, and task state are excluded from the source tree. Sanitized evidence is stored in `docs/validation/persistent_online_ui_validation.json`.

GitHub Actions uses Python 3.11 with the deterministic fake ColabFold executable. CI does not install GPU drivers, download model weights, or call the public MSA service.

## Real 8HFW Prediction

Input: public PDB 8HFW chain A sequence, 283 aa.

Preset: standard, five AlphaFold2-PTM models, three recycles, one seed, rank_001 Amber relaxation, templates off.

| Check | Observed | Required | Result |
|---|---:|---:|---|
| Ranked models | 5 | 5 | PASS |
| rank_001 mean pLDDT | 98.0674 | >= 70 | PASS |
| rank_001 pTM | 0.950 | valid | PASS |
| qTM vs experimental 8HFW | 0.9925 | >= 0.75 | PASS |
| tTM vs experimental 8HFW | 0.9925 | >= 0.75 | PASS |
| RMSD | 0.554 Å | reported | PASS |
| Catalytic mapping | S111->100, D241->230, H274->263; identities match | complete | PASS |
| Weight checksums | 5 | non-empty | PASS |

The mapped prediction positions are sequence indices. The reference labels retain the experimental PDB residue numbers 111, 241, and 274.

## Real Q9FVW3 Prediction

Input lineage: Q9FVW3 selected from a completed, annotated AlphaFold/Swiss-Prot mining run. The prediction manifest records both the source run and accession.

| Check | Observed | Required | Result |
|---|---:|---:|---|
| Ranked models | 5 | 5 | PASS |
| rank_001 mean pLDDT | 83.5407 | >= 70 | PASS |
| rank_001 pTM | 0.740 | valid | PASS |
| qTM vs AlphaFold DB v6 | 0.7815 | >= 0.75 | PASS |
| tTM vs AlphaFold DB v6 | 0.7815 | >= 0.75 | PASS |
| aligned coverage | 0.8247 / 0.8247 | reported | PASS |
| RMSD | 2.386 Å | reported | PASS |
| Candidate lineage | source run and Q9FVW3 recorded | traceable | PASS |

The lower confidence relative to 8HFW is concentrated by the full-length plant protein context, including an N-terminal targeting/low-complexity region. Acceptance is based on the stated global TM-score threshold, not on hiding that region.

## Predicted-Query Foldseek Benchmark

The rank_001 8HFW prediction was used as the real Foldseek query against the fixed seven-control database.

| Rank | Control | Role | Foldseek alignment score |
|---:|---|---|---:|
| 1 | 8HFW | self positive | 0.9975 |
| 2 | Q9FVW3 | related positive | 0.8450 |
| 3 | 6TJ2 | related positive | 0.8393 |
| 4 | Q6RYA0 | related positive | 0.8167 |
| 5 | 1XKL | related positive | 0.8091 |
| 6 | P00520 | remote negative | 0.2493 |
| 7 | P60174 | length-matched negative | 0.2443 |

Acceptance results:

- Self score `0.9975`, above the predicted-query requirement of `0.75`.
- Lowest positive `0.8091`, above highest negative `0.2493`.
- P00520 rejected by the default structural filter.
- Report identifies `query_mode=predicted` and links back to the prediction manifest.

## Complete Swiss-Prot Search

The same predicted rank_001 queried the complete pinned AlphaFold/Swiss-Prot v6 2025_03 Foldseek database.

| Metric | Observed |
|---|---:|
| Database structures | 590,183 |
| Raw Foldseek hits | 999 |
| Default-filter candidates | 12 |
| Valid annotations and sequences | 12 |
| Run warnings | 0 |
| Silent errors | 0 |
| Prediction/run lineage | bidirectional |

The earlier v0.2 experimental-query production result remains preserved; the 12-candidate count above is specifically for the v0.3 predicted query.

## UI Acceptance

The AF Prediction workbench was exercised with Playwright 1.61.0 against system Edge.

| Viewport | Horizontal overflow | Keyboard focus | pLDDT/PAE | py3Dmol pixels | Result |
|---|---|---|---|---|---|
| 1440x900 | none | visible | two nonblank charts | nonblank | PASS |
| 390x844 | none | visible | two nonblank charts | nonblank | PASS |

The interaction covered AF tab navigation, public-query FASTA entry, catalytic mapping preview, completed prediction selection, five-model ranking, pLDDT, PAE, and the confidence-colored 3D structure. Public screenshots contain only the public 8HFW query fixture and no production candidate sequence or machine path.

## Failure-Path Evidence And Fixes

1. ColabFold unified memory exhausted GPU allocation while loading five models. The run failed explicitly and did not fall back to CPU. Disabling unified memory reduced observed allocation to about 5.9 GB and preserved all scientific parameters.
2. A post-inference shell-variable quoting error prevented parameter-checksum collection. The completed model outputs were not promoted to success. Checksum collection now reads the verified parameter provenance JSON directly, with a shell-free hash fallback.
3. The public MSA service returned `too many open files`. The workflow detected ColabFold's semantic failure despite exit code zero, retained the upstream error, and failed after three attempts. A later immutable run succeeded and populated the sequence-hash cache.
4. Streamlit's default accent and deployment chrome did not match the specified workbench design. The fixed theme now uses structure green and hides the deployment menu; screenshots were regenerated and revalidated.

## Reproducibility, Privacy, And Packaging

- Prediction manifests include normalized input hash, MSA cache state, software versions, runtime patch provenance, parameter source and checksums, model records, confidence artifacts, catalytic mapping, events, warnings, errors, and linked run IDs.
- A cache miss requires explicit per-job consent before the sequence is sent to the public ColabFold MSA service.
- One public MSA request and one GPU task can run at a time.
- AlphaFold parameters remain under their bundled CC BY-NC-SA 4.0 license and are excluded from source and public archives.
- The public archive excludes databases, weights, caches, Obsidian data, real candidate results and sequences, logs, ZIP history, and machine-specific paths.

## Conclusion

All local v0.3 acceptance conditions pass. The result is a reproducible monomer-structure prediction and Bur01-like structure-mining workflow. It does not establish catalytic activity, S/H ratio, tunnel behavior, water-phase esterification, or a formal WASE model without the explicitly missing mechanistic and experimental evidence.
