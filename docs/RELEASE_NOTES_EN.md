# AquaEsterNet-StructMine v0.3.0 Release Notes

## Scope

v0.3.0 adds local AlphaFold2-PTM monomer prediction to the traceable v0.2 structure-mining workflow. A single FASTA or an existing candidate can be predicted alone or predicted and mined, with rank_001 passed to the same Foldseek service layer as a dynamic query.

## Main Changes

- Isolated LocalColabFold 1.6.1, AlphaFold2-PTM, JAX/CUDA, and OpenMM inside `Ubuntu-Foldseek`.
- Strict single-protein input: 20-1000 aa, no more than 5% `X`, with aggregated validation errors.
- Standard preset: five models, three recycles, one seed, and rank_001 Amber relaxation; templates remain off by default.
- Separate public MSA and local GPU stages, sequence-hash caching, serialized requests, per-job transfer consent, and up to three exponential retries.
- One job per GPU, cancellation, stale-lock recovery, explicit memory failures, and no CPU inference fallback.
- Global mapping of 8HFW `Ser111/Asp241/His274` with manual correction; identity failures never receive fabricated catalytic evidence.
- Bidirectional prediction/mining lineage with runtime versions, parameter provenance and SHA-256, ranking, pLDDT, pTM, PAE, events, warnings, and failures.
- A new AF Prediction web view with candidate re-prediction, rankings, pLDDT, PAE, and confidence-colored 3D structure.
- Python 3.11 non-GPU CI, deterministic fake-ColabFold end-to-end coverage, sanitized packaging, and repeatable Playwright UI acceptance.

## Local Acceptance

- Automated tests: 34 passed.
- Doctor: 31 passed, 0 warnings, 0 failures.
- 8HFW: five models; rank_001 mean pLDDT 98.07, pTM 0.950, minimum bidirectional TM-score 0.9925 against experiment.
- Q9FVW3: five models; rank_001 mean pLDDT 83.54, pTM 0.740, minimum bidirectional TM-score 0.7815 against AlphaFold DB v6.
- Predicted-8HFW benchmark: lowest positive 0.8091, above highest negative 0.2493; P00520 rejected.
- Full AlphaFold/Swiss-Prot: 590,183 structures, 999 raw hits, 12 filtered candidates, 12 valid annotations, no warnings.
- UI: screenshot, overflow, keyboard-focus, chart, and nonblank WebGL checks passed at 1440x900 and 390x844.

## Privacy And Licensing

A cache-miss sequence is sent to the public ColabFold MSA service only after explicit per-job confirmation. AlphaFold parameters remain under their bundled CC BY-NC-SA 4.0 license and are not included in the repository or Release assets. The public archive also excludes databases, caches, Obsidian data, real candidate results, model weights, logs, and machine-specific paths.

## Scientific Boundaries

AlphaFold prediction and structural similarity do not establish activity, S/H ratio, or a water-phase esterification mechanism. WASE-screen v0.1 remains transparent triage, not the formal WASE model. Interpret candidates with pLDDT, PAE, structural coverage, annotation quality, and experiments.

## Upgrade Notes

Existing v0.2 results, databases, Obsidian vault, and delivery archives remain unchanged. v0.3 predictions write to `results/predictions/`; linked mining runs continue under `results/runs/`.
