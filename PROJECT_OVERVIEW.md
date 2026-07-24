# AquaEsterNet-StructMine v0.4 Overview

AquaEsterNet-StructMine is the structure-search and evidence-integration module for Bur01-like enzyme discovery.

## Naming

- `AquaEsterNet`: broader water-phase ester mechanism and discovery framework.
- `StructMine`: structure mining, NR sequence expansion, annotation, triage, and visualization module.
- `bur01_structmine`: Python package name retained for compatibility.
- `WASE-screen`: transparent v0.1 candidate screen, distinct from future formal WASE.
- `NR expansion`: independent sequence search, clustering, and catalytic-identity mapping; it is not structural evidence by itself.
- `NR novelty`: nearest-neighbor audit for predicted or designed sequences; natural-sequence self-hits are not novelty evidence.

## Query

The fixed query is PDB 8HFW chain A, residues 12-294, with Ser111/Asp241/His274 as the catalytic reference.

## Main Routes

| Route | Command | Role |
|---|---|---|
| Persistent web workbench | `install_ui_autostart.bat`, then `start_ui.bat` | Fixed local port, logon startup, automatic recovery, and seven bilingual views |
| Web status / stop | `ui_status.bat` / `stop_ui.bat` | Inspect or stop the local background service |
| Online expert demo | `install_online_ui.bat` | Password-protected HTTPS access to the same local compute service |
| Online status / stop | `online_ui_status.bat` / `stop_online_ui.bat` | Inspect or immediately close public access |
| Doctor | `python start_here.py --action doctor` | Environment and database diagnosis |
| Formal benchmark | `python start_here.py --action benchmark` | Fixed positive/negative Foldseek acceptance |
| TM-align preview | `python -m bur01_structmine.cli benchmark --engine tmtools` | Small-set preview only |
| Production | `python start_here.py --action production-run` | Complete AlphaFold/Swiss-Prot Foldseek search |
| AF2 prediction | `python -m bur01_structmine.cli predict ...` | Public MSA plus local GPU monomer prediction |
| Predict and mine | `python -m bur01_structmine.cli predict-and-mine ...` | Use AF2 rank_001 as a dynamic Foldseek query |
| Local NR | `python -m bur01_structmine.cli nr-expand --source-run <id>` | Search the real query-derived NCBI NR hit snapshot |
| Full NR export | `python -m bur01_structmine.cli nr-export-job --source-run <id>` | Build a credential-free Bash/Slurm/PBS server package |
| Full NR import | `python -m bur01_structmine.cli nr-import-job --bundle <zip>` | Verify and import a complete-NR server result |
| NR novelty | `python -m bur01_structmine.cli nr-novelty --fasta <path> ...` | Audit the nearest NR sequence neighbor |
| Tests | `python start_here.py --action validate` | Isolated automated tests |

## Data Flow

Prediction route:

```text
FASTA/candidate -> validation -> public MSA cache -> local AF2-PTM -> rank_001 -> catalytic mapping -> Foldseek
```

Every prediction has an independent manifest under `results/predictions/`; prediction and mining manifests reference each other.

NR route:

```text
Structure run -> Foldseek target sequences -> MMseqs2/NR -> accession and SHA-256 deduplication
        -> 50% clustering -> direct 8HFW catalytic mapping -> AF2-ready queue
        -> AlphaFold2 -> Foldseek/TM-align -> WASE-screen
```

```text
8HFW:A -> Foldseek/TM-align -> structural filters -> catalytic mapping
        -> cached annotation -> WASE-screen -> run manifest and exports
```

## Result Contract

Each run under `results/runs/<run_id>/` includes a validated configuration snapshot, input hash, software and database versions, progress events, raw and filtered hits, annotations, WASE-screen fields, warnings, and failure or cancellation state.

Each NR run under `results/nr_runs/<nr_run_id>/` includes database scope, fixed MMseqs2 parameters, seed hashes, raw and deduplicated hits, clusters, direct 8HFW alignments, sequence-screen exports, AF2 lineage, warnings, and checksums. Sequence-only candidates keep WASE explicitly unevaluated.

## Verified Reference

- Standard 8HFW prediction: five models, mean pLDDT `98.07`, pTM `0.950`, experimental-structure minimum TM-score `0.9925`.
- Traceable Q9FVW3 prediction: five models, mean pLDDT `83.54`, pTM `0.740`, AlphaFold DB v6 minimum TM-score `0.7815`.
- Predicted-query Foldseek benchmark: PASS; lowest positive `0.8091`, highest negative `0.2493`, P00520 rejected.
- Full AlphaFold/Swiss-Prot v6 2025_03 search: 590,183 structures, 999 raw hits, 12 filtered candidates, 12 valid annotations, zero warnings.
- v0.4 automated evidence: 65 isolated tests pass; MMseqs2 `18.8cc5c` and NCBI BLAST+ `2.17.0` pass diagnostics.
- Complete-NR server acceptance remains pending, so the development branch is not eligible for merge, tag, or Release.

## Main Documents

- `README_CN.md` / `README.md`: operating guides.
- `MODEL_GUIDE_CN.md` / `MODEL_GUIDE_EN.md`: model purpose, principle, innovation, and boundaries.
- `VALIDATION_REPORT.md`: current test and benchmark evidence.
- `COMPLETION_AUDIT.md`: requirement-by-requirement audit.
- `PROJECT_MANIFEST.md`: file inventory.
- `CHANGELOG.md`: version history.
- `docs/NR_WORKFLOW_CN.md` / `docs/NR_WORKFLOW_EN.md`: NR evidence boundary and local/server operation.
