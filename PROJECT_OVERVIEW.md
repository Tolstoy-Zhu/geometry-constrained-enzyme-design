# AquaEsterNet-StructMine v0.3 Overview

AquaEsterNet-StructMine is the structure-search and evidence-integration module for Bur01-like enzyme discovery.

## Naming

- `AquaEsterNet`: broader water-phase ester mechanism and discovery framework.
- `StructMine`: structure mining, annotation, triage, and visualization module.
- `bur01_structmine`: Python package name retained for compatibility.
- `WASE-screen`: transparent v0.1 candidate screen, distinct from future formal WASE.

## Query

The fixed query is PDB 8HFW chain A, residues 12-294, with Ser111/Asp241/His274 as the catalytic reference.

## Main Routes

| Route | Command | Role |
|---|---|---|
| Persistent web workbench | `install_ui_autostart.bat`, then `start_ui.bat` | Fixed local port, logon startup, automatic recovery, and six bilingual views |
| Web status / stop | `ui_status.bat` / `stop_ui.bat` | Inspect or stop the local background service |
| Online expert demo | `install_online_ui.bat` | Password-protected HTTPS access to the same local compute service |
| Online status / stop | `online_ui_status.bat` / `stop_online_ui.bat` | Inspect or immediately close public access |
| Doctor | `python start_here.py --action doctor` | Environment and database diagnosis |
| Formal benchmark | `python start_here.py --action benchmark` | Fixed positive/negative Foldseek acceptance |
| TM-align preview | `python -m bur01_structmine.cli benchmark --engine tmtools` | Small-set preview only |
| Production | `python start_here.py --action production-run` | Complete AlphaFold/Swiss-Prot Foldseek search |
| AF2 prediction | `python -m bur01_structmine.cli predict ...` | Public MSA plus local GPU monomer prediction |
| Predict and mine | `python -m bur01_structmine.cli predict-and-mine ...` | Use AF2 rank_001 as a dynamic Foldseek query |
| Tests | `python start_here.py --action validate` | Isolated automated tests |

## Data Flow

Prediction route:

```text
FASTA/candidate -> validation -> public MSA cache -> local AF2-PTM -> rank_001 -> catalytic mapping -> Foldseek
```

Every prediction has an independent manifest under `results/predictions/`; prediction and mining manifests reference each other.

```text
8HFW:A -> Foldseek/TM-align -> structural filters -> catalytic mapping
        -> cached annotation -> WASE-screen -> run manifest and exports
```

## Result Contract

Each run under `results/runs/<run_id>/` includes a validated configuration snapshot, input hash, software and database versions, progress events, raw and filtered hits, annotations, WASE-screen fields, warnings, and failure or cancellation state.

## Verified Reference

- Standard 8HFW prediction: five models, mean pLDDT `98.07`, pTM `0.950`, experimental-structure minimum TM-score `0.9925`.
- Traceable Q9FVW3 prediction: five models, mean pLDDT `83.54`, pTM `0.740`, AlphaFold DB v6 minimum TM-score `0.7815`.
- Predicted-query Foldseek benchmark: PASS; lowest positive `0.8091`, highest negative `0.2493`, P00520 rejected.
- Full AlphaFold/Swiss-Prot v6 2025_03 search: 590,183 structures, 999 raw hits, 12 filtered candidates, 12 valid annotations, zero warnings.
- Automated evidence: 43 tests, 31 doctor checks, local/online recovery injection, and both desktop/mobile UI acceptance profiles passed.

## Main Documents

- `README_CN.md` / `README.md`: operating guides.
- `MODEL_GUIDE_CN.md` / `MODEL_GUIDE_EN.md`: model purpose, principle, innovation, and boundaries.
- `VALIDATION_REPORT.md`: current test and benchmark evidence.
- `COMPLETION_AUDIT.md`: requirement-by-requirement audit.
- `PROJECT_MANIFEST.md`: file inventory.
- `CHANGELOG.md`: version history.
