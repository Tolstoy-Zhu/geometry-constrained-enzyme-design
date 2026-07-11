# Fixed Bur01-like Benchmark

This directory contains the v0.1 fixed acceptance set.

## Controls

- Self positive: PDB 8HFW chain A.
- Related positives: PDB 1XKL, PDB 6TJ2, AlphaFold Q6RYA0, AlphaFold Q9FVW3.
- Length-matched negative: AlphaFold P60174.
- Remote negative: AlphaFold P00520.

`benchmark_manifest.yaml` defines roles and acceptance criteria. `download_manifest.json` records the downloaded file, source URL, model version, chains, residue counts, and SHA-256.

The `structures/` files are used for both the TM-align preview and the Foldseek formal benchmark. Multi-chain PDB entries are compared one chain at a time.

## Acceptance

- 8HFW self TM-score > 0.98.
- 8HFW self RMSD < 0.1 Å.
- All positives rank above both negatives.
- P00520 fails the default structural filter.
- Formal acceptance uses Foldseek.

Prepare or refresh sources:

```powershell
python scripts/prepare_benchmark.py
```

Run formal acceptance:

```powershell
python scripts/run_benchmark.py --engine foldseek
```

Run preview only:

```powershell
python scripts/run_benchmark.py --engine tmtools
```
