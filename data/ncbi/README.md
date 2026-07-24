# NCBI NR Data Routes

## v0.4 Supported Route

NCBI NR is a sequence database. v0.4 keeps its evidence independent from Foldseek/TM-align structural evidence.

The local demonstration database is created by one explicitly confirmed Bur01 remote `blastp nr` request:

```powershell
$env:AQUA_NCBI_EMAIL = "name@example.org"
python -m bur01_structmine.cli nr-demo-build `
  --query-fasta data/query/Bur01m_seed.fasta `
  --max-hits 2500 `
  --confirm-ncbi
```

The ignored output is stored under `databases/nr_demo/`. Its manifest must say `nr_remote_snapshot`; it is never described as complete NR.

Run sequence expansion from a structure run:

```powershell
python -m bur01_structmine.cli nr-expand `
  --source-run <run_id> `
  --config configs/nr_local_demo.yaml
```

For complete NR, export a credential-free server/HPC job:

```powershell
python -m bur01_structmine.cli nr-export-job `
  --source-run <run_id> `
  --config configs/nr_full_server.yaml
```

Complete NR is built only by the explicit server command `mmseqs databases NR`. It is not downloaded during a normal search. Import the returned result with `nr-import-job`; the importer validates scope, version, paths, and SHA-256 before exposing candidates to the web workbench.

See `docs/NR_WORKFLOW_CN.md` and `docs/NR_WORKFLOW_EN.md`.

## Legacy Small-FASTA Route

`configs/ncbi_prostt5.yaml` remains compatible with a small user-defined FASTA:

```powershell
python scripts/download_ncbi_protein_fasta.py `
  --query "esterase[Protein Name] AND bacteria[Organism]" `
  --max-records 200 `
  --out data/ncbi/ncbi_candidates.fasta

python prepare_databases.py `
  --config configs/ncbi_prostt5.yaml `
  --only ncbi_prostt5
```

This legacy route is not complete NR and must not be used for formal v0.4 NR acceptance.

## Structural Alternatives

Predicted or experimental PDB/mmCIF files may still be placed under `predicted_structures/` and searched with the existing Foldseek or TM-align configuration. TM-align rows contain real qTM/tTM/RMSD and leave Foldseek-only significance fields empty.

Databases, real NR candidates, caches, run logs, server configurations, and returned ZIP files are ignored by Git and excluded from public releases.
