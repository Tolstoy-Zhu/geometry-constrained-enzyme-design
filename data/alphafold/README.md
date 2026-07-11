# AlphaFold Data

`structures/` is reserved for small accession-selected AlphaFold PDB files used by TM-align previews or custom Foldseek databases.

Download selected accessions:

```powershell
python scripts/download_alphafold_structures.py --accessions Q6RYA0 Q9FVW3
```

Production search does not enumerate this folder. It uses the complete Foldseek `Alphafold/Swiss-Prot` database configured in `configs/alphafold_swissprot.yaml` and stored under `databases/afdb_swissprot`.

The validated local database is `v6 2025_03` with 590,183 structures. The pinned official archive is 1,557,648,318 bytes with SHA-256 `2ae623b90e74ab8a42b6840f61f2ed3edfa2d4e006e04d6afa0eec5284c73d35`. Installation evidence is stored under `databases/environment/`.

This 32 GB workstation intentionally does not download UniProt50. Its memory requirement is outside the supported local profile.

v0.3 local prediction is separate from AlphaFold DB search. LocalColabFold and AlphaFold2-PTM parameters live inside the WSL distribution under `/opt/colabfold`; generated MSAs and models live under ignored cache/result directories. Model parameters, A3M files, prediction logs, and real candidate outputs are never included in source or release archives.

The prediction runtime lock is `configs/prediction_runtime.lock.yaml`; the exact installed Python/CUDA package set is `configs/colabfold-pip-lock.txt`. Every prediction manifest records parameter-file checksums after the first verified download. The AlphaFold parameter archive and five PTM files are SHA-256 pinned and remain under their bundled CC BY-NC-SA 4.0 license.

AlphaFold hits must be interpreted with pLDDT, PAE, UniProt reviewed status, alignment coverage, and experimental evidence. A predicted structure is not proof of function.
