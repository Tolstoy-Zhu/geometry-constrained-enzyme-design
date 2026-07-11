# NCBI Candidate Data

NCBI Protein is sequence-based. Structural search requires one of two routes:

1. Predict or obtain PDB/mmCIF structures and place them in `predicted_structures/`.
2. Convert a candidate FASTA set into a Foldseek-searchable ProstT5 representation.

Download a sequence set:

```powershell
python scripts/download_ncbi_protein_fasta.py --query "esterase[Protein Name] AND bacteria[Organism]" --max-records 200 --out data/ncbi/ncbi_candidates.fasta
```

Prepare the ProstT5 route:

```powershell
python prepare_databases.py --config configs/ncbi_prostt5.yaml --only ncbi_prostt5
```

Preview a small predicted-structure folder with TM-align:

```powershell
python -m bur01_structmine.cli run --config configs/ncbi_predicted_structures_python.yaml
```

TM-align preview rows contain real qTM/tTM/RMSD and leave Foldseek-only significance fields empty. NCBI annotation failures are reported in `annotation_errors`.
