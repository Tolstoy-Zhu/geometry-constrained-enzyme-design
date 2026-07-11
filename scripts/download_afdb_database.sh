#!/usr/bin/env bash
set -e
mkdir -p databases tmp

# AlphaFold DB search spaces supported by current Foldseek builds.
# Swiss-Prot is the reviewed subset and is a practical first target.
foldseek databases Alphafold/Swiss-Prot databases/afdb_swissprot tmp

# UniProt50 is clustered and broader than Swiss-Prot. It requires more storage.
# Uncomment when you are ready to search beyond the reviewed Swiss-Prot subset.
# foldseek databases Alphafold/UniProt50 databases/afdb_uniprot50 tmp

# The full UniProt or proteome-scale AlphaFold DB is very large.
# Check local storage and RAM before enabling either line below.
# foldseek databases Alphafold/UniProt databases/afdb_uniprot tmp
# foldseek databases Alphafold/Proteome databases/afdb_proteome tmp
