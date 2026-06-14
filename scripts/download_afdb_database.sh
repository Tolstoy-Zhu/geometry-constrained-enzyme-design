#!/usr/bin/env bash
set -e
mkdir -p databases tmp
# This may require large disk space. Start from smaller databases if available.
# Run `foldseek databases` to check supported database names in your version.
foldseek databases Alphafold/Proteome databases/afdb_proteome tmp
