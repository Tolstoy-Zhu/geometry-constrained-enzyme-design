#!/usr/bin/env bash
set -e
mkdir -p databases tmp
# This command requires Foldseek to be installed.
# The database name may vary slightly by Foldseek version.
# Run `foldseek databases` to check supported names.
foldseek databases PDB100 databases/pdb100 tmp
