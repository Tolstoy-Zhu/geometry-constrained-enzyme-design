#!/usr/bin/env bash
set -e
# Recommended installation with conda/mamba:
#   conda create -n bur01_structmine python=3.10 -y
#   conda activate bur01_structmine
#   conda install -c conda-forge -c bioconda foldseek -y
# Or download Foldseek binary from the official GitHub release page.

foldseek version || echo "Foldseek is not installed yet. Please install it with conda or the official binary."
