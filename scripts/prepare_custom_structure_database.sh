#!/usr/bin/env bash
set -e

if [ "$#" -lt 1 ]; then
  echo "Usage: bash scripts/prepare_custom_structure_database.sh data/ncbi/predicted_structures [databases/ncbi_predicted_structures]"
  exit 1
fi

STRUCTURE_DIR="$1"
TARGET_DB="${2:-databases/ncbi_predicted_structures}"

mkdir -p databases "$(dirname "$TARGET_DB")"

if [ ! -d "$STRUCTURE_DIR" ]; then
  echo "Structure directory not found: $STRUCTURE_DIR"
  echo "Put AlphaFold/ColabFold PDB or mmCIF files in this directory first."
  exit 1
fi

foldseek createdb "$STRUCTURE_DIR" "$TARGET_DB"

echo "Created Foldseek database: $TARGET_DB"
