#!/usr/bin/env bash
set -e

if [ "$#" -lt 1 ]; then
  echo "Usage: bash scripts/prepare_ncbi_prostt5_database.sh data/ncbi/ncbi_candidates.fasta [databases/ncbi_prostt5]"
  exit 1
fi

INPUT_FASTA="$1"
TARGET_DB="${2:-databases/ncbi_prostt5}"
PROSTT5_DB="${PROSTT5_DB:-databases/prostt5_weights}"
TMP_DIR="${TMP_DIR:-tmp}"

mkdir -p databases "$TMP_DIR" "$(dirname "$TARGET_DB")"

if [ ! -f "$INPUT_FASTA" ]; then
  echo "NCBI FASTA not found: $INPUT_FASTA"
  echo "Export protein FASTA records from NCBI Protein or NCBI Datasets first."
  exit 1
fi

if [ ! -e "$PROSTT5_DB" ]; then
  foldseek databases ProstT5 "$PROSTT5_DB" "$TMP_DIR"
fi

foldseek createdb "$INPUT_FASTA" "$TARGET_DB" --prostt5-model "$PROSTT5_DB"

echo "Created Foldseek database: $TARGET_DB"
