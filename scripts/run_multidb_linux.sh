#!/usr/bin/env bash
set -euo pipefail

CONFIG="config.yaml"
PYTHON_BIN="${PYTHON_BIN:-python3}"
VENV_DIR="${VENV_DIR:-.venv}"
PREPARE=0
RUN_PIPELINE=0
VALIDATE=0
WITH_NCBI_NETWORK=0
ONLY_ARGS=()
REQUIRE_SOURCES=()

usage() {
  cat <<'EOF'
Usage: bash scripts/run_multidb_linux.sh [options]

Options:
  --config PATH              Config file, default: config.yaml
  --prepare                  Download/build enabled Foldseek databases
  --only NAME [NAME ...]     Prepare only selected enabled database names
  --run                      Run the full pipeline after optional preparation
  --validate                 Run local validation and real-output checks
  --require-source NAME      Require source_db in final annotated output; can repeat
  --with-ncbi-network        Include live NCBI smoke test in validation
  -h, --help                 Show this help

Examples:
  bash scripts/run_multidb_linux.sh --validate
  bash scripts/run_multidb_linux.sh --prepare --only alphafold_swissprot
  bash scripts/run_multidb_linux.sh --prepare --run --validate --require-source alphafold_swissprot
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG="$2"
      shift 2
      ;;
    --prepare)
      PREPARE=1
      shift
      ;;
    --only)
      shift
      while [[ $# -gt 0 && "$1" != --* ]]; do
        ONLY_ARGS+=("$1")
        shift
      done
      ;;
    --run)
      RUN_PIPELINE=1
      shift
      ;;
    --validate)
      VALIDATE=1
      shift
      ;;
    --require-source)
      REQUIRE_SOURCES+=("$2")
      shift 2
      ;;
    --with-ncbi-network)
      WITH_NCBI_NETWORK=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  echo "Python not found: $PYTHON_BIN" >&2
  exit 1
fi

if [[ ! -d "$VENV_DIR" ]]; then
  "$PYTHON_BIN" -m venv "$VENV_DIR"
fi

# shellcheck source=/dev/null
source "$VENV_DIR/bin/activate"
python -m pip install --upgrade pip
python -m pip install -r requirements.txt

echo "[INFO] Python: $(python --version)"

python prepare_databases.py --config "$CONFIG" --dry-run

if [[ "$PREPARE" -eq 1 ]]; then
  if ! command -v foldseek >/dev/null 2>&1; then
    echo "Foldseek is required for --prepare. Install it first, for example:" >&2
    echo "  conda install -c conda-forge -c bioconda foldseek -y" >&2
    exit 1
  fi
  if [[ "${#ONLY_ARGS[@]}" -gt 0 ]]; then
    python prepare_databases.py --config "$CONFIG" --only "${ONLY_ARGS[@]}"
  else
    python prepare_databases.py --config "$CONFIG"
  fi
fi

if [[ "$RUN_PIPELINE" -eq 1 ]]; then
  if ! command -v foldseek >/dev/null 2>&1; then
    echo "Foldseek is required for --run. Install it first." >&2
    exit 1
  fi
  python run_pipeline.py --config "$CONFIG"
fi

if [[ "$VALIDATE" -eq 1 ]]; then
  VALIDATE_ARGS=()
  if [[ "$WITH_NCBI_NETWORK" -eq 1 ]]; then
    VALIDATE_ARGS+=(--with-ncbi-network)
  fi
  python scripts/validate_multidb_integration.py "${VALIDATE_ARGS[@]}"

  if [[ "$RUN_PIPELINE" -eq 1 || "${#REQUIRE_SOURCES[@]}" -gt 0 ]]; then
    REAL_ARGS=(--annotated results/Bur01_vs_database.annotated.xlsx)
    for source in "${REQUIRE_SOURCES[@]}"; do
      REAL_ARGS+=(--require-source "$source")
    done
    python scripts/validate_real_search_outputs.py "${REAL_ARGS[@]}"
  fi
fi

echo "[DONE] AquaEsterNet-StructMine Linux multi-database workflow finished."
