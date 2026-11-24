#!/usr/bin/env bash
# Python 3.12 and a CUDA-13-compatible driver must already be available.
set -euo pipefail
AQUA_DEPLOY_ROOT="${1:?Provide the deployment directory on the data disk}"
AQUA_STRUCTURE_BASE_PYTHON="${AQUA_STRUCTURE_BASE_PYTHON:-python3.12}"
AQUA_SCRIPT_ROOT="$(cd "$(dirname "$0")" && pwd)"
"$AQUA_STRUCTURE_BASE_PYTHON" -c 'import sys; assert sys.version_info[:2] == (3,12), "Use Python 3.12"'
mkdir -p "$AQUA_DEPLOY_ROOT/bin" "$AQUA_DEPLOY_ROOT/envs" "$AQUA_DEPLOY_ROOT/models/alphafold" "$AQUA_DEPLOY_ROOT/logs"
"$AQUA_STRUCTURE_BASE_PYTHON" -m venv "$AQUA_DEPLOY_ROOT/envs/structure"
AQUA_STRUCTURE_PY="$AQUA_DEPLOY_ROOT/envs/structure/bin/python"
"$AQUA_STRUCTURE_PY" -m pip install -r "$AQUA_SCRIPT_ROOT/../../configs/unified/structure-linux-py312.lock.txt" --report "$AQUA_DEPLOY_ROOT/logs/structure-install.json"
"$AQUA_STRUCTURE_PY" -m pip check
install -m 755 "$AQUA_SCRIPT_ROOT/structure-python" "$AQUA_DEPLOY_ROOT/bin/structure-python"
install -m 755 "$AQUA_SCRIPT_ROOT/colabfold-batch" "$AQUA_DEPLOY_ROOT/bin/colabfold-batch"
"$AQUA_DEPLOY_ROOT/bin/structure-python" -c 'import jax,jax.numpy as j; x=j.ones((256,256)); assert jax.default_backend()=="gpu"; print(jax.devices(), float((x@x).mean()))'
"$AQUA_DEPLOY_ROOT/bin/structure-python" -c 'import sys; from pathlib import Path; from colabfold.download import download_alphafold_params; download_alphafold_params("alphafold2_ptm",Path(sys.argv[1]))' "$AQUA_DEPLOY_ROOT/models/alphafold"
"$AQUA_STRUCTURE_PY" "$AQUA_SCRIPT_ROOT/verify_af_parameters.py" --data "$AQUA_DEPLOY_ROOT/models/alphafold"
printf 'STRUCTURE_RUNTIME_READY\n'
