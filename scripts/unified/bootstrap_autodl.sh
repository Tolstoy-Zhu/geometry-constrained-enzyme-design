#!/usr/bin/env bash
# Install isolated runtimes on the data disk. Does not change the system driver.
set -euo pipefail
export PIP_INDEX_URL="${PIP_INDEX_URL:-https://pypi.org/simple}"
AQUA_DEPLOY_ROOT="${1:?Provide the deployment directory on the data disk}"
AQUA_BASE_PYTHON="${AQUA_BASE_PYTHON:-/root/miniconda3/bin/python}"
mkdir -p "$AQUA_DEPLOY_ROOT/envs" "$AQUA_DEPLOY_ROOT/logs" "$AQUA_DEPLOY_ROOT/models"
"$AQUA_BASE_PYTHON" -m venv "$AQUA_DEPLOY_ROOT/envs/aqua"
AQUA_PY="$AQUA_DEPLOY_ROOT/envs/aqua/bin/python"
"$AQUA_PY" -m pip install --upgrade pip
if [ -n "${AQUA_WHEELHOUSE:-}" ]; then
    # Operator-selected, previously checksum-verified wheels; resolve small dependencies below.
    "$AQUA_PY" -m pip install --no-deps "$AQUA_WHEELHOUSE"/*.whl
fi
"$AQUA_PY" -m pip install torch==2.9.1 --index-url "${AQUA_TORCH_INDEX:-https://download.pytorch.org/whl/cu128}" --report "$AQUA_DEPLOY_ROOT/logs/torch-install.json"
"$AQUA_PY" -m pip install numpy==1.26.4 pandas==2.2.3 scipy==1.15.3 scikit-learn==1.6.1 rdkit==2024.9.6 transformers==4.46.3 peft==0.13.2 pydantic==2.14.4 biopython==1.88 PyYAML==6.0.3 requests==2.32.5 tmtools==0.3.0 openpyxl==3.1.5 pytest==9.1.1 tqdm==4.67.1
"$AQUA_PY" -m pip install openmm==8.2.0 OpenMM-CUDA-12==8.2.0 nvidia-cuda-nvcc-cu12==12.8.93 MDAnalysis==2.9.0 propka==3.5.1
"$AQUA_PY" -m pip install tensorflow-cpu==2.15.1
"$AQUA_PY" -m pip freeze > "$AQUA_DEPLOY_ROOT/logs/runtime-freeze.txt"
"$AQUA_PY" -m openmm.testInstallation > "$AQUA_DEPLOY_ROOT/logs/openmm-installation.txt" 2>&1 || true
"$AQUA_PY" -c 'import torch; x=torch.ones((1024,1024),device="cuda"); print(torch.__version__,torch.version.cuda,torch.cuda.get_device_name(0),float((x@x).mean())); assert torch.cuda.is_available()'
printf 'BOOTSTRAP_COMPLETE\n'
