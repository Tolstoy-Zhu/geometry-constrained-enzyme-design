#!/usr/bin/env bash
set -euo pipefail

COLABFOLD_ROOT="${COLABFOLD_ROOT:-/opt/colabfold}"
ENV_PREFIX="${COLABFOLD_ENV:-${COLABFOLD_ROOT}/env}"
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-${COLABFOLD_ROOT}/mamba}"
MICROMAMBA_BIN="${COLABFOLD_ROOT}/bin/micromamba"
MICROMAMBA_URL="https://micro.mamba.pm/api/micromamba/linux-64/latest"
PYPI_INDEX_URL="${PYPI_INDEX_URL:-https://pypi.tuna.tsinghua.edu.cn/simple}"
PYPI_FALLBACK_URL="${PYPI_FALLBACK_URL:-https://pypi.org/simple}"
RUN_USER="${SUDO_USER:-aqua}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ "${EUID}" -ne 0 ]]; then
  echo "Run this installer as root inside Ubuntu-Foldseek." >&2
  exit 2
fi

install -d -m 0755 "${COLABFOLD_ROOT}/bin" "${COLABFOLD_ROOT}/downloads" \
  "${COLABFOLD_ROOT}/data" "${MAMBA_ROOT_PREFIX}"

apt-get update
DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
  bzip2 ca-certificates curl git libgomp1 wget

if [[ ! -x "${MICROMAMBA_BIN}" ]]; then
  archive="${COLABFOLD_ROOT}/downloads/micromamba-linux-64.tar.bz2"
  curl --fail --location --retry 3 --output "${archive}" "${MICROMAMBA_URL}"
  tar -xjf "${archive}" -C "${COLABFOLD_ROOT}/downloads" bin/micromamba
  install -m 0755 "${COLABFOLD_ROOT}/downloads/bin/micromamba" "${MICROMAMBA_BIN}"
  sha256sum "${archive}" > "${archive}.sha256"
fi

export MAMBA_ROOT_PREFIX
if [[ ! -x "${ENV_PREFIX}/bin/python" ]]; then
  "${MICROMAMBA_BIN}" create --yes --prefix "${ENV_PREFIX}" \
    --channel conda-forge --channel bioconda --strict-channel-priority \
    python=3.13.14 pip kalign2=2.04 hhsuite=3.3.0 mmseqs2=18.8cc5c
fi

"${ENV_PREFIX}/bin/python" -m pip install \
  --index-url "${PYPI_INDEX_URL}" \
  --upgrade "pip==26.1.2"
PIP_LOCK="${SCRIPT_DIR}/../configs/colabfold-pip-lock.txt"
if [[ ! -f "${PIP_LOCK}" ]]; then
  echo "Pinned ColabFold package lock not found: ${PIP_LOCK}" >&2
  exit 5
fi
if ! "${ENV_PREFIX}/bin/python" -m pip install \
  --index-url "${PYPI_INDEX_URL}" --requirement "${PIP_LOCK}"; then
  "${ENV_PREFIX}/bin/python" -m pip install \
    --index-url "${PYPI_FALLBACK_URL}" --requirement "${PIP_LOCK}"
fi

"${ENV_PREFIX}/bin/python" "${SCRIPT_DIR}/patch_colabfold_result_timeout.py"

cat > "${COLABFOLD_ROOT}/bin/colabfold_batch" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export PATH="${ENV_PREFIX}/bin:\${PATH}"
export COLABFOLD_DATA="${COLABFOLD_ROOT}/data"
export XLA_PYTHON_CLIENT_PREALLOCATE="false"
exec "${ENV_PREFIX}/bin/colabfold_batch" --data "${COLABFOLD_ROOT}/data" "\$@"
EOF
chmod 0755 "${COLABFOLD_ROOT}/bin/colabfold_batch"

{
  echo "installed_at_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "installer=https://github.com/sokrypton/ColabFold/tree/v1.6.1"
  echo "micromamba_url=${MICROMAMBA_URL}"
  echo "pypi_index_url=${PYPI_INDEX_URL}"
  echo "pypi_fallback_url=${PYPI_FALLBACK_URL}"
  echo "micromamba_sha256=$(cut -d' ' -f1 "${COLABFOLD_ROOT}/downloads/micromamba-linux-64.tar.bz2.sha256")"
  echo "python=$(${ENV_PREFIX}/bin/python --version 2>&1)"
  echo "colabfold=$(${ENV_PREFIX}/bin/python -c 'import importlib.metadata as m; print(m.version("colabfold"))')"
  echo "jax=$(${ENV_PREFIX}/bin/python -c 'import jax; print(jax.__version__)')"
  echo "openmm=$(${ENV_PREFIX}/bin/python -c 'import openmm; print(openmm.__version__)')"
  echo "kalign=$(${ENV_PREFIX}/bin/kalign -h 2>&1 | grep -m1 'Kalign version')"
  echo "hhblits=$(${ENV_PREFIX}/bin/hhblits -h 2>&1 | head -n 1)"
  echo "mmseqs=$(${ENV_PREFIX}/bin/mmseqs version 2>&1 | head -n 1)"
  echo "colabfold_patch_provenance=${COLABFOLD_ROOT}/COLABFOLD_PATCH_PROVENANCE.json"
} > "${COLABFOLD_ROOT}/INSTALL_PROVENANCE.txt"

# `pip freeze` can preserve a conda build-host file URL for pip itself. The
# package-list form records portable name/version pins instead.
"${ENV_PREFIX}/bin/python" -m pip list --format=freeze > "${COLABFOLD_ROOT}/pip-freeze.txt"
"${MICROMAMBA_BIN}" list --prefix "${ENV_PREFIX}" --json > "${COLABFOLD_ROOT}/conda-lock.json"

chown -R "${RUN_USER}:${RUN_USER}" "${COLABFOLD_ROOT}"

"${ENV_PREFIX}/bin/python" - <<'PY'
import jax
import openmm

devices = jax.devices()
print("JAX devices:", devices)
if not any(device.platform == "gpu" for device in devices):
    raise SystemExit("No JAX GPU device detected; CPU fallback is not accepted.")
print("OpenMM platforms:", [openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())])
PY

echo "LocalColabFold installation completed at ${COLABFOLD_ROOT}."
