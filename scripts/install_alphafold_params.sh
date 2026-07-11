#!/usr/bin/env bash
set -euo pipefail

COLABFOLD_ROOT="${COLABFOLD_ROOT:-/opt/colabfold}"
ENV_PREFIX="${COLABFOLD_ENV:-${COLABFOLD_ROOT}/env}"
PARAMS_DIR="${COLABFOLD_ROOT}/data/params"
DOWNLOAD_DIR="${COLABFOLD_ROOT}/downloads"
ARCHIVE="${DOWNLOAD_DIR}/alphafold_params_2021-07-14.tar"
SOURCE_URL="https://storage.googleapis.com/alphafold/alphafold_params_2021-07-14.tar"
RUN_USER="${SUDO_USER:-aqua}"

if [[ "${EUID}" -ne 0 ]]; then
  echo "Run this parameter installer as root inside Ubuntu-Foldseek." >&2
  exit 2
fi
if [[ "$(realpath -m "${PARAMS_DIR}")" != "/opt/colabfold/data/params" ]]; then
  echo "Refusing to modify unexpected parameter path: ${PARAMS_DIR}" >&2
  exit 3
fi

install -d -m 0755 "${PARAMS_DIR}" "${DOWNLOAD_DIR}"
if ! command -v aria2c >/dev/null 2>&1; then
  apt-get update
  DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends aria2
fi

aria2c \
  --allow-overwrite=true \
  --auto-file-renaming=false \
  --continue=true \
  --file-allocation=none \
  --max-connection-per-server=16 \
  --min-split-size=5M \
  --split=16 \
  --dir="${DOWNLOAD_DIR}" \
  --out="$(basename "${ARCHIVE}")" \
  "${SOURCE_URL}"

archive_sha256="$(sha256sum "${ARCHIVE}" | cut -d' ' -f1)"
archive_size="$(stat -c '%s' "${ARCHIVE}")"

"${ENV_PREFIX}/bin/python" - "${ARCHIVE}" "${PARAMS_DIR}" <<'PY'
from __future__ import annotations

import shutil
import sys
import tarfile
from pathlib import Path

archive = Path(sys.argv[1]).resolve()
destination = Path(sys.argv[2]).resolve()
expected = Path("/opt/colabfold/data/params")
if destination != expected:
    raise SystemExit(f"Refusing unexpected destination: {destination}")
for child in destination.iterdir():
    if child.is_dir() and not child.is_symlink():
        shutil.rmtree(child)
    else:
        child.unlink()
with tarfile.open(archive, "r:") as bundle:
    bundle.extractall(destination, filter="data")
PY

ptm_count="$(find "${PARAMS_DIR}" -maxdepth 1 -type f -name 'params_model_*_ptm.npz' | wc -l)"
if [[ "${ptm_count}" -ne 5 ]]; then
  echo "Expected five AlphaFold2-PTM parameter files, found ${ptm_count}." >&2
  exit 4
fi

find "${PARAMS_DIR}" -maxdepth 1 -type f -name 'params_model_*_ptm.npz' -print0 \
  | sort -z | xargs -0 sha256sum > "${COLABFOLD_ROOT}/model-parameter-sha256.tsv"
touch "${PARAMS_DIR}/download_finished.txt"

"${ENV_PREFIX}/bin/python" - "${SOURCE_URL}" "${archive_sha256}" "${archive_size}" <<'PY'
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

checksums = {}
for line in Path("/opt/colabfold/model-parameter-sha256.tsv").read_text(encoding="utf-8").splitlines():
    checksum, path = line.split(maxsplit=1)
    checksums[Path(path).name] = checksum
payload = {
    "source_url": sys.argv[1],
    "archive_sha256": sys.argv[2],
    "archive_size_bytes": int(sys.argv[3]),
    "installed_at": datetime.now(timezone.utc).isoformat(),
    "model_type": "alphafold2_ptm",
    "parameter_files": checksums,
}
path = Path("/opt/colabfold/ALPHAFOLD_PARAMETER_PROVENANCE.json")
temporary = path.with_suffix(".json.tmp")
temporary.write_text(json.dumps(payload, indent=2), encoding="utf-8")
temporary.replace(path)
PY

chown -R "${RUN_USER}:${RUN_USER}" "${COLABFOLD_ROOT}"
echo "AlphaFold2-PTM parameters installed and verified."
cat "${COLABFOLD_ROOT}/ALPHAFOLD_PARAMETER_PROVENANCE.json"
