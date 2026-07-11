#!/usr/bin/env bash
set -euo pipefail

url="https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz"
archive="${1:-foldseek-linux-avx2.tar.gz}"
install_root="${FOLDSEEK_HOME:-/opt/foldseek}"

grep -qw avx2 /proc/cpuinfo || {
  echo "AVX2 is required for this build." >&2
  exit 1
}

if [[ ! -f "$archive" ]]; then
  curl -fL --retry 3 "$url" -o "$archive"
fi

echo "Foldseek archive SHA-256: $(sha256sum "$archive" | awk '{print $1}')"
mkdir -p "$install_root"
tar -xzf "$archive" -C "$install_root" --strip-components=1
chmod 0755 "$install_root/bin/foldseek"
"$install_root/bin/foldseek" version
