#!/usr/bin/env bash
set -euo pipefail

VERSION="2.17.0"
ARCHIVE_NAME="ncbi-blast-${VERSION}+-x64-linux.tar.gz"
URL="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${VERSION}/${ARCHIVE_NAME}"
EXPECTED_MD5="bdec166721de3b55f90a3badc83538e8"
INSTALL_ROOT="/opt/ncbi-blast"
VERSION_DIR="${INSTALL_ROOT}/${VERSION}"
ARCHIVE="${TMPDIR:-/tmp}/${ARCHIVE_NAME}"

DOWNLOADED_ARCHIVE=0
if [[ -n "${AQUA_BLAST_ARCHIVE:-}" ]]; then
  ARCHIVE="$AQUA_BLAST_ARCHIVE"
else
  curl --fail --location --retry 3 --output "$ARCHIVE" "$URL"
  DOWNLOADED_ARCHIVE=1
fi
printf '%s  %s\n' "$EXPECTED_MD5" "$ARCHIVE" | md5sum --check -
SHA256="$(sha256sum "$ARCHIVE" | awk '{print $1}')"

sudo mkdir -p "$INSTALL_ROOT"
if [[ ! -x "$VERSION_DIR/bin/blastp" ]]; then
  STAGING="$(sudo mktemp -d "${INSTALL_ROOT}/.${VERSION}.installing.XXXXXX")"
  sudo tar -xzf "$ARCHIVE" --strip-components=1 -C "$STAGING"
  sudo mv "$STAGING" "$VERSION_DIR"
fi
sudo chmod -R a+rX "$VERSION_DIR"
if [[ -e "$INSTALL_ROOT/bin" && ! -L "$INSTALL_ROOT/bin" ]]; then
  sudo mv "$INSTALL_ROOT/bin" "$INSTALL_ROOT/bin.backup.$(date -u +%Y%m%dT%H%M%SZ)"
fi
sudo ln -sfn "$VERSION_DIR/bin" "$INSTALL_ROOT/bin"
sudo tee "$VERSION_DIR/INSTALL_PROVENANCE.json" >/dev/null <<EOF
{
  "version": "$VERSION",
  "source_url": "$URL",
  "official_md5": "$EXPECTED_MD5",
  "archive_sha256": "$SHA256"
}
EOF
if [[ "$DOWNLOADED_ARCHIVE" == "1" ]]; then
  rm -f "$ARCHIVE"
fi
"$INSTALL_ROOT/bin/blastp" -version
printf 'archive_sha256=%s\n' "$SHA256"
