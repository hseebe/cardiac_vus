#!/bin/zsh
set -euo pipefail

CACHE_DIR="${VEP_CACHE:-$HOME/.vep}"
mkdir -p "$CACHE_DIR"

run_local_install() {
  local INSTALL=""
  if command -v vep_install >/dev/null 2>&1; then
    INSTALL=vep_install
  elif command -v INSTALL.pl >/dev/null 2>&1; then
    INSTALL=INSTALL.pl
  elif command -v conda >/dev/null 2>&1; then
    INSTALL="$(conda run which vep_install 2>/dev/null || true)"
    [[ -z "$INSTALL" ]] && INSTALL="$(conda run which INSTALL.pl 2>/dev/null || true)"
  fi
  if [[ -n "$INSTALL" ]]; then
    echo "Installing VEP cache locally using: $INSTALL"
    $INSTALL -a cf -s homo_sapiens -y GRCh38 -c "$CACHE_DIR" --CONVERT --force --CACHE_VERSION 111 || \
    $INSTALL -a cf -s homo_sapiens -y GRCh38 -c "$CACHE_DIR" --force
    return $?
  fi
  return 1
}

run_docker_install() {
  local TAG="ensemblorg/ensembl-vep:release_115.1"
  echo "Installing VEP cache via Docker image $TAG"
  if ! docker image inspect "$TAG" >/dev/null 2>&1; then
    docker pull "$TAG"
  fi
  # Use container path /opt/vep/.vep for cache, mount host cache there
  docker run --rm -v "$CACHE_DIR:/opt/vep/.vep" "$TAG" \
    vep_install -a cf -s homo_sapiens -y GRCh38 -c /opt/vep/.vep --CONVERT --force || \
  docker run --rm -v "$CACHE_DIR:/opt/vep/.vep" "$TAG" \
    vep_install -a cf -s homo_sapiens -y GRCh38 -c /opt/vep/.vep --force
}

if ! run_local_install; then
  echo "Local vep_install not available; trying Docker-based cache install..."
  run_docker_install
fi

echo "VEP cache installed in $CACHE_DIR"
