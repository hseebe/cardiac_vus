#!/bin/zsh
set -euo pipefail

IN_VCF="$1"
OUT_VCF="$2"
DBNSFP_FILE="${3:-}"

# Resolve absolute paths and repo root for Docker mount
ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
ABS_IN="$(cd "$(dirname "$IN_VCF")" && pwd)/$(basename "$IN_VCF")"
ABS_OUT_DIR="$(cd "$(dirname "$OUT_VCF")" && pwd)"
ABS_OUT="$ABS_OUT_DIR/$(basename "$OUT_VCF")"
ABS_DBNSFP=""
if [[ -n "${DBNSFP_FILE}" ]]; then
  ABS_DBNSFP="$(cd "$(dirname "$DBNSFP_FILE")" 2>/dev/null && pwd 2>/dev/null)/$(basename "$DBNSFP_FILE")"
fi

# Determine VEP path (avoid broken local wrapper if its runtime is missing)
VEP_BIN=""
if [[ -x "$(dirname "$0")/vep" ]]; then
  # Verify the wrapper's conda exists; otherwise ignore it
  if grep -q "/.mambaforge/bin/conda" "$(dirname "$0")/vep" && [[ ! -x "$(cd "$(dirname "$0")/.." && pwd)/.mambaforge/bin/conda" ]]; then
    VEP_BIN=""
  else
    VEP_BIN="$(dirname "$0")/vep"
  fi
fi
if [[ -z "$VEP_BIN" ]] && command -v vep >/dev/null 2>&1; then
  VEP_BIN="vep"
fi

# Cache dir (host). If using Docker, this will be mounted to /opt/vep/.vep
# Cache directory; use offline mode only if cache present
CACHE_DIR="${VEP_CACHE:-$ROOT_DIR/.vep_cache}"
mkdir -p "$CACHE_DIR"
USE_CACHE=0
if [[ -d "$CACHE_DIR/homo_sapiens" ]] || [[ -n "$(find "$CACHE_DIR" -maxdepth 2 -type d -name 'homo_sapiens' 2>/dev/null)" ]]; then
  USE_CACHE=1
fi

build_args_local() {
  local args=(
    --species homo_sapiens
    --assembly GRCh38
    --vcf
    --format vcf
    --force_overwrite
    --canonical
    --symbol
    --fork 4
    --input_file "$ABS_IN"
    --output_file "$ABS_OUT"
    --everything
  )
  if [[ $USE_CACHE -eq 1 ]]; then
    args+=( --offline --cache --dir_cache "$CACHE_DIR" )
  else
    args+=( --database --no_stats )
  fi
  echo "${args[@]}"
}
ARGS=( $(build_args_local) )

# dbNSFP plugin if available
if [[ -n "${ABS_DBNSFP}" && -f "${ABS_DBNSFP}" ]]; then
  ARGS+=( --plugin dbNSFP,"${ABS_DBNSFP}",Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,SIFT_pred,MutationTaster_pred,GERP++_RS,phyloP100way_vertebrate )
fi

run_local() {
  [[ -z "$VEP_BIN" ]] && return 1
  echo "Running VEP locally: $VEP_BIN"
  "$VEP_BIN" ${ARGS[@]}
}

run_docker() {
  local TAG="ensemblorg/ensembl-vep:release_115.1"
  echo "Falling back to Docker: $TAG"
  # Ensure image exists
  if ! docker image inspect "$TAG" >/dev/null 2>&1; then
    docker pull "$TAG"
  fi
  # Map host paths into container
  local WORK_MNT="$ROOT_DIR:/work"
  local CACHE_MNT="$CACHE_DIR:/opt/vep/.vep"
  # Convert host absolute paths to container paths under /work
  local CIN="$ABS_IN"
  local COUT="$ABS_OUT"
  local CDB="$ABS_DBNSFP"
  if [[ "$CIN" == "$ROOT_DIR"* ]]; then CIN="/work${CIN#$ROOT_DIR}"; fi
  if [[ "$COUT" == "$ROOT_DIR"* ]]; then COUT="/work${COUT#$ROOT_DIR}"; fi
  if [[ -n "$CDB" && "$CDB" == "$ROOT_DIR"* ]]; then CDB="/work${CDB#$ROOT_DIR}"; fi

  # Build container args; use container cache path and updated IO paths
  local CARGS=(
    --species homo_sapiens --assembly GRCh38
    --vcf --format vcf --force_overwrite --canonical --symbol --fork 4
    --input_file "$CIN" --output_file "$COUT" --everything
  )
  if [[ $USE_CACHE -eq 1 ]]; then
    CARGS+=( --offline --cache --dir_cache /opt/vep/.vep )
  else
    CARGS+=( --database --no_stats )
  fi
  if [[ -n "$CDB" && -f "$ABS_DBNSFP" ]]; then
    CARGS+=( --plugin dbNSFP,"$CDB",Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,SIFT_pred,MutationTaster_pred,GERP++_RS,phyloP100way_vertebrate )
  fi

  # Use current user to avoid root-owned outputs
  local UIDGID
  UIDGID="$(id -u):$(id -g)"
  docker run --rm -u "$UIDGID" -v "$WORK_MNT" -v "$CACHE_MNT" "$TAG" \
    vep ${CARGS[@]}
}

# Try local first; on failure, fall back to Docker
set +e
if [[ -n "$VEP_BIN" ]]; then
  run_local
  rc=$?
else
  rc=1
fi
set -e
if [[ $rc -ne 0 ]]; then
  echo "Local VEP unavailable/failed (exit $rc); attempting Docker..."
  run_docker
fi

echo "VEP annotation complete: $OUT_VCF"
