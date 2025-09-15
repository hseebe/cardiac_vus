#!/bin/zsh
set -euo pipefail

# Usage: scripts/download_dbnsfp.sh [OUT_DIR]
# Env opts:
#   DBNSFP_FORCE_REBUILD=1   Force rebuild of merged txt.gz even if it exists

OUT_DIR="${1:-data/dbnsfp}"
mkdir -p "$OUT_DIR"
ZIP_PATH="$OUT_DIR/dbNSFP4.4a.zip"
URL="ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFP4.4a.zip"

if [[ ! -f "$ZIP_PATH" ]]; then
  echo "Downloading dbNSFP to $ZIP_PATH"
  if ! wget -O "$ZIP_PATH" "$URL"; then
    echo "Failed to download dbNSFP (FTP may require updated credentials)."
    echo "Please download manually from: $URL or official dbNSFP site, then place dbNSFP4.4a.zip into $OUT_DIR and re-run."
    exit 1
  fi
else
  echo "dbNSFP zip already exists at $ZIP_PATH"
fi

# Unzip (idempotent)
echo "Unzipping (idempotent) into $OUT_DIR"
unzip -qn "$ZIP_PATH" -d "$OUT_DIR" || true

# Discover per-chromosome files (either root or nested folder)
typeset -a CHR_FILES
for c in {1..22} X Y M; do
  if [[ -f "$OUT_DIR/dbNSFP4.4a_variant.chr${c}.gz" ]]; then
    CHR_FILES+=("$OUT_DIR/dbNSFP4.4a_variant.chr${c}.gz")
  elif [[ -f "$OUT_DIR/dbNSFP4.4a/dbNSFP4.4a_variant.chr${c}.gz" ]]; then
    CHR_FILES+=("$OUT_DIR/dbNSFP4.4a/dbNSFP4.4a_variant.chr${c}.gz")
  fi
done

if (( ${#CHR_FILES[@]} == 0 )); then
  echo "Could not locate per-chromosome dbNSFP files after unzip. Contents:" >&2
  ls -lh "$OUT_DIR" >&2 || true
  exit 1
fi

MERGED="$OUT_DIR/dbNSFP4.4a.txt.gz"

validate_merged() {
  local f="$1"
  [[ -f "$f" ]] || return 1
  # Must be larger than a trivial header and contain the expected first columns
  if [[ $(stat -f%z "$f") -lt 1000000 ]]; then
    return 1
  fi
  local hdr
  if ! hdr=$(zcat -f "$f" 2>/dev/null | head -n 1); then
    return 1
  fi
  # Expect columns like: chr  pos(ref)  ref  alt ... (header names may vary slightly, accept 'chr' and 'pos')
  echo "$hdr" | grep -qiE '(^|\t)chr(\t|$)' || return 1
  echo "$hdr" | grep -qiE '(^|\t)pos(\t|$)' || return 1
  return 0
}

rebuild_merged() {
  local merged="$1"
  echo "Rebuilding merged dbNSFP: $merged"
  rm -f "$merged" "$merged.tbi"
  local first_file="${CHR_FILES[1]}"
  local header
  header=$(gzip -dc "$first_file" | head -n 1)
  # Stream header once, then append all body lines skipping their headers
  (
    print -r -- "$header"
    for f in "${CHR_FILES[@]}"; do
      echo "Appending $(basename "$f")" >&2
      gzip -dc "$f" | tail -n +2
    done
  ) | bgzip -@ 4 -c > "$merged"
  echo "Indexing $merged with tabix"
  tabix -f -s 1 -b 2 -e 2 "$merged"
}

if [[ -n "${DBNSFP_FORCE_REBUILD:-}" ]]; then
  rebuild_merged "$MERGED"
elif validate_merged "$MERGED"; then
  echo "Existing merged dbNSFP looks valid: $MERGED"
else
  echo "Preparing merged dbNSFP (this may take time)."
  rebuild_merged "$MERGED"
fi

if ! validate_merged "$MERGED"; then
  echo "ERROR: Merged dbNSFP at $MERGED failed validation." >&2
  exit 1
fi

echo "dbNSFP download/prep complete: $MERGED"
