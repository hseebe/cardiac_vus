#!/bin/zsh
set -euo pipefail
OUT_DIR="${1:-data/conservation}"
mkdir -p "$OUT_DIR"

# UCSC URLs for 100-way vertebrate conservation (hg38)
PHYLOP_URL="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw"
PHASTCONS_URL="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.phastCons100way.bw"

for url in "$PHYLOP_URL" "$PHASTCONS_URL"; do
  fname=$(basename "$url")
  if [[ ! -f "$OUT_DIR/$fname" ]]; then
    echo "Downloading $fname"
    wget -O "$OUT_DIR/$fname" "$url"
  else
    echo "$fname already exists"
  fi
done
