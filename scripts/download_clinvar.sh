#!/bin/zsh
set -euo pipefail
OUT_VCF="${1:-clinvar.vcf.gz}"
BASE_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"

if [[ ! -f "$OUT_VCF" ]]; then
  echo "Downloading ClinVar VCF to $OUT_VCF"
  wget -O "$OUT_VCF" "$BASE_URL"
else
  echo "ClinVar VCF already exists at $OUT_VCF"
fi

if [[ ! -f "${OUT_VCF}.tbi" ]]; then
  echo "Downloading ClinVar index ${OUT_VCF}.tbi"
  wget -O "${OUT_VCF}.tbi" "${BASE_URL}.tbi"
else
  echo "ClinVar index already exists at ${OUT_VCF}.tbi"
fi
