#!/bin/zsh
set -euo pipefail
OUT_GCT="${1:-data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz}"
ATTR_URLS=(
  "https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
  "https://gtexportal.org/home/datasets/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
)
mkdir -p "$(dirname "$OUT_GCT")"
if [[ ! -s "$OUT_GCT" ]]; then
  # Remove zero-byte or partial files
  [[ -f "$OUT_GCT" ]] && rm -f "$OUT_GCT"
  echo "Attempting to download GTEx gene median TPM (v8) to $OUT_GCT"
  urls=(
    "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
    "https://gtexportal.org/static/datasets/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
    "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.3.6_gene_median_tpm.gct.gz"
    "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9.gene_median_tpm.gct.gz"
    "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/gene_median_tpm.gct.gz"
    "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz" # per-sample fallback
  )
  success=0
  for u in ${urls[@]}; do
    echo "Trying $u"
    if wget -O "$OUT_GCT" "$u"; then
      success=1
      break
    fi
  done
  if [[ "$success" -ne 1 ]]; then
    echo "Failed to download GTEx GCT. Please update URL or download manually to $OUT_GCT"
    exit 1
      "https://gtexportal.org/static/datasets/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
  fi
else
  echo "GTEx GCT already exists at $OUT_GCT"
fi

ATTR_PATH="$(dirname "$OUT_GCT")/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
if [[ ! -s "$ATTR_PATH" ]]; then
  [[ -f "$ATTR_PATH" ]] && rm -f "$ATTR_PATH"
  echo "Downloading GTEx Sample Attributes to $ATTR_PATH"
  success=0
  for u in ${ATTR_URLS[@]}; do
    echo "Trying $u"
    if wget -O "$ATTR_PATH" "$u"; then
      success=1
      break
    fi
  done
  if [[ "$success" -ne 1 ]]; then
    echo "Failed to download GTEx sample attributes. You can download manually to $ATTR_PATH"
  fi
else
  echo "GTEx Sample Attributes already exists at $ATTR_PATH"
fi
