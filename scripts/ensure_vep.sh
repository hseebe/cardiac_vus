#!/bin/zsh
set -euo pipefail
if command -v vep >/dev/null 2>&1; then
  echo "VEP is installed: $(vep --help | head -n 1)"
  exit 0
fi
if command -v brew >/dev/null 2>&1; then
  echo "Trying Homebrew for VEP (formula may be unavailable)..."
  brew install ensembl-vep || true
  if command -v vep >/dev/null 2>&1; then
    echo "VEP installed via Homebrew."
    exit 0
  fi
fi

if command -v conda >/dev/null 2>&1; then
  echo "Installing VEP via conda (bioconda)..."
  conda install -y -c bioconda ensembl-vep || true
  if command -v vep >/dev/null 2>&1; then
    echo "VEP installed via conda."
    exit 0
  fi
fi

if command -v mamba >/dev/null 2>&1; then
  echo "Installing VEP via mamba (bioconda)..."
  mamba install -y -c bioconda ensembl-vep || true
  if command -v vep >/dev/null 2>&1; then
    echo "VEP installed via mamba."
    exit 0
  fi
fi

cat <<EOF
VEP not installed automatically.
Options:
  1) Install via conda: conda install -c bioconda ensembl-vep
  2) Install manually from Ensembl docs: https://www.ensembl.org/info/docs/tools/vep/index.html
After installing, run: make vep-cache
EOF
exit 1
