#!/bin/zsh
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
VEP_SRC_DIR="$ROOT_DIR/.vep_src"
VEP_VER="111.0"
VEP_TGZ="https://github.com/Ensembl/ensembl-vep/archive/refs/tags/release/${VEP_VER}.tar.gz"
CACHE_DIR="${VEP_CACHE:-$HOME/.vep}"

mkdir -p "$VEP_SRC_DIR"
cd "$VEP_SRC_DIR"

if [[ ! -f "vep_${VEP_VER}.tar.gz" ]]; then
  echo "Downloading VEP source ${VEP_VER}..."
  curl -fsSL "$VEP_TGZ" -o "vep_${VEP_VER}.tar.gz"
fi

if [[ ! -d "ensembl-vep-release-${VEP_VER}" ]]; then
  echo "Extracting VEP..."
  tar -xzf "vep_${VEP_VER}.tar.gz"
fi

cd "ensembl-vep-release-${VEP_VER}"

# Ensure cpanminus
if ! command -v cpanm >/dev/null 2>&1; then
  echo "Bootstrapping cpanminus..."
  curl -fsSL https://cpanmin.us | perl - App::cpanminus
fi

echo "Installing VEP Perl dependencies (this may take a while)..."
perl INSTALL.pl -a a -q || true

echo "Installing VEP cache (GRCh38) into $CACHE_DIR ..."
perl INSTALL.pl -a cf -s homo_sapiens -y GRCh38 -c "$CACHE_DIR" -q || true

# Create wrapper to the local VEP script
WRAPPER="$ROOT_DIR/scripts/vep"
BIN_VEP="$VEP_SRC_DIR/ensembl-vep-release-${VEP_VER}/vep"
if [[ -f "$BIN_VEP" ]]; then
  cat > "$WRAPPER" <<EOF
#!/bin/zsh
  export PERL5LIB="$VEP_SRC_DIR/ensembl-vep-release-${VEP_VER}/modules:
  \$HOME/.vep:\$PERL5LIB"
  exec perl "$BIN_VEP" "$@"
EOF
  chmod +x "$WRAPPER"
  echo "VEP wrapper installed: $WRAPPER"
else
  echo "ERROR: VEP script not found at $BIN_VEP"
  exit 1
fi

echo "VEP source installation attempted. If modules are missing at runtime, rerun this script to let INSTALL.pl resolve them."
