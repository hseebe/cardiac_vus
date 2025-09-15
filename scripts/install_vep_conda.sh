#!/bin/zsh
set -euo pipefail

# Install Mambaforge locally and a local VEP environment, then create a wrapper at scripts/vep
ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
CONDA_DIR="$ROOT_DIR/.mambaforge"
ENV_DIR="$CONDA_DIR/envs/vep"
WRAPPER="$ROOT_DIR/scripts/vep"

# Always (re)create wrapper after ensuring env

if [[ ! -x "$CONDA_DIR/bin/conda" ]]; then
  echo "Installing Mambaforge into $CONDA_DIR ..."
  KERNEL=$(uname -s)
  ARCH=$(uname -m)
  TMP_SH="/tmp/Mambaforge.sh"
  if [[ "$KERNEL" == "Darwin" ]]; then
    if [[ "$ARCH" == "arm64" ]]; then
      URL="https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-arm64.sh"
    else
      URL="https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-x86_64.sh"
    fi
  else
    if [[ "$ARCH" == "aarch64" || "$ARCH" == "arm64" ]]; then
      URL="https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-aarch64.sh"
    else
      URL="https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh"
    fi
  fi
  echo "Downloading $URL"
  if curl -fsSL "$URL" -o "$TMP_SH"; then
    if [[ $(stat -f%z "$TMP_SH" 2>/dev/null || stat -c%s "$TMP_SH" 2>/dev/null) -lt 1000000 ]]; then
      echo "Mambaforge installer seems too small; falling back to Miniconda."
      USE_MINI=1
    else
      bash "$TMP_SH" -b -p "$CONDA_DIR" || USE_MINI=1
    fi
  else
    echo "Failed to download Mambaforge; falling back to Miniconda."
    USE_MINI=1
  fi
  if [[ "${USE_MINI:-0}" -eq 1 ]]; then
    echo "Installing Miniconda into $CONDA_DIR ..."
    if [[ "$KERNEL" == "Darwin" ]]; then
      if [[ "$ARCH" == "arm64" ]]; then
        MINI_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh"
      else
        MINI_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
      fi
    else
      if [[ "$ARCH" == "aarch64" || "$ARCH" == "arm64" ]]; then
        MINI_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh"
      else
        MINI_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
      fi
    fi
    MINI_SH="/tmp/Miniconda.sh"
    curl -fsSL "$MINI_URL" -o "$MINI_SH"
    bash "$MINI_SH" -b -p "$CONDA_DIR"
  fi
fi

MAMBA="$CONDA_DIR/bin/mamba"
CONDA="$CONDA_DIR/bin/conda"
if [[ ! -x "$MAMBA" ]]; then
  MAMBA="$CONDA_DIR/bin/conda"
fi

echo "Creating VEP conda env at $ENV_DIR ..."
"$CONDA" create -y --override-channels -p "$ENV_DIR" -c conda-forge -c bioconda ensembl-vep || true

VEP_BIN="$ENV_DIR/bin/vep"

cat > "$WRAPPER" <<'EOF'
#!/bin/zsh
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
ENV_DIR="$ROOT_DIR/.mambaforge/envs/vep"
"$ROOT_DIR/.mambaforge/bin/conda" run -p "$ENV_DIR" vep "$@"
EOF
chmod +x "$WRAPPER"

if [[ -x "$VEP_BIN" ]]; then
  echo "VEP wrapper installed: $WRAPPER"
else
  echo "VEP environment creation may have failed; you can still try running the wrapper which will use conda run."
fi
