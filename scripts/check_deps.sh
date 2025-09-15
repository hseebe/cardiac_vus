#!/bin/zsh
set -euo pipefail

need() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Missing dependency: $1"
    return 1
  fi
}

BREW_CMD=${BREW_CMD:-brew}
MISSING=()
for bin in wget bcftools tabix bgzip unzip; do
  if ! command -v "$bin" >/dev/null 2>&1; then
    MISSING+=($bin)
  fi
done

if (( ${#MISSING[@]} > 0 )); then
  echo "Missing: ${MISSING[@]}"
  if command -v "$BREW_CMD" >/dev/null 2>&1; then
    echo "Attempting to install via Homebrew..."
    case " ${MISSING[@]} " in
      *" bcftools "*) $BREW_CMD install bcftools ;; esac
    case " ${MISSING[@]} " in
      *" wget "*) $BREW_CMD install wget ;; esac
    case " ${MISSING[@]} " in
      *" unzip "*) $BREW_CMD install unzip ;; esac
    # tabix and bgzip come with htslib
    case " ${MISSING[@]} " in
      *" tabix "*|*" bgzip "*) $BREW_CMD install htslib ;; esac
  else
    echo "Homebrew not found; please install missing deps manually."
    exit 1
  fi
fi

echo "All required CLI deps are available."
