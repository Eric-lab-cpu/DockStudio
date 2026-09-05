#!/usr/bin/env bash
# ============================================================
# DockStudio (Eric Studio) - build Windows online installer
# Works on Linux/macOS/Windows with NSIS available.
# Usage: bash packaging/online_installer/build_online_installer.sh
# Output: dist/DockStudio_Setup_1.1.0.exe
# ============================================================
set -e
cd "$(dirname "$0")/../.."

# 1. stage payload (app + micromamba + vina)
python3 packaging/online_installer/prepare_payload.py

# 2. find makensis
MAKENSIS=$(command -v makensis || true)
if [ -z "$MAKENSIS" ]; then
  if [ -n "$NSISDIR" ] && [ -x "$NSISDIR/../bin/makensis" ]; then
    MAKENSIS="$NSISDIR/../bin/makensis"
  else
    echo "ERROR: makensis not found. Install NSIS (apt install nsis / brew install nsis / chocolatey) or set NSISDIR."
    exit 1
  fi
fi

# 3. compile
"$MAKENSIS" packaging/online_installer/DockStudio_online.nsi
echo "DONE: dist/DockStudio_Setup_1.1.0.exe"
