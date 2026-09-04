#!/usr/bin/env bash
# DockStudio PyInstaller build (Linux/macOS, advanced)
set -e
cd "$(dirname "$0")/../.."
python -m PyInstaller --noconfirm packaging/pyinstaller/dockstudio.spec
echo "Build finished: dist/DockStudio/"
