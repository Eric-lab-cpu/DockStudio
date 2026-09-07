#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# DockStudio v2.0 release orchestration (Linux/CI host)
#   Windows NSIS installer + Linux conda-pack portable + MD5 list + Release notes
#
# This runs on a Linux CI host that already has:
#   - a conda env with packaging/environment.yml (name = ENV_NAME below)
#   - conda-forge 'nsis' installed (so makensis is available)
#   - gh (GitHub CLI) authenticated when PUBLISH=1
#
# Usage:
#   bash packaging/ci/release_all.sh [ENV_NAME] [PUBLISH=0|1]
# ---------------------------------------------------------------------------
set -euo pipefail

ENV_NAME="${1:-dockstudio}"
PUBLISH="${2:-0}"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
VERSION="$(python -c 'from dockstudio._version import __version__;print(__version__)')"
mkdir -p dist

echo "==> DockStudio v${VERSION} release build"

# 1) Linux conda-pack portable
bash packaging/cross/build_condapack_unix.sh "$ENV_NAME" dist

# 2) Windows NSIS online installer (makensis on Linux compiles the .exe)
#    prepare payload then compile the .nsi
python packaging/online_installer/prepare_payload.py
makensis packaging/online_installer/DockStudio_online.nsi
mv dist/DockStudio_Setup_*.exe dist/DockStudio_Setup_${VERSION}.exe 2>/dev/null || true

# 3) checksums + release notes file
( cd dist && md5sum DockStudio_Setup_${VERSION}.exe DockStudioPortable.tar.gz \
    > DockStudio_${VERSION}_MD5.txt && cat DockStudio_${VERSION}_MD5.txt )

cat > dist/Release_v${VERSION}.md <<EOF
# DockStudio v${VERSION}
See docs/定版记录_v${VERSION}.md for the changelog, honest known limitations and test record.
Assets: Windows NSIS installer (exe), Linux conda-pack portable (tar.gz), MD5 list.
EOF

echo "==> release assets ready in dist/"
ls -la dist | grep -E "DockStudio_(Setup|Portable|${VERSION})" || true

if [ "${PUBLISH}" = "1" ]; then
  echo "==> publish (gh release) - manual step template:"
  echo "  gh release create v${VERSION} dist/DockStudio_Setup_${VERSION}.exe dist/DockStudioPortable.tar.gz dist/DockStudio_${VERSION}_MD5.txt -F dist/Release_v${VERSION}.md"
fi
