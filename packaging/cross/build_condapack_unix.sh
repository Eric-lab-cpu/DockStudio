#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# DockStudio (C) 2026 Eric Studio - cross-platform conda-pack green directory
# for Linux / macOS (v2.0). Windows NSIS installer is produced separately by
# packaging/online_installer/build_online_installer.bat / .sh
#
# What it does
#   Given an already-created conda env (see packaging/environment.yml), pack it
#   with conda-pack into a relocatable green directory that already contains the
#   whole engine toolchain (PyMOL/RDKit/Meeko/PLIP/Vina...), then add the
#   DockStudio source, an ASCII launcher and this repo's README.
#
# Usage
#   bash packaging/cross/build_condapack_unix.sh [conda-env-name] [out-dir]
#   default env : dockstudio
#   default out : dist
# ---------------------------------------------------------------------------
set -euo pipefail

ENV_NAME="${1:-dockstudio}"
OUT_DIR="${2:-dist}"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

cd "$ROOT"
mkdir -p "$OUT_DIR"

echo "==> conda-pack env '$ENV_NAME' -> $OUT_DIR/DockStudioPortable.tar.gz"
conda pack -n "$ENV_NAME" -o "$OUT_DIR/DockStudioPortable.tar.gz" -q

echo "==> inflate green directory"
rm -rf "$OUT_DIR/DockStudioPortable"
mkdir -p "$OUT_DIR/DockStudioPortable"
tar -xzf "$OUT_DIR/DockStudioPortable.tar.gz" -C "$OUT_DIR/DockStudioPortable"

echo "==> add DockStudio source + docs"
mkdir -p "$OUT_DIR/DockStudioPortable/app"
cp -r dockstudio "$OUT_DIR/DockStudioPortable/app/"
cp -r examples "$OUT_DIR/DockStudioPortable/app/" 2>/dev/null || true
cp -r docs "$OUT_DIR/DockStudioPortable/app/" 2>/dev/null || true
cp README.md "$OUT_DIR/DockStudioPortable/app/" 2>/dev/null || true
cp THIRD_PARTY_NOTICES.md "$OUT_DIR/DockStudioPortable/app/" 2>/dev/null || true

echo "==> write launcher DockStudio.sh"
cat > "$OUT_DIR/DockStudio.sh" <<'EOF'
#!/usr/bin/env bash
# DockStudio launcher (Linux/macOS). Uses the bundled python (relocatable).
HERE="$(cd "$(dirname "$0")" && pwd)"
export PYTHONPATH="$HERE/app${PYTHONPATH:+:$PYTHONPATH}"
exec "$HERE/bin/python" -m dockstudio "$@"
EOF
chmod +x "$OUT_DIR/DockStudio.sh"

echo "==> checksum"
( cd "$OUT_DIR" && md5sum DockStudioPortable.tar.gz > DockStudioPortable.tar.gz.md5 && cat DockStudioPortable.tar.gz.md5 )

echo "==> done."
echo "    Linux/macOS portable dir : $OUT_DIR/DockStudioPortable  (launcher: $OUT_DIR/DockStudio.sh)"
echo "    Compressed archive       : $OUT_DIR/DockStudioPortable.tar.gz"
