# -*- mode: python ; coding: utf-8 -*-
# PyInstaller spec for DockStudio GUI.
# NOTE: PyMOL / Vina / meeko are large native packages; the RECOMMENDED
# distribution is conda-pack (see packaging/condapack + 打包说明.md).
# This spec is provided for advanced users who want a single .exe.
# Build:  pyinstaller packaging/pyinstaller/dockstudio.spec

import os
from PyInstaller.utils.hooks import collect_all, collect_submodules

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
ICON = os.path.join(ROOT, "assets", "icon.ico")
VERSION_FILE = os.path.join(os.path.dirname(__file__), "version_info.txt")

hidden = collect_submodules("vina") + collect_submodules("meeko") + \
    collect_submodules("plip") + collect_submodules("gemmi")

datas = [(os.path.join(ROOT, "dockstudio", "resources"), "dockstudio/resources")]
binaries = []
for pkg in ("vina", "meeko", "openbabel", "pymol"):
    try:
        d, b, _ = collect_all(pkg)
        datas += d
        binaries += b
    except Exception:
        pass

a = Analysis(
    [os.path.join(ROOT, "dockstudio", "gui", "app.py")],
    pathex=[ROOT],
    binaries=binaries,
    datas=datas,
    hiddenimports=hidden,
    hookspath=[],
    runtime_hooks=[],
    excludes=["PyQt5.QtWebEngineWidgets"],
    noarchive=False,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="DockStudio",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    icon=ICON,
    version=VERSION_FILE,
    disable_windowed_traceback=False,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    name="DockStudio",
)
