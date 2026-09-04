#!/usr/bin/env python3
"""Stage the online-installer payload for DockStudio (Eric Studio).

Layout produced (default: <repo>/dist/payload):

    DockStudio.bat            - launcher; creates env on first run then starts GUI
    icon.ico                  - installer / shortcut icon
    scripts/
        setup_env.bat         - creates conda-forge env (PyMOL/RDKit/etc) via micromamba
        micromamba.exe        - Windows micromamba (self-contained, ~11 MB)
    vina/
        vina.exe              - official AutoDock Vina 1.2.7 Windows binary
    app/                      - full application source tree
        dockstudio/  assets/  examples/demo/  README.md  USER_MANUAL.md  LICENSE  pyproject.toml
"""

from __future__ import annotations

import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RUNTIME = os.path.join(HERE, "runtime")
DEFAULT_OUT = os.path.join(ROOT, "dist", "payload")

COPY_APP_DIRS = ["dockstudio", "assets", "examples"]
COPY_APP_FILES = ["README.md", "用户手册.md", "LICENSE", "THIRD_PARTY_NOTICES.md", "pyproject.toml"]


def _copy_tree(src: str, dst: str, ignore_pycache: bool = True):
    if ignore_pycache:
        def ign(d, names):
            return {n for n in names if n == "__pycache__" or n.endswith(".pyc")}
        shutil.copytree(src, dst, ignore=ign)
    else:
        shutil.copytree(src, dst)


def main(out: str = DEFAULT_OUT) -> str:
    if os.path.exists(out):
        shutil.rmtree(out)
    os.makedirs(out)
    os.makedirs(os.path.join(out, "scripts"))
    os.makedirs(os.path.join(out, "vina"))
    app = os.path.join(out, "app")
    os.makedirs(app)

    # launcher + icon
    shutil.copy(os.path.join(HERE, "DockStudio.bat"), os.path.join(out, "DockStudio.bat"))
    shutil.copy(os.path.join(HERE, "DockStudio_debug.bat"), os.path.join(out, "DockStudio_debug.bat"))
    shutil.copy(os.path.join(HERE, "Install_PyMOL.bat"), os.path.join(out, "Install_PyMOL.bat"))
    shutil.copy(os.path.join(RUNTIME, "icon.ico"), os.path.join(out, "icon.ico"))
    # scripts
    shutil.copy(os.path.join(HERE, "setup_env.bat"), os.path.join(out, "scripts", "setup_env.bat"))
    shutil.copy(os.path.join(HERE, "fetch_micromamba.ps1"), os.path.join(out, "scripts", "fetch_micromamba.ps1"))
    shutil.copy(os.path.join(RUNTIME, "micromamba.exe"), os.path.join(out, "scripts", "micromamba.exe"))
    # VC++ runtime DLLs (msvcp140.dll etc.) - needed by micromamba.exe on machines
    # without the Visual C++ Redistributable; Windows loads them from this folder.
    vcrt_src = os.path.join(RUNTIME, "vcrt")
    for dll in os.listdir(vcrt_src):
        if dll.lower().endswith(".dll"):
            shutil.copy(os.path.join(vcrt_src, dll), os.path.join(out, "scripts", dll))
    # vina
    shutil.copy(os.path.join(RUNTIME, "vina_1.2.7_win.exe"), os.path.join(out, "vina", "vina.exe"))
    # VC++ runtime DLLs also next to vina.exe (in case system redistributable absent)
    for dll in os.listdir(vcrt_src):
        if dll.lower().endswith(".dll"):
            shutil.copy(os.path.join(vcrt_src, dll), os.path.join(out, "vina", dll))
    # app
    for d in COPY_APP_DIRS:
        src = os.path.join(ROOT, d)
        if os.path.isdir(src):
            _copy_tree(src, os.path.join(app, d))
    for f in COPY_APP_FILES:
        src = os.path.join(ROOT, f)
        if os.path.exists(src):
            shutil.copy(src, os.path.join(app, f))
    # keep examples/demo only (remove demo_out if it slipped in)
    demo_out = os.path.join(app, "examples", "demo_out")
    if os.path.isdir(demo_out):
        shutil.rmtree(demo_out)
    print("payload staged at:", out)
    return out


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else DEFAULT_OUT)
