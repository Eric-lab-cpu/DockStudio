"""Environment & tool-version introspection (section 3 of the protocol)."""

from __future__ import annotations

import importlib
import os
import platform
import sys
from typing import Optional

from . import utils

#: full tool inventory (version reporting)
TOOLS = [
    "numpy", "scipy", "pandas", "PIL", "rdkit", "gemmi",
    "meeko", "plip", "vina", "pymol", "openbabel",
]

#: tools that DockStudio requires at runtime (GUI warns if missing)
REQUIRED = [
    "numpy", "pandas", "PIL", "rdkit", "gemmi",
    "meeko", "plip", "openbabel", "vina",
]

#: tools that are optional (their absence only degrades specific features)
OPTIONAL = ["scipy", "pymol", "prody"]

BINARIES = {
    "vina": ["vina"],
    "pymol": ["pymol"],
    "plip": ["plip"],
    "mk_prepare_ligand": ["mk_prepare_ligand.py", "mk_prepare_ligand"],
    "mk_prepare_receptor": ["mk_prepare_receptor.py", "mk_prepare_receptor"],
    "obabel": ["obabel"],
}


def version_of(module_name: str) -> Optional[str]:
    try:
        mod = importlib.import_module(module_name)
        v = getattr(mod, "__version__", None)
        if v is None and module_name == "PIL":
            v = getattr(mod, "PILLOW_VERSION", None)
        if module_name == "pymol":
            try:
                v = mod.cmd.get_version()[0]
            except Exception:
                v = None
        return str(v) if v else "?"
    except Exception:
        return None


def find_vina_binary() -> Optional[str]:
    """Locate a usable Vina engine: python bindings OR CLI executable."""
    for key in ("DOCKSTUDIO_VINA", "VINA_BIN"):
        p = os.environ.get(key)
        if p and os.path.exists(p):
            return p
    return utils.which("vina", "vina.exe")


def vina_status() -> Optional[str]:
    """Return a human readable Vina engine status: python version, CLI path,
    or None if completely unavailable."""
    py = version_of("vina")
    if py:
        return f"python {py}"
    cli = find_vina_binary()
    if cli:
        return f"CLI {os.path.basename(cli)}"
    return None


def python_info() -> dict:
    return {
        "python": sys.version.split()[0],
        "executable": sys.executable,
        "platform": platform.platform(),
        "system": platform.system(),
        "release": platform.release(),
    }


def gather_environment() -> dict:
    env = {"python": python_info(), "tools": {}}
    for m in TOOLS:
        env["tools"][m] = version_of(m)
    env["tools"]["vina"] = vina_status()
    env["binaries"] = {k: utils.which(*names) for k, names in BINARIES.items()}
    env["missing_required"] = [t for t in REQUIRED if not env["tools"].get(t)]
    env["missing_optional"] = [t for t in OPTIONAL if not env["tools"].get(t)]
    env["collected_at"] = utils.now_str()
    return env


def write_environment_json(path: str) -> dict:
    env = gather_environment()
    utils.write_json(env, path)
    return env
