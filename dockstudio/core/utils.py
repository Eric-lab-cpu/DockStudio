"""Small shared helpers (no heavy third-party imports)."""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import unicodedata
import zlib
from datetime import datetime
from typing import Optional, Tuple


def crc32_seed(receptor: str, ligand: str) -> int:
    """Deterministic random seed = crc32('<RECEPTOR>|<LIGAND>') & 0x7fffffff."""
    return zlib.crc32(f"{receptor}|{ligand}".encode("utf-8")) & 0x7FFFFFFF


def short_path(path: str) -> str:
    """Return the Windows short (8.3) name for ``path``.

    RDKit / meeko / Vina use narrow (ANSI) file APIs on Windows and can fail
    on paths that contain non-ASCII characters (e.g. Chinese folder names).
    The 8.3 short path is pure ASCII and points to the same folder, so all
    engine file writes become path-safe while the user still sees their own
    (possibly Chinese) output directory. No-op on non-Windows platforms.
    """
    if not path or os.name != "nt":
        return path
    try:
        import ctypes
        buf = ctypes.create_unicode_buffer(1024)
        r = ctypes.windll.kernel32.GetShortPathNameW(str(path), buf, len(buf))
        if r and r < len(buf):
            short = buf.value
            if short:
                return short
    except Exception:
        pass
    return path


def now_str() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def safe_name(text: str) -> str:
    """Turn arbitrary input name into a filesystem-safe ASCII-ish token."""
    text = unicodedata.normalize("NFKD", str(text))
    text = text.encode("ascii", "ignore").decode("ascii")
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("._")
    return text or "mol"


def ascii_only(text: str) -> str:
    return unicodedata.normalize("NFKD", str(text)).encode("ascii", "ignore").decode("ascii")


def run_cmd(cmd, timeout_s: Optional[int] = 3600, cwd=None) -> Tuple[int, str, str]:
    """Run a subprocess; returns (returncode, stdout, stderr)."""
    env = dict(os.environ)
    env.setdefault("PYTHONUNBUFFERED", "1")
    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=timeout_s,
        cwd=cwd,
        env=env,
    )
    return proc.returncode, proc.stdout or "", proc.stderr or ""


def which(*names: str) -> Optional[str]:
    """Find the first executable among names on PATH."""
    for n in names:
        p = shutil.which(n)
        if p:
            return p
    return None


def write_json(obj, path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, ensure_ascii=False, indent=2, sort_keys=False)
    os.replace(tmp, path)


def read_json(path: str):
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def md5_file(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def distance(a, b) -> float:
    return math.dist(a, b)


def pdb_atom_name_for_element(atom_name: str, element: str) -> str:
    """PDB atom-name field is 4 chars; align by PDB convention."""
    if len(element) == 1:
        return f" {atom_name:>3s}"[:4]
    return f"{atom_name:<4s}"[:4]


def python_executable() -> str:
    return sys.executable or "python"


def truncate(text: str, n: int = 500) -> str:
    return text if len(text) <= n else text[: n - 3] + "..."
