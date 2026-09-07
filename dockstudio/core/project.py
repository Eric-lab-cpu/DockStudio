"""Project save/load and data-package export (v2.0, features #4-part / #18).

Two jobs live here:

1. **``.dsproj`` project files** - a thin JSON wrapper around the existing
   ``RunConfig`` (which is still saved as ``config.json`` inside every output
   directory for full backward compatibility). A project file stores the whole
   ``RunConfig`` plus format/version metadata, so a run can be re-opened and
   re-executed from the command line or from a future GUI "open project".

2. **One-click reproducible data-package export** - packs *all* real source
   files of an output directory (config/state, prepared structures, docking
   results, reports, CSV/JSON/HTML, environment) into a timestamped ``.zip``.
   The zip name contains the version and time so published attachments stay
   unique. No files are invented: only files that exist under ``out_dir`` are
   packed.

No GUI/tkinter dependency.
"""

from __future__ import annotations

import os
import zipfile
from typing import Optional

from .._version import APP_NAME_ASCII, __version__ as VERSION
from . import models, utils

#: project-file format marker
DSPROJ_MARKER = "dockstudio.dsproj"
DSPROJ_FORMAT_VERSION = 1


def save_project(cfg: models.RunConfig, path: str, out_dir: Optional[str] = None) -> str:
    """Write a ``.dsproj`` JSON wrapper (full config + version metadata)."""
    doc = {
        "format": DSPROJ_MARKER,
        "format_version": DSPROJ_FORMAT_VERSION,
        "app": APP_NAME_ASCII,
        "version": VERSION,
        "saved_at": utils.now_str(),
        "out_dir": out_dir or cfg.out_dir,
        "config": cfg.to_dict(),
    }
    utils.write_json(doc, path)
    return path


def is_dsproj(path: str) -> bool:
    try:
        d = utils.read_json(path)
        return isinstance(d, dict) and d.get("format") == DSPROJ_MARKER
    except Exception:
        return False


def load_run_config(path: str) -> models.RunConfig:
    """Load a RunConfig from a ``.dsproj`` file or a legacy ``config.json``."""
    if is_dsproj(path):
        d = utils.read_json(path)
        cfg = models.RunConfig.from_dict(d.get("config") or {})
        out = d.get("out_dir")
        if out and not cfg.out_dir:
            cfg.out_dir = out
        return cfg
    return models.RunConfig.load(path)


def default_export_name(out_dir: str, version: str = VERSION) -> str:
    import re
    import time
    ts = time.strftime("%Y%m%d_%H%M%S")
    stem = utils.safe_name(os.path.basename(os.path.abspath(out_dir))) or "DockStudio"
    return f"DockStudio_{stem}_{version.replace('.', '')}_{ts}.zip"


def export_run_zip(out_dir: str, dest_zip: Optional[str] = None,
                   version: str = VERSION) -> str:
    """Pack the whole output directory into one ``.zip`` (reproducible export).

    Files are added under the top-level folder ``<out_basename>/`` so unpacking
    never spills files into the current directory. The zip itself (when placed
    inside ``out_dir``) is excluded to avoid recursion.
    """
    out_dir = os.path.abspath(out_dir)
    if not os.path.isdir(out_dir):
        raise NotADirectoryError(f"output directory not found: {out_dir}")
    if dest_zip is None:
        dest_zip = os.path.join(os.path.dirname(out_dir), default_export_name(out_dir, version))
    dest_zip = os.path.abspath(dest_zip)
    top = os.path.basename(out_dir) or "DockStudio"
    n = 0
    with zipfile.ZipFile(dest_zip, "w", zipfile.ZIP_DEFLATED, allowZip64=True) as zf:
        for dirpath, dirnames, filenames in os.walk(out_dir):
            dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
            for f in sorted(filenames):
                if f.endswith((".pyc", ".tmp")):
                    continue
                full = os.path.join(dirpath, f)
                if os.path.abspath(full) == dest_zip:
                    continue
                rel = os.path.relpath(full, out_dir)
                arc = os.path.join(top, rel).replace(os.sep, "/")
                try:
                    zf.write(full, arc)
                    n += 1
                except OSError:
                    continue
    return dest_zip


def zip_file_count(zip_path: str) -> int:
    try:
        with zipfile.ZipFile(zip_path, "r") as zf:
            return len(zf.namelist())
    except Exception:
        return -1
