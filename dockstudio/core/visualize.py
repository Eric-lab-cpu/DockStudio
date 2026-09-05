"""PyMOL-based publication rendering (section 11).

Files produced per complex under results/:
  3D_poses/<REC>__<LIG>.png        - 3D binding mode render
  composite/<REC>__<LIG>.png       - 3D render + 2D structure side by side
  pymol_data/<REC>__<LIG>/*.pml    - pure-ASCII PyMOL scripts (+ local data)
  pse/<REC>__<LIG>.pse             - self-contained PyMOL session (optional)

The .pml scripts only ever reference local files with ASCII names inside
their own directory; non-ASCII characters never appear in the script.
"""

from __future__ import annotations

import os
from typing import Dict, List, Optional

from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdMolDescriptors

from . import interactions as plip_mod
from .utils import which

LIG_CHAIN = "Z"


def _find_pymol():
    """Return a command prefix (list[str]) that runs PyMOL, or None.

    Search order: ``pymol`` executable on PATH; else a ``pymol`` module in the
    current interpreter (run as ``python -m pymol``); else None.
    """
    import shutil
    exe = which("pymol")
    if exe:
        return [exe]
    try:
        import pymol  # noqa: F401
        return [shutil.which("python") or "python", "-m", "pymol"]
    except Exception:
        return None


def write_complex_pml(workdir: str, name: str, ligand_chain: str = LIG_CHAIN,
                      include_pse: bool = True) -> str:
    """Write an ASCII PyMOL script into ``workdir`` and return its path."""
    pml_path = os.path.join(workdir, "scene.pml")
    lines = [
        "# DockStudio PyMOL scene (pure ASCII)",
        "bg_color white",
        "set ray_opaque_background, 1",
        "set max_threads, 0",
        "load complex.pdb",
        "hide everything",
        "show cartoon, polymer",
        "color gray70, polymer",
        f"select lig, chain {ligand_chain}",
        "show sticks, lig",
        "util.cnc lig",
        "set stick_radius, 0.18, lig",
        "set two_sided_lighting, 1",
        "set specular, 0.3",
        # pocket residues within 5.5 A of ligand
        "select pocket, (polymer within 5.5 of lig)",
        "show sticks, pocket",
        "util.cnc pocket",
        # transparent surface around binding site
        "select sur, (polymer within 6.5 of lig)",
        "show surface, sur",
        "set surface_color, gray80, sur",
        "set transparency, 0.75, sur",
        "orient lig",
        "zoom lig, 12",
        "set ray_shadows, 1",
        "ray 1600, 1200",
        "png scene_3d.png, dpi=300",
    ]
    if include_pse:
        lines += ["", f"save {name}.pse"]
    with open(pml_path, "w", encoding="ascii") as fh:
        fh.write("\n".join(lines) + "\n")
    return pml_path


def render_complex(workdir: str, name: str, out_png: str, include_pse: bool = True) -> dict:
    """Run PyMOL on the scene in workdir; returns stats."""
    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
    pml_path = write_complex_pml(workdir, name, include_pse=include_pse)
    pymol = _find_pymol()
    stats = {"pml": pml_path, "rc": None, "log_tail": "", "png": out_png, "ok": False}
    if not pymol:
        stats["log_tail"] = "PyMOL executable/module not found; 3D ray render skipped."
        return stats
    from . import utils
    rc, out, err = utils.run_cmd(pymol + ["-cq", pml_path], timeout_s=900, cwd=workdir)
    stats["rc"] = rc
    stats["log_tail"] = (out + err)[-1500:]
    src = os.path.join(workdir, "scene_3d.png")
    if os.path.exists(src):
        os.replace(src, out_png)
        stats["ok"] = True
    pse = os.path.join(workdir, f"{name}.pse")
    if include_pse and os.path.exists(pse):
        stats["pse"] = pse
    return stats


def ligand_2d_png(sdf_path: str, out_png: str) -> bool:
    try:
        mol = Chem.MolFromMolFile(sdf_path, sanitize=True, removeHs=True)
        if mol is None:
            return False
        AllChem.Compute2DCoords(mol)
        d = Draw.MolToImage(mol, size=(600, 450))
        d.save(out_png)
        return True
    except Exception:
        return False


def composite_png(png3d: str, png2d: str, out_png: str) -> None:
    a = Image.open(png3d).convert("RGB")
    b = Image.open(png2d).convert("RGB")
    # match heights
    h = min(a.height, b.height)
    a = a.resize((int(a.width * h / a.height), h))
    b = b.resize((int(b.width * h / b.height), h))
    w = a.width + b.width + 10
    canvas = Image.new("RGB", (w, h), "white")
    canvas.paste(a, (0, 0))
    canvas.paste(b, (a.width + 10, 0))
    canvas.save(out_png)


def tile_pngs(png_paths: List[str], out_png: str, cols: int = 3) -> None:
    if not png_paths:
        return
    imgs = [Image.open(p).convert("RGB") for p in png_paths]
    th = min(i.height for i in imgs)
    imgs = [i.resize((int(i.width * th / i.height), th)) for i in imgs]
    tw = max(i.width for i in imgs)
    rows = (len(imgs) + cols - 1) // cols
    canvas = Image.new("RGB", (tw * cols + (cols + 1) * 5, th * rows + (rows + 1) * 5), "white")
    for k, img in enumerate(imgs):
        r, c = divmod(k, cols)
        canvas.paste(img, (5 + c * (tw + 5), 5 + r * (th + 5)))
    canvas.save(out_png)
