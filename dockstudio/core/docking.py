"""AutoDock Vina docking engine (protocol section 6).

Two interchangeable backends are supported:
  * ``python`` - the ``vina`` Python package (Linux/macOS wheels)
  * ``cli``    - the ``vina`` command-line executable (used automatically on
                 platforms without python bindings, e.g. Windows; the binary is
                 located via $DOCKSTUDIO_VINA / $VINA_BIN / PATH / app layout)

Each receptor-ligand pair writes a self-contained directory:

    result.json   - machine readable result (seed, exhaustiveness, box,
                    affinities, backend, files)
    poses.pdbqt   - all requested poses (MODEL blocks)
    best.pdbqt    - first (lowest-energy) pose
    done.flag     - written only after successful completion (resume marker)
    error.json    - written on failure
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import List, Optional

from . import utils

MODEL_RE = re.compile(r"^MODEL\s+(\d+)")
ENDMDL_RE = re.compile(r"^ENDMDL")
VINA_RESULT_RE = re.compile(r"REMARK\s+VINA RESULT:\s+(-?[\d.]+)")


@dataclass
class PoseAtom:
    name: str
    element: str
    x: float
    y: float
    z: float

    def coord(self) -> List[float]:
        return [self.x, self.y, self.z]


@dataclass
class DockedPose:
    index: int
    affinity: float
    atoms: List[PoseAtom] = field(default_factory=list)

    def heavy_atoms(self) -> List[PoseAtom]:
        return [a for a in self.atoms if a.element != "H"]


def _element_from_name(name: str) -> str:
    # PDBQT names may be like 'C1', 'N2', 'O', 'Cl', 'Br'
    m = re.match(r"^([A-Z][a-z]?|[A-Z])", name)
    if not m:
        return "C"
    sym = m.group(1)
    if sym in ("C", "N", "O", "S", "P", "H", "F", "Cl", "Br", "I", "B", "Se"):
        return sym
    return sym[0]


# AutoDock (PDBQT) atom-type -> real chemical element.  PDBQT atom lines end
# with an AutoDock *type* token (not necessarily an element symbol):
#   C  -> carbon, A -> aromatic carbon, NA -> hydrogen-bond acceptor N,
#   OA -> acceptor O, SA -> acceptor S, HD -> hydrogen-donor H, ...
# Treating "A"/"NA"/"OA"/"SA"/"HD" as elements would silently corrupt
# self-docking RMSD element-matching, so we map them to real elements.
_AUTODOCK_TYPE_ELEMENT = {
    "C": "C", "A": "C",
    "N": "N", "NA": "N",
    "O": "O", "OA": "O",
    "S": "S", "SA": "S",
    "H": "H", "HD": "H",
    "P": "P", "F": "F", "CL": "Cl", "BR": "Br", "I": "I",
    "B": "B", "SI": "Si", "SE": "Se", "MG": "Mg", "MN": "Mn",
    "ZN": "Zn", "CA": "Ca", "FE": "Fe", "CO": "Co", "CU": "Cu",
    "MO": "Mo", "CD": "Cd", "HG": "Hg", "NI": "Ni", "PT": "Pt",
}


def _element_of_pdbqt_atom(line: str, name: str) -> str:
    """Robust element detection for a PDBQT ATOM/HETATM line.

    Prefers a real element symbol if one is present, then the AutoDock type
    mapping, then a name-based fallback.
    """
    toks = line.split()
    if toks:
        last = toks[-1]
        # A valid element symbol at the end wins (Vina may write element here).
        if re.fullmatch(r"(C|N|O|S|P|H|F|Cl|Br|I|B|Se|Si|Mg|Mn|Zn|Ca|Fe|Co|Cu|Mo)", last, re.I):
            return last[:1].upper() + last[1:].lower()
        mapped = _AUTODOCK_TYPE_ELEMENT.get(last.upper())
        if mapped:
            return mapped
    # columns 77-79 (1-based 77) may hold an element or short type code
    el = line[76:79].strip()
    if len(el) > 1 and el[0] == " ":
        el = el[1:]
    if el and not el.isdigit():
        mapped = _AUTODOCK_TYPE_ELEMENT.get(el.upper())
        if mapped:
            return mapped
        if re.fullmatch(r"(C|N|O|S|P|H|F|Cl|Br|I|B|Se)", el, re.I):
            return el[:1].upper() + el[1:].lower()
    return _element_from_name(name)


def parse_pdbqt_poses(path: str) -> List[DockedPose]:
    """Parse a multi-MODEL PDBQT file written by vina into DockedPose objects."""
    poses: List[DockedPose] = []
    cur: Optional[DockedPose] = None
    for line in open(path, "r", errors="replace"):
        s = line.strip()
        mm = MODEL_RE.match(s)
        if mm:
            cur = DockedPose(index=int(mm.group(1)), affinity=float("nan"))
            poses.append(cur)
            continue
        if cur is None:
            continue
        if ENDMDL_RE.match(s):
            cur = None
            continue
        vr = VINA_RESULT_RE.search(s)
        if vr and cur and cur.affinity != cur.affinity:  # nan
            cur.affinity = float(vr.group(1))
            continue
        if s.startswith(("ATOM", "HETATM")):
            name = line[12:16].strip()
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            # element: real element / AutoDock type token / name-based fallback
            el = _element_of_pdbqt_atom(line, name)
            cur.atoms.append(PoseAtom(name=name, element=el, x=x, y=y, z=z))
    return poses


def _models_block_text(pdbqt_path: str, model_numbers: Optional[list] = None) -> str:
    """Extract MODEL block text (model numbers optional)."""
    blocks = []
    cur = []
    in_model = False
    for line in open(pdbqt_path, errors="replace"):
        if MODEL_RE.match(line):
            in_model = True
            cur = [line]
            continue
        if ENDMDL_RE.match(line) and in_model:
            cur.append(line)
            blocks.append(cur)
            in_model = False
            continue
        if in_model:
            cur.append(line)
    if model_numbers is None:
        return "".join("".join(b) for b in blocks)
    keep = set(model_numbers)
    return "".join("".join(b) for b in blocks if int(MODEL_RE.match(b[0]).group(1)) in keep)


def _find_vina_cli() -> Optional[str]:
    """Locate a Vina command-line executable (used on platforms without the
    python-vina bindings, e.g. Windows). Search order:
      1) $DOCKSTUDIO_VINA / $VINA_BIN
      2) 'vina' on PATH
      3) <app_root>/vina/vina[.exe]  (installer layout)
    """
    for key in ("DOCKSTUDIO_VINA", "VINA_BIN"):
        env = os.environ.get(key)
        if env and os.path.exists(env):
            return env
    p = utils.which("vina", "vina.exe")
    if p:
        return p
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    for cand in (os.path.join(here, "vina", "vina.exe"),
                 os.path.join(here, "vina", "vina"),
                 os.path.join(os.path.dirname(os.path.dirname(here)), "vina", "vina.exe")):
        if os.path.exists(cand):
            return cand
    return None


def _python_vina_available() -> bool:
    try:
        import vina  # noqa: F401
        return True
    except Exception:
        return False


def resolve_backend(force_cli: bool = False) -> str:
    env_cli = os.environ.get("DOCKSTUDIO_VINA") or os.environ.get("VINA_BIN")
    if force_cli or env_cli or os.environ.get("DOCKSTUDIO_FORCE_CLI") == "1":
        if _find_vina_cli():
            return "cli"
        raise RuntimeError("DOCKSTUDIO_VINA set but no Vina CLI executable found")
    if _python_vina_available():
        return "python"
    if _find_vina_cli():
        return "cli"
    raise RuntimeError(
        "no Vina engine available: python 'vina' package not importable and no "
        "vina CLI executable found (set DOCKSTUDIO_VINA to the vina binary path).")


def _write_result_files(pair_dir, pair, rec_name, lig_name, stage, seed, exhaustiveness,
                        n_poses, energy_range, center, size, poses_pdbqt,
                        affinities, backend, vina_version="1.2.x") -> dict:
    """Shared file/result writing for python and CLI backends."""
    poses = parse_pdbqt_poses(poses_pdbqt)
    if not affinities:
        affinities = [p.affinity for p in poses]
    best_text = _models_block_text(poses_pdbqt, [1])
    best_pdbqt = os.path.join(pair_dir, "best.pdbqt")
    with open(best_pdbqt, "w", encoding="utf-8") as fh:
        fh.write(best_text)
    result = {
        "pair": pair, "receptor": rec_name, "ligand": lig_name, "stage": stage,
        "seed": seed, "exhaustiveness": exhaustiveness, "n_poses": n_poses,
        "energy_range": energy_range,
        "box_center": [round(c, 3) for c in center],
        "box_size": [round(s, 3) for s in size],
        "affinities": affinities,
        "affinity_best": affinities[0] if affinities else None,
        "backend": backend,
        "vina_version": vina_version,
        "files": {"poses_pdbqt": os.path.basename(poses_pdbqt),
                  "best_pdbqt": os.path.basename(best_pdbqt)},
        "completed_at": utils.now_str(),
    }
    utils.write_json(result, os.path.join(pair_dir, "result.json"))
    with open(os.path.join(pair_dir, "done.flag"), "w", encoding="utf-8") as fh:
        fh.write(utils.now_str() + "\n")
    return result


def _dock_with_python(receptor_pdbqt, ligand_pdbqt, pair_dir, pair, rec_name, lig_name,
                      center, size, exhaustiveness, n_poses, energy_range, cpu, seed,
                      stage) -> dict:
    from vina import Vina
    v = Vina(sf_name="vina", cpu=cpu, seed=seed, verbosity=0)
    v.set_receptor(rigid_pdbqt_filename=receptor_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)
    v.compute_vina_maps(center=center, box_size=size)
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    poses_pdbqt = os.path.join(pair_dir, "poses.pdbqt")
    v.write_poses(poses_pdbqt, n_poses=n_poses, energy_range=energy_range, overwrite=True)
    return _write_result_files(pair_dir, pair, rec_name, lig_name, stage, seed,
                               exhaustiveness, n_poses, energy_range, center, size,
                               poses_pdbqt, [], backend="python")


def _dock_with_cli(vina_bin, receptor_pdbqt, ligand_pdbqt, pair_dir, pair, rec_name,
                   lig_name, center, size, exhaustiveness, n_poses, energy_range,
                   cpu, seed, stage) -> dict:
    conf = os.path.join(pair_dir, "vina.conf.txt")
    poses_pdbqt = os.path.join(pair_dir, "poses.pdbqt")
    log_txt = os.path.join(pair_dir, "vina.log")
    # Align box sizes to the Vina grid spacing (0.375 A) so the number of
    # grid points per axis is an integer (avoids .map / voxel-count errors).
    spacing = 0.375
    sz = [max(6.0, round(s / spacing) * spacing) for s in size]
    lines = [
        f"receptor = {receptor_pdbqt}",
        f"ligand = {ligand_pdbqt}",
        f"center_x = {center[0]}", f"center_y = {center[1]}", f"center_z = {center[2]}",
        f"size_x = {sz[0]}", f"size_y = {sz[1]}", f"size_z = {sz[2]}",
        f"exhaustiveness = {exhaustiveness}",
        f"num_modes = {n_poses}",
        f"energy_range = {energy_range}",
        f"cpu = {max(1, int(cpu))}",
        f"seed = {seed}",
        f"out = {poses_pdbqt}",
    ]
    with open(conf, "w", encoding="ascii") as fh:
        fh.write("\n".join(lines) + "\n")

    def _run(args):
        return utils.run_cmd(args, timeout_s=7200, cwd=pair_dir)

    rc, out, err = _run([vina_bin, "--config", conf])
    method = "config"
    if rc != 0:
        # fallback: pass everything directly on the command line
        direct = [
            vina_bin,
            "--receptor", receptor_pdbqt, "--ligand", ligand_pdbqt,
            "--center_x", str(center[0]), "--center_y", str(center[1]),
            "--center_z", str(center[2]),
            "--size_x", str(sz[0]), "--size_y", str(sz[1]), "--size_z", str(sz[2]),
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes", str(n_poses), "--energy_range", str(energy_range),
            "--cpu", str(max(1, int(cpu))), "--seed", str(seed),
            "--out", poses_pdbqt,
        ]
        rc2, out2, err2 = _run(direct)
        if rc2 == 0:
            rc, out, err = rc2, out2, err2
            method = "direct"
        else:
            rc, out, err = rc2, out2, err2
            method = "direct(also failed)"
    if rc != 0:
        detail = ""
        if os.path.exists(log_txt):
            detail = "\n[log]\n" + open(log_txt, encoding="utf-8",
                                        errors="replace").read()
        combined = (out or "") + "\n" + (err or "")
        lines = [ln for ln in combined.splitlines() if ln.strip()]
        head = "\n".join(lines[:20]) if lines else "(no output)"
        tail = "\n".join(lines[-10:]) if len(lines) > 20 else ""
        try:
            with open(os.path.join(pair_dir, "vina_stderr.txt"), "w",
                      encoding="utf-8", errors="replace") as fh:
                fh.write(combined)
        except Exception:
            pass
        raise RuntimeError(
            f"vina CLI failed ({method}, rc={rc})\n--- head ---\n{head}\n"
            f"--- tail ---\n{tail}\n{detail[-2000:]}")
    if not os.path.exists(poses_pdbqt):
        raise RuntimeError("vina CLI did not produce an output poses file")
    try:
        with open(log_txt, "w", encoding="utf-8", errors="replace") as fh:
            fh.write((out or "") + "\n" + (err or ""))
    except Exception:
        pass
    aff = [p.affinity for p in parse_pdbqt_poses(poses_pdbqt)]
    return _write_result_files(pair_dir, pair, rec_name, lig_name, stage, seed,
                               exhaustiveness, n_poses, energy_range, center, sz,
                               poses_pdbqt, aff, backend="cli")



def dock_pair(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    out_dir: str,
    rec_name: str,
    lig_name: str,
    center: List[float],
    size: List[float],
    exhaustiveness: int = 16,
    n_poses: int = 9,
    energy_range: float = 4.0,
    cpu: int = 1,
    seed: Optional[int] = None,
    stage: str = "docking",
    overwrite: bool = False,
) -> dict:
    """Dock one ligand into one receptor; write result files. Returns result dict."""
    pair = f"{rec_name}__{lig_name}"
    pair_dir = os.path.join(out_dir, pair)
    os.makedirs(pair_dir, exist_ok=True)
    done_flag = os.path.join(pair_dir, "done.flag")
    result_json = os.path.join(pair_dir, "result.json")

    if os.path.exists(done_flag) and not overwrite:
        try:
            return utils.read_json(result_json)
        except Exception:
            pass

    seed = seed if seed is not None else utils.crc32_seed(rec_name, lig_name)
    backend = resolve_backend()
    try:
        if backend == "python":
            return _dock_with_python(
                receptor_pdbqt, ligand_pdbqt, pair_dir, pair, rec_name, lig_name,
                center, size, exhaustiveness, n_poses, energy_range, cpu, seed, stage)
        vina_bin = _find_vina_cli()
        if not vina_bin:
            raise RuntimeError("Vina CLI path not resolved")
        return _dock_with_cli(
            vina_bin, receptor_pdbqt, ligand_pdbqt, pair_dir, pair, rec_name, lig_name,
            center, size, exhaustiveness, n_poses, energy_range, cpu, seed, stage)
    except Exception as e:
        err = {"pair": pair, "stage": stage,
               "error": f"{type(e).__name__}: {e}", "backend": backend, "log_tail": ""}
        utils.write_json(err, os.path.join(pair_dir, "error.json"))
        raise

