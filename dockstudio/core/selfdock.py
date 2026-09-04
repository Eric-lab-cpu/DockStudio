"""Self-docking (crystal-ligand redock) validation (section 9)."""

from __future__ import annotations

import csv
import math
import os
from typing import Dict, List, Optional, Tuple

import numpy as np
from rdkit import Chem

from . import docking, models, utils
from .models import project_subdirs


def load_heavy_coords_sdf(sdf_path: str) -> List[Tuple[str, List[float]]]:
    mol = Chem.MolFromMolFile(sdf_path, sanitize=False, removeHs=True)
    if mol is None:
        return []
    conf = mol.GetConformer()
    out = []
    for a in mol.GetAtoms():
        if a.GetSymbol() == "H":
            continue
        p = conf.GetAtomPosition(a.GetIdx())
        out.append((a.GetSymbol(), [p.x, p.y, p.z]))
    return out


def _match_and_rmsd(query: List[Tuple[str, List[float]]],
                    target: List[Tuple[str, List[float]]]) -> Optional[dict]:
    """Greedy nearest-neighbour matching by element (no re-superposition)."""
    tgt = list(target)
    total = 0.0
    used = set()
    pairs = []
    for el, q in query:
        best = None
        bestd = float("inf")
        for j, (tel, tp) in enumerate(tgt):
            if j in used or tel != el:
                continue
            d = math.dist(q, tp)
            if d < bestd:
                bestd = d
                best = j
        if best is None:
            continue
        used.add(best)
        total += bestd ** 2
        pairs.append((bestd, el))
    if not pairs:
        return None
    n = len(pairs)
    return {"rmsd": math.sqrt(total / n), "n_matched": n, "n_query": len(query)}


def pose_rmsd_vs_crystal(pose: docking.DockedPose, crystal_sdf: str) -> Optional[float]:
    crystal = load_heavy_coords_sdf(crystal_sdf)
    if not crystal:
        return None
    pose_heavy = [(a.element, a.coord()) for a in pose.heavy_atoms()]
    res = _match_and_rmsd(crystal, pose_heavy)
    return res["rmsd"] if res else None


def evaluate_selfdock(
    pair_dir: str,
    crystal_sdf: str,
) -> dict:
    """Compute RMSD of docked poses vs crystal coordinates (same frame)."""
    poses_path = os.path.join(pair_dir, "poses.pdbqt")
    if not os.path.exists(poses_path):
        return {"error": "no poses.pdbqt"}
    poses = docking.parse_pdbqt_poses(poses_path)
    crystal = load_heavy_coords_sdf(crystal_sdf)
    out = {"poses_available": len(poses), "mode1_rmsd": None,
           "min_rmsd": None, "min_rmsd_pose": None, "energy_min_rmsd": None,
           "energy_delta_min_rmsd": None, "n_crystal_atoms": len(crystal),
           "pass": False}
    per_pose = []
    for p in poses:
        rmsd = pose_rmsd_vs_crystal(p, crystal_sdf)
        per_pose.append({"pose": p.index, "affinity": round(p.affinity, 2),
                         "rmsd": round(rmsd, 3) if rmsd is not None else None})
        if out["mode1_rmsd"] is None and p.index == 1:
            out["mode1_rmsd"] = round(rmsd, 3) if rmsd is not None else None
        if rmsd is not None and (out["min_rmsd"] is None or rmsd < out["min_rmsd"]):
            out["min_rmsd"] = round(rmsd, 3)
            out["min_rmsd_pose"] = p.index
            out["energy_min_rmsd"] = round(p.affinity, 2)
    if out["mode1_rmsd"] is not None:
        out["pass"] = out["mode1_rmsd"] < 2.0
    if out["energy_min_rmsd"] is not None and out["mode1_rmsd"] is not None and len(poses) > 0:
        out["energy_delta_min_rmsd"] = round(out["energy_min_rmsd"] - poses[0].affinity, 2)
    out["poses"] = per_pose
    return out


def run_selfdock_validation(cfg: models.RunConfig, cocrystal: Dict[str, dict],
                            out_dir: str, progress=None, log=None) -> List[dict]:
    """Evaluate every completed selfdock pair and write summary files."""
    results_dir = os.path.join(out_dir, "results")
    os.makedirs(results_dir, exist_ok=True)
    selfdock_dir = os.path.join(out_dir, "selfdock")
    rows = []
    for rec, meta in cocrystal.items():
        pair = f"{rec}__{meta['lig_id']}"
        pair_dir = os.path.join(selfdock_dir, pair)
        if not os.path.exists(os.path.join(pair_dir, "done.flag")):
            continue
        ev = evaluate_selfdock(pair_dir, meta["crystal_sdf"])
        ev.update({"receptor": rec, "cocrystal_ligand": meta["lig_id"],
                   "cocrystal_resname": meta.get("resname", ""),
                   "box": cfg.boxes.get(rec)})
        rows.append(ev)
        if log:
            log(f"[selfdock] {rec}: mode-1 RMSD = {ev.get('mode1_rmsd')} A "
                f"({'PASS' if ev.get('pass') else 'FAIL'})")
    # write csv + json
    utils.write_json({"selfdock": rows}, os.path.join(results_dir, "selfdock_summary.json"))
    cols = ["receptor", "cocrystal_ligand", "cocrystal_resname", "mode1_rmsd", "min_rmsd",
            "min_rmsd_pose", "energy_min_rmsd", "energy_delta_min_rmsd", "pass", "n_crystal_atoms"]
    with open(os.path.join(results_dir, "selfdock_summary.csv"), "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return rows
