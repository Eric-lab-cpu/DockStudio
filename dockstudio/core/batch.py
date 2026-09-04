"""Batch execution with resume, error capture and optional time-slicing (section 7)."""

from __future__ import annotations

import os
import time
from typing import Callable, Dict, List, Optional

from . import docking, models, utils
from .models import project_subdirs

# progress callback signature: cb(dict)
ProgressCb = Optional[Callable[[dict], None]]
LogCb = Optional[Callable[[str], None]]


def pair_dir_name(rec: str, lig: str) -> str:
    return f"{rec}__{lig}"


def _emit(cb, msg: str):
    if cb:
        try:
            cb(msg)
        except Exception:
            pass


def run_docking_stage(
    cfg: models.RunConfig,
    stage: str,  # screening | refine | selfdock
    prepared_receptors: Dict[str, dict],
    ligand_pdbqt_map: Dict[str, str],     # lig_id -> pdbqt path
    cocrystal_pdbqt_map: Dict[str, dict] = None,  # rec -> {pdbqt, lig_id} (selfdock)
    pairs: Optional[List[tuple]] = None,  # (rec, lig); default all receptor x ligand
    time_slice_s: Optional[float] = None,
    progress: ProgressCb = None,
    log: LogCb = None,
    stage_subdir: Optional[str] = None,
    exhaustiveness_override: Optional[int] = None,
    overwrite: bool = False,
) -> dict:
    """Run a stage with resume-by-``done.flag``.

    stage determines the sub-directory and default exhaustiveness:
      screening -> docking/  with cfg.exhaustiveness
      refine    -> refine/   with cfg.refine_exhaustiveness
      selfdock  -> selfdock/ with cfg.refine_exhaustiveness
    """
    out_root = cfg.out_dir
    sub = project_subdirs(out_root)
    if stage_subdir is None:
        stage_subdir = {
            "screening": sub["screening"],
            "refine": sub["refine"],
            "selfdock": os.path.join(out_root, "selfdock"),
        }[stage]
    os.makedirs(stage_subdir, exist_ok=True)

    if pairs is None:
        recs = list(prepared_receptors.keys())
        if stage == "selfdock":
            pairs = [(r, d["lig_id"]) for r, d in (cocrystal_pdbqt_map or {}).items() if r in recs]
        else:
            pairs = [(r, l) for r in recs for l in ligand_pdbqt_map]

    if stage == "selfdock":
        exhaustiveness = cfg.refine_exhaustiveness if cfg.run_refine else cfg.exhaustiveness
    elif stage == "refine":
        exhaustiveness = cfg.refine_exhaustiveness
    else:
        exhaustiveness = cfg.exhaustiveness
    if exhaustiveness_override is not None:
        exhaustiveness = exhaustiveness_override

    todo = list(pairs)
    done = skipped = errors = 0
    error_list = []
    started = time.time()
    remaining = todo[:]

    for idx, (rec, lig) in enumerate(todo, 1):
        flag = os.path.join(stage_subdir, pair_dir_name(rec, lig), "done.flag")
        if os.path.exists(flag) and not overwrite:
            skipped += 1
            remaining.pop(0)
            continue
        if rec not in prepared_receptors:
            error_list.append({"pair": f"{rec}__{lig}", "error": "receptor not prepared"})
            errors += 1
            remaining.pop(0)
            continue
        if stage == "selfdock":
            coc = (cocrystal_pdbqt_map or {}).get(rec, {})
            lig_pdbqt = coc.get("pdbqt") or ligand_pdbqt_map.get(lig)
        else:
            lig_pdbqt = ligand_pdbqt_map.get(lig)
        if not lig_pdbqt:
            error_list.append({"pair": f"{rec}__{lig}", "error": "ligand not prepared"})
            errors += 1
            remaining.pop(0)
            continue
        box = cfg.boxes.get(rec)
        if not box:
            error_list.append({"pair": f"{rec}__{lig}", "error": "no search box for receptor"})
            errors += 1
            remaining.pop(0)
            continue
        try:
            _emit(log, f"[{stage}] docking {rec} x {lig} ({idx}/{len(todo)}) exh={exhaustiveness}")
            docking.dock_pair(
                prepared_receptors[rec]["pdbqt"], lig_pdbqt, stage_subdir,
                rec_name=rec, lig_name=lig,
                center=box["center"], size=box["size"],
                exhaustiveness=exhaustiveness, n_poses=cfg.n_poses,
                energy_range=cfg.energy_range, cpu=max(1, int(cfg.cpu)),
                stage=stage, overwrite=overwrite,
            )
            done += 1
        except Exception as e:
            errors += 1
            error_list.append({"pair": f"{rec}__{lig}", "error": str(e)})
            _emit(log, f"[{stage}] ERROR {rec} x {lig}: {e}")
        remaining.pop(0)
        if progress:
            try:
                progress({"done": done, "skipped": skipped, "errors": errors,
                          "total": len(todo), "current": f"{rec}__{lig}", "stage": stage,
                          "remaining": len(remaining)})
            except Exception:
                pass
        if time_slice_s is not None and (time.time() - started) > time_slice_s:
            _emit(log, f"[{stage}] time-slice reached; {len(remaining)} pairs still pending")
            break

    summary = {"stage": stage, "done": done, "skipped": skipped, "errors": errors,
               "total": len(todo), "remaining": remaining,
               "error_list": error_list}
    return summary


def collect_results(stage_dir: str) -> List[dict]:
    """Read result.json from every completed pair directory."""
    results = []
    if not os.path.isdir(stage_dir):
        return results
    for d in sorted(os.listdir(stage_dir)):
        rj = os.path.join(stage_dir, d, "result.json")
        if os.path.exists(rj):
            try:
                results.append(utils.read_json(rj))
            except Exception:
                continue
    return results
