"""Batch execution with resume, error capture and optional time-slicing (section 7).

v2.0 additions
--------------
* process-pool parallel execution (``workers > 1``) for virtual-screening
  throughput. Every worker is an independent process so the python-vina backend
  stays safe; each dock still uses ``cfg.cpu`` threads, so results are
  reproducible no matter how many ligands run concurrently.
* streaming result iteration (:func:`iter_results`) so ranking / shortlist /
  matrix logic never loads tens of thousands of ``result.json`` documents into
  memory at once.
* run summaries now include wall-clock time and throughput.
"""

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
    workers: int = 1,
) -> dict:
    """Run a stage with resume-by-``done.flag``.

    stage determines the sub-directory and default exhaustiveness:
      screening -> docking/  with cfg.exhaustiveness
      refine    -> refine/   with cfg.refine_exhaustiveness
      selfdock  -> selfdock/ with cfg.refine_exhaustiveness

    When ``workers > 1`` the valid pairs are docked through a process pool
    (each worker is an independent process; ``cfg.cpu`` threads are used per
    dock, so per-pair results do not depend on ``workers``).
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

    if workers and int(workers) > 1:
        return _run_stage_pool(
            cfg, stage, prepared_receptors, ligand_pdbqt_map, cocrystal_pdbqt_map,
            pairs, time_slice_s, progress, log, stage_subdir, exhaustiveness,
            overwrite, int(workers))

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

    elapsed = time.time() - started
    processed = done + skipped + errors
    summary = {"stage": stage, "done": done, "skipped": skipped, "errors": errors,
               "total": len(todo), "remaining": remaining,
               "error_list": error_list,
               "workers": 1, "elapsed_s": round(elapsed, 2),
               "throughput_per_s": round(processed / elapsed, 2) if elapsed else None}
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


def iter_results(stage_dir: str):
    """Yield result.json dicts one at a time (memory-safe for large screens)."""
    if not os.path.isdir(stage_dir):
        return
    for d in sorted(os.listdir(stage_dir)):
        rj = os.path.join(stage_dir, d, "result.json")
        if os.path.exists(rj):
            try:
                yield utils.read_json(rj)
            except Exception:
                continue


# ---------------------------------------------------------------------------
# Process-pool parallel docking (v2.0)
# ---------------------------------------------------------------------------
def _task_args_for(cfg, stage, prepared_receptors, ligand_pdbqt_map,
                   cocrystal_pdbqt_map, rec, lig, stage_subdir, exhaustiveness,
                   overwrite) -> Optional[dict]:
    """Pre-flight checks for one pair -> task-args dict or None if invalid."""
    if rec not in prepared_receptors:
        return {"invalid": f"receptor not prepared: {rec}"}
    if stage == "selfdock":
        coc = (cocrystal_pdbqt_map or {}).get(rec, {})
        lig_pdbqt = coc.get("pdbqt") or ligand_pdbqt_map.get(lig)
    else:
        lig_pdbqt = ligand_pdbqt_map.get(lig)
    if not lig_pdbqt:
        return {"invalid": f"ligand not prepared: {lig}"}
    box = cfg.boxes.get(rec)
    if not box:
        return {"invalid": f"no search box for receptor: {rec}"}
    return {
        "receptor_pdbqt": prepared_receptors[rec]["pdbqt"],
        "ligand_pdbqt": lig_pdbqt,
        "stage_subdir": stage_subdir,
        "rec_name": rec,
        "lig_name": lig,
        "center": list(box["center"]),
        "size": list(box["size"]),
        "exhaustiveness": exhaustiveness,
        "n_poses": cfg.n_poses,
        "energy_range": cfg.energy_range,
        "cpu": max(1, int(cfg.cpu)),
        "stage": stage,
        "overwrite": overwrite,
    }


def _dock_one_task(args: dict) -> dict:
    """Top-level worker entry (must stay importable for process-pool spawn)."""
    rec = args["rec_name"]
    lig = args["lig_name"]
    try:
        res = docking.dock_pair(
            args["receptor_pdbqt"], args["ligand_pdbqt"], args["stage_subdir"],
            rec_name=rec, lig_name=lig,
            center=args["center"], size=args["size"],
            exhaustiveness=args["exhaustiveness"], n_poses=args["n_poses"],
            energy_range=args["energy_range"], cpu=args["cpu"],
            stage=args["stage"], overwrite=args["overwrite"],
        )
        return {"pair": f"{rec}__{lig}", "ok": True, "error": None,
                "affinity_best": res.get("affinity_best")}
    except Exception as e:  # docking.dock_pair already writes error.json
        return {"pair": f"{rec}__{lig}", "ok": False,
                "error": f"{type(e).__name__}: {e}"}


def _run_stage_pool(cfg, stage, prepared_receptors, ligand_pdbqt_map,
                    cocrystal_pdbqt_map, pairs, time_slice_s, progress, log,
                    stage_subdir, exhaustiveness, overwrite, workers) -> dict:
    """Parallel dock many pairs through a process pool (bounded in-flight)."""
    from concurrent.futures import ProcessPoolExecutor

    todo = list(pairs)
    started = time.time()
    done = skipped = errors = 0
    error_list = []
    pending_tasks = []          # valid pairs not yet dispatched
    preflight_errors = []
    for (rec, lig) in todo:
        flag = os.path.join(stage_subdir, pair_dir_name(rec, lig), "done.flag")
        if os.path.exists(flag) and not overwrite:
            skipped += 1
            continue
        t = _task_args_for(cfg, stage, prepared_receptors, ligand_pdbqt_map,
                           cocrystal_pdbqt_map, rec, lig, stage_subdir,
                           exhaustiveness, overwrite)
        if t is None:
            errors += 1
            error_list.append({"pair": f"{rec}__{lig}",
                               "error": "task planning failed"})
        elif "invalid" in t:
            errors += 1
            preflight_errors.append((rec, lig))
            error_list.append({"pair": f"{rec}__{lig}", "error": t["invalid"]})
        else:
            pending_tasks.append((rec, lig, t))

    _emit(log, f"[{stage}] parallel pool: workers={workers}, "
               f"pending={len(pending_tasks)}, skipped={skipped}, "
               f"preflight_errors={len(preflight_errors)}")
    done_pairs = set()
    if pending_tasks:
        from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
        # bound the number of in-flight futures to avoid 100k+ queued pickles
        window = max(workers * 2, workers + 2)
        pending = list(pending_tasks)
        idx = 0
        inflight = {}
        stop_dispatch = False
        executor = ProcessPoolExecutor(max_workers=workers)
        try:
            def _submit_next() -> bool:
                nonlocal idx
                if idx >= len(pending):
                    return False
                rec, lig, t = pending[idx]
                idx += 1
                f = executor.submit(_dock_one_task, t)
                inflight[f] = (rec, lig)
                return True

            # prime the window
            while len(inflight) < window and _submit_next():
                pass
            while inflight:
                finished, _ = wait(list(inflight), return_when=FIRST_COMPLETED)
                for fut in finished:
                    rec, lig = inflight.pop(fut)
                    try:
                        r = fut.result()
                    except Exception as e:  # worker crashed
                        r = {"ok": False, "error": f"worker crashed: {e}",
                             "pair": f"{rec}__{lig}"}
                    if r.get("ok"):
                        done += 1
                    else:
                        errors += 1
                        error_list.append({"pair": r.get("pair") or f"{rec}__{lig}",
                                           "error": r.get("error") or "unknown"})
                        _emit(log, f"[{stage}] ERROR {rec} x {lig}: {r.get('error')}")
                    done_pairs.add((rec, lig))
                    if progress:
                        try:
                            progress({"done": done, "skipped": skipped,
                                      "errors": errors, "total": len(todo),
                                      "current": f"{rec}__{lig}", "stage": stage,
                                      "remaining": len(pending) - len(done_pairs)})
                        except Exception:
                            pass
                # stop dispatching when the time slice expires
                if time_slice_s is not None and (time.time() - started) > time_slice_s:
                    stop_dispatch = True
                # replenish the window
                if not stop_dispatch:
                    while len(inflight) < window and _submit_next():
                        pass
                elif len(inflight) == 0:
                    break
            remaining = [p for p in pending if p[:2] not in done_pairs]
            if stop_dispatch:
                _emit(log, f"[{stage}] time-slice reached; {len(remaining)} "
                           f"pairs still pending")
        finally:
            executor.shutdown(wait=True)
    else:
        remaining = []
    elapsed = time.time() - started
    processed = done + skipped + errors
    summary = {"stage": stage, "done": done, "skipped": skipped, "errors": errors,
               "total": len(todo), "remaining": [(r, l) for r, l in remaining],
               "error_list": error_list,
               "workers": workers, "elapsed_s": round(elapsed, 2),
               "throughput_per_s": round(processed / elapsed, 2) if elapsed else None}
    _emit(log, f"[{stage}] pool done: {done} ok, {skipped} skipped, "
               f"{errors} errors in {elapsed:.1f}s")
    return summary
