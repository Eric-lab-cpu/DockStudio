"""Automated QC cross-checks (protocol sections 8.5 and 12)."""

from __future__ import annotations

import csv
import json
import os
from typing import Dict, List

from . import batch, models, utils
from .models import project_subdirs


def count_done_flags(stage_dir: str) -> int:
    if not os.path.isdir(stage_dir):
        return 0
    return sum(1 for d in os.listdir(stage_dir)
               if os.path.isfile(os.path.join(stage_dir, d, "done.flag")))


def check_matrix_vs_results(matrix: Dict[str, dict], results: List[dict], tol: float = 0.02) -> List[str]:
    issues = []
    for res in results:
        r = res["receptor"]; l = res["ligand"]
        val = res.get("affinity_best")
        if r not in matrix or l not in matrix[r]:
            issues.append(f"result {r}__{l} missing from matrix")
            continue
        mval = matrix[r][l]
        if val is None or mval is None:
            continue
        if abs(float(val) - float(mval)) > tol:
            issues.append(f"matrix mismatch for {r}__{l}: result={val} matrix={mval}")
    return issues


def check_csv_vs_json(csv_path: str, json_path: str) -> List[str]:
    issues = []
    if not (os.path.exists(csv_path) and os.path.exists(json_path)):
        return [f"missing {os.path.basename(csv_path)} or {os.path.basename(json_path)}"]
    with open(json_path, encoding="utf-8") as fh:
        j = json.load(fh)
    with open(csv_path, encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        rec = row["receptor"]; lig = row["ligand"]
        top = j.get(rec, {}).get("top", [])
        if not any(t["ligand"] == lig for t in top):
            issues.append(f"CSV row {rec}__{lig} missing from {os.path.basename(json_path)}")
    for rec, val in j.items():
        if len(val.get("top", [])) != sum(1 for r in rows if r["receptor"] == rec):
            issues.append(f"{os.path.basename(json_path)} receptor {rec} count mismatch vs CSV")
    return issues


def run_qc(cfg: models.RunConfig, out_dir: str) -> List[dict]:
    """Execute all automated consistency checks; returns list of {check, ok, detail}."""
    sub = project_subdirs(out_dir)
    checks = []
    screening_dir = sub["screening"]
    refine_dir = sub["refine"]
    selfdock_dir = os.path.join(out_dir, "selfdock")

    n_screen = count_done_flags(screening_dir)
    n_refine = count_done_flags(refine_dir)
    n_selfdock = count_done_flags(selfdock_dir)
    checks.append({"check": "all screening done.flags present",
                   "ok": n_screen >= len(cfg.receptors) * len([l for l in _ligand_names(cfg)]),
                   "detail": f"screening done flags={n_screen}"})

    # matrix vs screening results
    results = batch.collect_results(screening_dir)
    if results:
        recs = list({r["receptor"] for r in results})
        ligs = list({r["ligand"] for r in results})
        mat = {}
        for r in recs:
            mat[r] = {}
        for res in results:
            mat.setdefault(res["receptor"], {})[res["ligand"]] = res.get("affinity_best")
        iss = check_matrix_vs_results(mat, results)
        checks.append({"check": "screening matrix vs result.json", "ok": not iss,
                       "detail": "; ".join(iss[:5]) or "consistent"})
    # refinement done flags vs the expected shortlist size (per receptor,
    # top_refine ligands per scored set), computed from the screening matrix.
    if results:
        expected_refine = 0
        for rec in cfg.receptor_names:
            row = mat.get(rec, {})
            n_scored = sum(1 for v in row.values() if v is not None)
            expected_refine += min(max(1, int(cfg.top_refine)), n_scored)
        if cfg.run_refine is False:
            expected_refine = 0
        ok_refine = (n_refine >= expected_refine) if expected_refine else True
        detail = f"refine done flags={n_refine} (expected>={expected_refine})"
        if expected_refine and n_refine < expected_refine:
            detail += f"; missing {expected_refine - n_refine}"
        checks.append({"check": "refinement done.flags present", "ok": ok_refine,
                       "detail": detail})
    k = max(1, int(cfg.top_k))
    csv_path = os.path.join(out_dir, "results", f"final_top{k}.csv")
    json_path = os.path.join(out_dir, "results", f"_final_top{k}.json")
    iss2 = check_csv_vs_json(csv_path, json_path)
    checks.append({"check": f"final_top{k}.csv vs _final_top{k}.json", "ok": not iss2,
                   "detail": "; ".join(iss2[:5]) or "consistent"})
    return checks


def _ligand_names(cfg: models.RunConfig) -> List[str]:
    out = []
    for l in cfg.ligands:
        nm = l["name"] or os.path.splitext(os.path.basename(l["path"]))[0]
        out.append(nm)
    return out
