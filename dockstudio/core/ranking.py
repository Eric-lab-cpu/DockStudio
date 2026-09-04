"""Ranking, refinement and final Top-K logic (protocol section 8)."""

from __future__ import annotations

import csv
import os
from typing import Dict, List

from . import utils
from .models import RunConfig, project_subdirs


def affinity_matrix(results: List[dict], receptors: List[str], ligands: List[str]) -> Dict[str, dict]:
    """Screening mode-1 affinity per receptor x ligand."""
    rows: Dict[str, dict] = {r: {"receptor": r} for r in receptors}
    for l in ligands:
        for r in receptors:
            rows[r].setdefault(l, None)
    for res in results:
        r = res["receptor"]; l = res["ligand"]
        if r in rows:
            rows[r][l] = res.get("affinity_best")
    return rows


def write_affinity_matrix(matrix: Dict[str, dict], path: str, ligands: List[str]) -> None:
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["receptor"] + ligands)
        for r, row in matrix.items():
            w.writerow([r] + ["" if row.get(l) is None else round(row[l], 2) for l in ligands])


def rank_ligands_for_receptor(matrix: Dict[str, dict], receptor: str) -> List[tuple]:
    row = matrix[receptor]
    ranked = []
    for l, aff in row.items():
        if l == "receptor" or aff is None:
            continue
        ranked.append((l, float(aff)))
    ranked.sort(key=lambda x: x[1])
    return ranked


def choose_refine_pairs(cfg: RunConfig, matrix: Dict[str, dict],
                        shortlist: int = None) -> List[tuple]:
    n = shortlist or cfg.top_refine
    pairs = []
    for rec in cfg.receptor_names:
        for lig, aff in rank_ligands_for_receptor(matrix, rec)[:n]:
            pairs.append((rec, lig))
    return pairs


def write_final_topk(cfg: RunConfig, refine_results: List[dict],
                     screening_matrix: Dict[str, dict],
                     results_dir: str,
                     ligands: List[str]) -> dict:
    """Derive final Top-K from refined (mode-1) energies.

    Writes final_top5.csv and _final_top5.json; returns a container dict.
    """
    best_refine: Dict[tuple, float] = {}
    for res in refine_results:
        best_refine[(res["receptor"], res["ligand"])] = res.get("affinity_best")

    output_rows = []
    top_json = {}
    for rec in cfg.receptor_names:
        scored = []
        for lig in ligands:
            aff = best_refine.get((rec, lig))
            source = "refine"
            if aff is None:
                aff = screening_matrix.get(rec, {}).get(lig)
                source = "screening-fallback"
            if aff is None:
                continue
            scored.append((lig, float(aff), source))
        scored.sort(key=lambda x: x[1])
        k = cfg.top_k
        for lig, aff, source in scored[:k]:
            output_rows.append({"receptor": rec, "ligand": lig,
                                "best_affinity_kcal_mol": round(aff, 2),
                                "source": source})
        top_json[rec] = {"top": [{"ligand": lig, "affinity_kcal_mol": round(aff, 2)}
                                 for lig, aff, source in scored[:k]]}

    csv_path = os.path.join(results_dir, "final_top5.csv")
    utils.write_json(top_json, os.path.join(results_dir, "_final_top5.json"))
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as fh:
        cols = ["receptor", "ligand", "best_affinity_kcal_mol", "source"]
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(output_rows)
    return {"rows": output_rows, "json": top_json, "csv_path": csv_path}
