"""PLIP interaction analysis (section 10)."""

from __future__ import annotations

import csv
import os
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional

from . import docking, utils
from .utils import which

TYPE_MAP = {
    "hydrophobic_interactions": "hydrophobic",
    "hydrogen_bonds": "hydrogen_bond",
    "water_bridges": "water_bridge",
    "salt_bridges": "salt_bridge",
    "pi_stacks": "pi_stack",
    "pi_cation_interactions": "pi_cation",
    "halogen_bonds": "halogen_bond",
    "metal_complexes": "metal_complex",
}

DETAIL_TAGS = ["sidechain", "protispos", "donortype", "acceptortype", "dist_h-a",
               "dist_d-a", "don_angle", "lig_group", "protispos", "ligcarbonidx",
               "protcarbonidx", "lig_idx_list", "prot_idx_list"]


def _write_complex_pdb(receptor_pdb: str, best_pdbqt: str, out_path: str) -> None:
    """Merge a receptor PDB and a Vina best pose PDBQT into one PDB file."""
    lines = []
    serial = 1
    for line in open(receptor_pdb, errors="replace"):
        if line.startswith("ATOM"):
            lines.append(line[:6] + f"{serial:5d}" + line[11:])
            serial += 1
    poses = docking.parse_pdbqt_poses(best_pdbqt)
    if poses:
        for a in poses[0].atoms:
            nl = (f"HETATM{serial:5d} {a.name:>4s} LIG Z 161    "
                  f"{a.x:8.3f}{a.y:8.3f}{a.z:8.3f}{1.00:6.2f}{0.00:6.2f}          {a.element:>2s}")
            lines.append(nl + "\n")
            serial += 1
    with open(out_path, "w", encoding="utf-8") as fh:
        for l in lines:
            fh.write(l if l.endswith("\n") else l + "\n")
        fh.write("END\n")


def _parse_plip_xml(path: str) -> List[dict]:
    rows = []
    try:
        root = ET.parse(path).getroot()
    except Exception:
        return rows
    for bs in root.findall("bindingsite"):
        inter = bs.find("interactions")
        if inter is None:
            continue
        for group_tag, label in TYPE_MAP.items():
            group = inter.find(group_tag)
            if group is None:
                continue
            for it in group.findall(label):
                d = {"type": label}
                for sub in it:
                    d[sub.tag] = (sub.text or "").strip()
                rows.append(d)
    return rows


def _residue_label(d: dict) -> str:
    return f"{d.get('restype', '?')}{d.get('resnr', '?')}{d.get('reschain', '')}"


def run_plip_on_complex(complex_pdb: str, out_dir: str, name: str) -> dict:
    """Run PLIP for one complex; returns parsed interaction dict or error."""
    os.makedirs(out_dir, exist_ok=True)
    plip = which("plip")
    if not plip:
        return {"error": "plip executable not found on PATH"}
    name_xml = os.path.join(out_dir, f"{name}.xml")
    try:
        utils.run_cmd([plip, "-f", complex_pdb, "-o", out_dir, "--name", name,
                       "-q", "-x"], timeout_s=600)
    except Exception as e:
        return {"error": f"plip run failed: {e}"}
    if not os.path.exists(name_xml):
        # PLIP writes report.xml when --name handling differs; try default
        alt = os.path.join(out_dir, "report.xml")
        if os.path.exists(alt):
            name_xml = alt
        else:
            return {"error": "plip produced no report"}
    rows = _parse_plip_xml(name_xml)
    if not rows:
        return {"interactions": [], "summary": {}}
    summary = {}
    for r in rows:
        summary[r["type"]] = summary.get(r["type"], 0) + 1
    # rename into generic fields
    flat = []
    for r in rows:
        flat.append({
            "type": r["type"],
            "residue": _residue_label(r),
            "distance_A": r.get("dist", "") or r.get("dist_h-a", ""),
            "details": "; ".join(f"{k}={r[k]}" for k in DETAIL_TAGS if r.get(k)),
        })
    return {"interactions": flat, "summary": summary, "xml": name_xml}


def run_interactions_for_topk(cfg, receptor_clean_map: Dict[str, str],
                              top_pairs: List[dict], out_dir: str,
                              log=None, progress=None) -> dict:
    """Run PLIP on the final top-K complexes and write summary files."""
    results_dir = os.path.join(out_dir, "results", "plip")
    complexes_dir = os.path.join(out_dir, "results", "complexes")
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(complexes_dir, exist_ok=True)
    all_rows = []
    summaries = {}
    for item in top_pairs:
        rec = item["receptor"]; lig = item["ligand"]
        name = f"{rec}__{lig}"
        receptor_pdb = receptor_clean_map.get(rec)
        best_pdbqt = item.get("best_pdbqt")
        if not receptor_pdb or not best_pdbqt or not os.path.exists(best_pdbqt):
            summaries[name] = {"error": "missing receptor or pose file"}
            continue
        complex_pdb = os.path.join(complexes_dir, f"{name}.pdb")
        try:
            _write_complex_pdb(receptor_pdb, best_pdbqt, complex_pdb)
        except Exception as e:
            summaries[name] = {"error": f"complex build failed: {e}"}
            continue
        if log:
            log(f"[plip] analyzing {name}")
        res = run_plip_on_complex(complex_pdb, results_dir, name)
        summaries[name] = res
        for it in res.get("interactions", []):
            all_rows.append({"complex": name, "interaction_type": it["type"],
                             "residue": it["residue"], "distance_A": it["distance_A"],
                             "details": it["details"]})
        if progress:
            progress({"current": name, "stage": "plip"})
    # write flat CSV
    cols = ["complex", "interaction_type", "residue", "distance_A", "details"]
    with open(os.path.join(results_dir, "plip_interactions_all.csv"), "w", newline="",
              encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(all_rows)
    # summary json with counts merged
    merged = {}
    for name, res in summaries.items():
        merged[name] = {"counts": res.get("summary", {}),
                        "n_interactions": len(res.get("interactions", [])),
                        "error": res.get("error")}
        if res.get("interactions"):
            merged[name]["interactions"] = res["interactions"]
    utils.write_json(merged, os.path.join(results_dir, "plip_summary.json"))
    return {"rows": all_rows, "summary": merged}
