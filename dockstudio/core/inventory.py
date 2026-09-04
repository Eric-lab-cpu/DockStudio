"""Input inventory and CSV export (protocol section 4)."""

from __future__ import annotations

import csv
import os
from typing import Dict, List

from . import structure as st
from .models import RunConfig
from .ligand import LigandEntry, iter_molecules, analyze_molecule
from .utils import md5_file


def inventory_receptor(pdb_path: str, name: str, selected_chains: List[str]) -> dict:
    inv = st.inventory_pdb(pdb_path, name=name)
    rows = []
    for chain, ci in inv.chains.items():
        # co-crystal ligand candidates in this chain
        for lig in ci.het_ligands:
            rows.append({
                "receptor": name,
                "file": os.path.basename(pdb_path),
                "md5": inv.md5,
                "chain": chain,
                "n_polymer_res": ci.polymer_res,
                "res_range": f"{ci.res_range[0]}-{ci.res_range[1]}",
                "ligand": lig["resname"],
                "ligand_resseq": lig["resseq"],
                "ligand_n_heavy": lig["n_heavy"],
                "waters": ci.waters,
                "modified": ";".join(f"{m['resname']}{m['resseq']}" for m in ci.modified_residues),
                "altloc_residues": ";".join(ci.altloc_residues),
                "selected_for_docking": "yes" if chain in selected_chains else "no",
            })
    # rows for chains with no cocrystal ligand still summarized once
    seen = set()
    for chain, ci in inv.chains.items():
        if not ci.het_ligands and chain not in seen:
            seen.add(chain)
            rows.append({
                "receptor": name, "file": os.path.basename(pdb_path), "md5": inv.md5,
                "chain": chain, "n_polymer_res": ci.polymer_res,
                "res_range": f"{ci.res_range[0]}-{ci.res_range[1]}",
                "ligand": "", "ligand_resseq": "", "ligand_n_heavy": "",
                "waters": ci.waters,
                "modified": ";".join(f"{m['resname']}{m['resseq']}" for m in ci.modified_residues),
                "altloc_residues": ";".join(ci.altloc_residues),
                "selected_for_docking": "yes" if chain in selected_chains else "no",
            })
    return rows


def inventory_ligand_sdf(sdf_path: str, name: str) -> List[dict]:
    out = []
    for i, mol in iter_molecules(sdf_path):
        nm = name if name else os.path.splitext(os.path.basename(sdf_path))[0]
        if mol is None:
            out.append({"ligand_id": f"{nm}_mol{i}", "source_file": os.path.basename(sdf_path),
                        "index_in_sdf": i, "formula": "", "MW": "", "formal_charge": "",
                        "rotatable_bonds": "", "heavy_atoms": "", "fragments": "",
                        "has_3D": "", "sanitized": "no", "warnings": "unparseable"})
            continue
        entry = analyze_molecule(mol, f"{nm}_mol{i}", sdf_path, i)
        d = entry.to_row()
        d["warnings"] = "; ".join(entry.warnings)
        out.append(d)
    return out


def write_csv(rows: List[dict], path: str) -> None:
    if not rows:
        rows = [{"note": "no data"}]
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    cols = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)


def run_inventory(cfg: RunConfig, out_dir: str) -> dict:
    """Write receptor_inventory.csv and ligand_list.csv; returns inventories."""
    os.makedirs(out_dir, exist_ok=True)
    rec_rows: List[dict] = []
    lig_rows: List[dict] = []
    for rec in cfg.receptors:
        r = rec
        rec_rows += inventory_receptor(r["path"], r["name"], r.get("chains", []))
    for lig in cfg.ligands:
        lig_rows += inventory_ligand_sdf(lig["path"], lig["name"] or None)
    write_csv(rec_rows, os.path.join(out_dir, "receptor_inventory.csv"))
    write_csv(lig_rows, os.path.join(out_dir, "ligand_list.csv"))
    return {"receptor_rows": rec_rows, "ligand_rows": lig_rows}
