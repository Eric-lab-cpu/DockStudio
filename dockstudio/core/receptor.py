"""Receptor preparation: cleaning + PDBQT generation (sections 5.1, 3)."""

from __future__ import annotations

import dataclasses
import math
import os
import re
from typing import List, Optional

from . import models, structure as st, utils
from .utils import run_cmd, which

_BAD_RES_LINE = re.compile(r"residue_key='([A-Za-z0-9]+):(-?\d+)'")
_BAD_RES_FAIL = re.compile(r"failed for:\s*\[([^\]]*)\]")


def _find_meeko_receptor_script() -> str:
    script = which("mk_prepare_receptor.py", "mk_prepare_receptor")
    if script:
        return script
    # frozen / pip layout: scripts live inside the meeko package
    try:
        import meeko
        for cand in (os.path.join(os.path.dirname(meeko.__file__), "cli", "mk_prepare_receptor.py"),
                     os.path.join(os.path.dirname(meeko.__file__), "mk_prepare_receptor.py")):
            if os.path.exists(cand):
                return cand
    except Exception:
        pass
    raise RuntimeError("meeko not importable; receptor preparation requires meeko")


def _parse_bad_residues(log: str) -> List[tuple]:
    keys = set()
    for m in _BAD_RES_LINE.finditer(log):
        keys.add((m.group(1), int(m.group(2))))
    for m in _BAD_RES_FAIL.finditer(log):
        for tok in m.group(1).split(","):
            tok = tok.strip()
            mm = re.match(r"'?([A-Za-z0-9]+):(-?\d+)'?", tok)
            if mm:
                keys.add((mm.group(1), int(mm.group(2))))
    return sorted(keys)


def _classify_bad_residues(clean_pdb: str, bad_keys: List[tuple],
                           box_center: Optional[List[float]],
                           box_size: Optional[List[float]],
                           reference_points: Optional[List[list]] = None,
                           threshold: float = 10.0):
    """Return (far_keys, near_keys).

    A residue counts as 'near' if any heavy atom lies within ``threshold`` A of
    any reference point (cocrystal ligand heavy atoms) - the criterion required
    by the protocol. Without reference points we fall back to distance to the
    box *surface* plus the box half-extent so that large boxes do not make
    distant residues look 'near'."""
    by_key: dict = {}
    for r in st.parse_pdb(clean_pdb):
        by_key.setdefault((r.chain, r.resseq), []).append(r)
    far, near = [], []
    for (chain, seq) in bad_keys:
        recs = by_key.get((chain, seq), [])
        if not recs:
            continue
        if reference_points:
            dmin = min(min(math.dist(r.coord(), p) for p in reference_points) for r in recs)
        elif box_center is not None and box_size is not None:
            # distance to box surface + half extent
            surf = min(st.dist_to_box(r.coord(), box_center, box_size) for r in recs)
            half = max(box_size) / 2.0
            dmin = surf + half
        else:
            near.append((chain, seq))
            continue
        if dmin > threshold:
            far.append((chain, seq))
        else:
            near.append((chain, seq))
    return far, near


def _delete_residues(clean_pdb: str, keys: List[tuple], out_path: str) -> None:
    drop = {(c, s) for (c, s) in keys}
    recs = [r for r in st.parse_pdb(clean_pdb) if (r.chain, r.resseq) not in drop]
    crystal = st.first_crystal_line(clean_pdb)
    st.write_pdb(recs, out_path, crystal_line=crystal,
                 remark=[f"DockStudio deleted {len(drop)} residues with incomplete templates"])
    return


def prepare_receptor(
    pdb_path: str,
    out_dir: str,
    name: str,
    chains: Optional[List[str]] = None,
    convert_modified: bool = True,
    box_center: Optional[List[float]] = None,
    box_size: Optional[List[float]] = None,
    reference_points: Optional[List[list]] = None,
    delete_bad_res_far: bool = True,
) -> dict:
    """Full receptor preparation: clean PDB -> PDBQT via meeko.

    Returns a dict with keys: cleaned_pdb, pdbqt, action (CleanAction),
    log, atom_counts, warnings, deleted_bad_residues.
    """
    os.makedirs(out_dir, exist_ok=True)
    clean_pdb = os.path.join(out_dir, f"{name}.clean.pdb")
    clean_pdb2 = os.path.join(out_dir, f"{name}.clean.deleted.pdb")
    pdbqt = os.path.join(out_dir, f"{name}.pdbqt")

    action = st.clean_receptor(
        pdb_path, clean_pdb, keep_chains=chains, convert_modified=convert_modified,
        box_center=box_center, box_size=box_size, name=name,
    )
    script = _find_meeko_receptor_script()
    base_clean = clean_pdb
    deleted_bad = []
    log = ""
    rc, out, err = run_cmd([script, "--read_pdb", clean_pdb, "-o",
                            os.path.splitext(pdbqt)[0], "-p"])
    log = out + err
    if rc != 0:
        bad = _parse_bad_residues(log)
        far, near = _classify_bad_residues(clean_pdb, bad, box_center, box_size,
                                           reference_points=reference_points)
        if near:
            raise RuntimeError(
                "meeko receptor preparation failed for residues near the search box:\n"
                + ", ".join(f"{c}:{s}" for c, s in near)
                + "\nThese cannot be deleted silently - please inspect the structure.")
        if bad and delete_bad_res_far and far:
            deleted_bad = far
            _delete_residues(clean_pdb, far, clean_pdb2)
            base_clean = clean_pdb2
            rc, out, err = run_cmd([script, "--read_pdb", clean_pdb2, "-o",
                                    os.path.splitext(pdbqt)[0], "-p"])
            log += "\n[second attempt after deleting far incomplete residues]\n" + out + err
            if rc != 0:
                bad2 = _parse_bad_residues(out + err)
                raise RuntimeError("receptor preparation still failed after deleting far bad residues: "
                                   + ", ".join(f"{c}:{s}" for c, s in bad2))
        elif bad:
            raise RuntimeError("receptor preparation failed for residues "
                               + ", ".join(f"{c}:{s}" for c, s in bad))
    if not os.path.exists(pdbqt):
        raise RuntimeError(f"receptor PDBQT not produced for {name}; meeko log:\n{log[-2000:]}")

    atom_counts = _count_atom_types(pdbqt)
    return {
        "name": name,
        "cleaned_pdb": base_clean,
        "pdbqt": pdbqt,
        "action": action,
        "log": log,
        "atom_counts": atom_counts,
        "deleted_bad_residues": [{"chain": c, "resseq": s} for (c, s) in deleted_bad],
        "n_atoms_pdbqt": sum(atom_counts.values()),
    }


def _count_atom_types(pdbqt_path: str) -> dict:
    counts = {}
    for line in open(pdbqt_path, errors="replace"):
        if line.startswith(("ATOM", "HETATM")):
            t = line[76:79].strip()
            counts[t] = counts.get(t, 0) + 1
    return counts
