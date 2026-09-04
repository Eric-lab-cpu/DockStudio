"""PDB parsing, inventory and receptor-cleaning helpers.

gemmi is used only for *reading*; writing uses an explicit PDB line
formatter that was validated byte-for-byte against gemmi output so that
meeko / PyMOL / PLIP accept the files without reformatting surprises.
"""

from __future__ import annotations

import dataclasses
import math
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from . import models
from .utils import md5_file, now_str


# ---------------------------------------------------------------------------
# PDB record parsing
# ---------------------------------------------------------------------------
@dataclass
class AtomRec:
    record: str
    serial: int
    name: str
    altloc: str
    resname: str
    chain: str
    resseq: int
    icode: str
    x: float
    y: float
    z: float
    occ: float
    b: float
    element: str

    def coord(self) -> List[float]:
        return [self.x, self.y, self.z]

    def reskey(self) -> Tuple[str, str, int, str]:
        return (self.chain, self.resname, self.resseq, self.icode)


def parse_pdb(path: str, keep_chains: Optional[List[str]] = None) -> List[AtomRec]:
    recs: List[AtomRec] = []
    for line in open(path, "r", errors="replace"):
        if len(line) < 54:
            continue
        rec = line[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue
        chain = line[21] if len(line) > 21 else " "
        if keep_chains and chain not in keep_chains:
            continue
        name = line[12:16].strip()
        alt = line[16] if len(line) > 16 else " "
        resn = line[17:20].strip()
        try:
            seq = int(line[22:26])
        except ValueError:
            seq = 0
        icode = line[26] if len(line) > 26 else " "
        try:
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        except ValueError:
            continue
        try:
            occ = float(line[54:60]) if line[54:60].strip() else 1.0
        except ValueError:
            occ = 1.0
        try:
            b = float(line[60:66]) if line[60:66].strip() else 0.0
        except ValueError:
            b = 0.0
        el = line[76:78].strip()
        if not el:
            # infer element from atom name
            el = name[0] if name else "C"
        recs.append(AtomRec(rec, len(recs) + 1, name, alt, resn, chain, seq, icode,
                            x, y, z, occ, b, el))
    return recs


def _name_field(name: str) -> str:
    return (" " + name.ljust(3))[:4]


def format_pdb_atom(serial: int, a: AtomRec, altloc: str = " ") -> str:
    el = a.element
    return (
        f"{a.record:<6s}{serial:5d} {_name_field(a.name)}{altloc}{a.resname:>3s} "
        f"{a.chain}{a.resseq:4d}{a.icode}   {a.x:8.3f}{a.y:8.3f}{a.z:8.3f}"
        f"{a.occ:6.2f}{a.b:6.2f}          {el:>2s}  "
    )


def write_pdb(recs: List[AtomRec], out_path: str, crystal_line: Optional[str] = None,
              remark: Optional[List[str]] = None) -> None:
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        if crystal_line:
            fh.write(crystal_line if crystal_line.endswith("\n") else crystal_line + "\n")
        for r in remark or []:
            fh.write(f"REMARK  {r}\n")
        for i, a in enumerate(recs, start=1):
            fh.write(format_pdb_atom(i, a) + "\n")
        fh.write("END\n")


def first_crystal_line(path: str) -> Optional[str]:
    for line in open(path, "r", errors="replace"):
        if line.startswith("CRYST1"):
            return line.rstrip("\n")
        if line.startswith(("ATOM", "HETATM")):
            break
    return None


# ---------------------------------------------------------------------------
# geometry helpers
# ---------------------------------------------------------------------------
def heavy_centroid(recs: List[AtomRec]) -> List[float]:
    h = [r for r in recs if r.element != "H"]
    if not h:
        h = recs
    n = len(h)
    return [sum(r.x for r in h) / n, sum(r.y for r in h) / n, sum(r.z for r in h) / n]


def coord_extent(recs: List[AtomRec], pad: float = 0.0) -> Tuple[List[float], List[float]]:
    """Return (min, max) per axis over atoms (heavy only if present)."""
    h = [r for r in recs if r.element != "H"] or recs
    xs = [r.x for r in h]; ys = [r.y for r in h]; zs = [r.z for r in h]
    lo = [min(xs) - pad, min(ys) - pad, min(zs) - pad]
    hi = [max(xs) + pad, max(ys) + pad, max(zs) + pad]
    return lo, hi


def dist_to_box(p: List[float], center: List[float], size: List[float]) -> float:
    """Minimal distance from point p to the box defined by center/size."""
    d = 0.0
    for c, s, v in zip(center, size, p):
        lo = c - s / 2.0
        hi = c + s / 2.0
        if v < lo:
            d += (lo - v) ** 2
        elif v > hi:
            d += (v - hi) ** 2
    return math.sqrt(d)


# ---------------------------------------------------------------------------
# inventory (section 4)
# ---------------------------------------------------------------------------
@dataclass
class ChainInventory:
    chain: str
    polymer_res: int = 0
    res_range: Tuple[int, int] = (0, 0)
    het_ligands: List[dict] = field(default_factory=list)
    waters: int = 0
    altloc_residues: List[str] = field(default_factory=list)
    modified_residues: List[dict] = field(default_factory=list)
    nucleic_res: int = 0


@dataclass
class StructureInventory:
    path: str
    name: str
    md5: str = ""
    chains: Dict[str, ChainInventory] = field(default_factory=dict)
    oligomer_hint: str = ""
    parsed_atoms: int = 0

    def to_dict(self) -> dict:
        d = dataclasses.asdict(self)
        d.pop("path", None)
        return d


def inventory_pdb(path: str, name: str = "", fetch_online_md5: Optional[str] = None) -> StructureInventory:
    """Line-based inventory of a PDB file (section 4.1)."""
    recs = parse_pdb(path)
    inv = StructureInventory(path=path, name=name or os.path.splitext(os.path.basename(path))[0])
    if os.path.exists(path):
        try:
            inv.md5 = md5_file(path)
        except Exception:
            inv.md5 = ""
    inv.parsed_atoms = len(recs)

    by_chain: Dict[str, List[AtomRec]] = {}
    for r in recs:
        by_chain.setdefault(r.chain, []).append(r)

    for chain, crecs in sorted(by_chain.items()):
        ci = ChainInventory(chain=chain)
        resseqs: List[int] = []
        ligand_map: Dict[tuple, List[AtomRec]] = {}
        mod_map: Dict[tuple, List[AtomRec]] = {}
        alt_res: set = set()
        for r in crecs:
            if r.record == "HETATM":
                resn = r.resname
                if resn in ("HOH", "WAT"):
                    ci.waters += 1
                elif resn in models.MODIFIED_TO_STANDARD:
                    mod_map.setdefault(r.reskey(), []).append(r)
                else:
                    ligand_map.setdefault(r.reskey(), []).append(r)
            else:
                if r.resname in models.NUCLEIC_RES:
                    ci.nucleic_res += 1
                    resseqs.append(r.resseq)
                else:
                    resseqs.append(r.resseq)
                if r.altloc not in ("", " "):
                    alt_res.add(r.resname + str(r.resseq))
        seqs = sorted(set(resseqs))
        ci.res_range = (seqs[0], seqs[-1]) if seqs else (0, 0)
        ci.polymer_res = len(seqs)
        ci.altloc_residues = sorted(alt_res)
        for key, rr in ligand_map.items():
            heavy = [a for a in rr if a.element != "H"]
            if len(heavy) >= 4:
                ci.het_ligands.append({
                    "resname": key[1], "chain": key[0], "resseq": key[2],
                    "icode": key[3], "n_heavy": len(heavy),
                })
        for key, rr in mod_map.items():
            ci.modified_residues.append({
                "resname": key[1], "chain": key[0], "resseq": key[2],
                "icode": key[3], "std": models.MODIFIED_TO_STANDARD.get(key[1], "?"),
                "n_atoms": len(rr),
            })
        inv.chains[chain] = ci

    n_poly_chains = sum(1 for c in inv.chains.values() if c.polymer_res > 0)
    inv.oligomer_hint = ("multimer" if n_poly_chains > 1 else
                         "monomer" if n_poly_chains == 1 else "none")
    return inv


# ---------------------------------------------------------------------------
# cleaning (section 5.1)
# ---------------------------------------------------------------------------
@dataclass
class CleanAction:
    dropped: Dict[str, int] = field(default_factory=dict)
    converted: List[dict] = field(default_factory=list)
    kept: int = 0
    warnings: List[str] = field(default_factory=list)


def group_residues(recs: List[AtomRec]) -> Dict[tuple, List[AtomRec]]:
    out: Dict[tuple, List[AtomRec]] = {}
    for r in recs:
        key = (r.chain, r.resname, r.resseq, r.icode)
        out.setdefault(key, []).append(r)
    return out


def clean_receptor(
    pdb_path: str,
    out_path: str,
    keep_chains: Optional[List[str]] = None,
    convert_modified: bool = True,
    drop_hetligands: bool = True,
    drop_ions_solvents: bool = True,
    drop_nucleic: bool = True,
    box_center: Optional[List[float]] = None,
    box_size: Optional[List[float]] = None,
    name: str = "",
) -> CleanAction:
    """Clean a receptor PDB for docking.

    Rules (protocol section 5.1): keep only chosen chains, standard amino
    acids; drop waters / ions / co-crystal ligands / nucleic acids; keep
    altloc 'A' (blanking the altloc column); convert modified residues to
    their standard amino acid (recording the action). Residues whose
    deletion is scientifically questionable inside the box are *not*
    silently removed; they are returned in ``action.warnings`` for the
    caller to decide.
    """
    recs = parse_pdb(pdb_path, keep_chains=keep_chains)
    crystal = first_crystal_line(pdb_path)
    action = CleanAction()

    # decide which residues to keep
    groups = group_residues(recs)
    kept: List[AtomRec] = []
    for key, rr in groups.items():
        chain, resn, seq, icode = key
        record = rr[0].record
        heavy = [a for a in rr if a.element != "H"]
        center = heavy_centroid(rr)

        if resn in ("HOH", "WAT") or (drop_ions_solvents and resn in models.IONS_SOLVENTS):
            action.dropped[resn] = action.dropped.get(resn, 0) + 1
            continue

        if record == "HETATM":
            if resn in models.MODIFIED_TO_STANDARD:
                if convert_modified:
                    # convert to standard residue
                    new_resn = models.MODIFIED_TO_STANDARD[resn]
                    kept_atoms = []
                    for a in rr:
                        if a.element == "H":
                            continue
                        if resn in ("TPO", "SEP", "PTR") and a.name in models.PHOSPHO_DROP_ATOMS:
                            continue
                        na = dataclasses.replace(a, resname=new_resn, record="ATOM", altloc=" ")
                        kept_atoms.append(na)
                    action.converted.append({"from": resn, "to": new_resn,
                                             "chain": chain, "resseq": seq,
                                             "n_atoms_removed": len(rr) - len(kept_atoms)})
                    if box_center and box_size and dist_to_box(center, box_center, box_size) < 6.0:
                        action.warnings.append(
                            f"modified residue {chain}:{resn}{seq} converted to {new_resn} "
                            f"but lies within 6 A of the search box (distance to box = "
                            f"{dist_to_box(center, box_center, box_size):.1f} A).")
                    kept.extend(kept_atoms)
                else:
                    action.dropped[resn] = action.dropped.get(resn, 0) + 1
                    action.warnings.append(
                        f"modified residue {chain}:{resn}{seq} dropped (convert_modified=False).")
                continue
            # true HET ligand / cofactor
            if drop_hetligands:
                action.dropped[resn] = action.dropped.get(resn, 0) + 1
                continue
            kept.extend(rr)
            continue

        # ATOM records
        if drop_nucleic and resn in models.NUCLEIC_RES:
            action.dropped[resn] = action.dropped.get(resn, 0) + 1
            continue
        # keep only A / blank altloc conformers
        for a in rr:
            if a.altloc in ("", " ", "A"):
                na = dataclasses.replace(a, altloc=" ")
                kept.append(na)

    action.kept = len(kept)
    remark = [
        f"DockStudio cleaned receptor from {os.path.basename(pdb_path)}",
        f"chains={','.join(keep_chains) if keep_chains else 'all'} "
        f"atoms_kept={action.kept} dropped={action.dropped} converted={len(action.converted)}",
        "altloc conformer A retained; waters/ions/ligands/nucleic acids removed.",
    ]
    write_pdb(kept, out_path, crystal_line=crystal, remark=remark)
    return action


def coords_of_ligand(pdb_path: str, resname: str, chain: str = "", resseq: Optional[int] = None) -> List[AtomRec]:
    """Return atom records of a specific HETATM residue (e.g. co-crystal ligand)."""
    out = []
    for r in parse_pdb(pdb_path):
        if r.record == "HETATM" and r.resname == resname:
            if chain and r.chain != chain:
                continue
            if resseq is not None and r.resseq != resseq:
                continue
            out.append(r)
    return out
