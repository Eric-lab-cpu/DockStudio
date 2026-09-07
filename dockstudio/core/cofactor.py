"""Co-crystal ligand reference structures (section 5.3).

Connectivity is taken from the RCSB Chemical Component Dictionary (CCD)
``.cif``; atomic coordinates are taken from the *local* PDB (HETATM), as the
protocol demands. CCD model coordinates are used only for a sanity check.
"""

from __future__ import annotations

import dataclasses
import os
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops

from . import structure as st
from . import utils


@dataclass
class CCDAtom:
    atom_id: str
    element: str
    charge: int
    x: float
    y: float
    z: float
    aromatic: bool = False


@dataclass
class CCDBond:
    a1: str
    a2: str
    order: int
    aromatic: bool = False


@dataclass
class CCD:
    comp_id: str
    atoms: List[CCDAtom] = field(default_factory=list)
    bonds: List[CCDBond] = field(default_factory=list)

    def centroid_model(self) -> List[float]:
        if not self.atoms:
            return [0.0, 0.0, 0.0]
        n = len(self.atoms)
        return [sum(a.x for a in self.atoms) / n,
                sum(a.y for a in self.atoms) / n,
                sum(a.z for a in self.atoms) / n]

    def element(self, atom_id: str) -> Optional[str]:
        for a in self.atoms:
            if a.atom_id == atom_id:
                return a.element
        return None


def parse_ccd(path: str) -> CCD:
    """Parse the atom/bond loops of a CCD mmCIF file (simple whitespace parser)."""
    lines = open(path, encoding="utf-8").read().splitlines()
    cc = CCD(comp_id=os.path.splitext(os.path.basename(path))[0])
    mode = None
    atom_cols = {}
    bond_cols = {}
    for line in lines:
        ls = line.strip()
        if ls.startswith("#"):
            mode = None
            continue
        if ls.startswith("loop_"):
            mode = "header"
            continue
        if ls.startswith("_"):
            if "_chem_comp_atom." in ls:
                key = ls.split(".")[1]
                atom_cols[key] = len(atom_cols)
            elif "_chem_comp_bond." in ls:
                key = ls.split(".")[1]
                bond_cols[key] = len(bond_cols)
            elif "_chem_comp.name" in ls and not cc.comp_id:
                pass
            mode = "header"
            continue
        # data row
        toks = ls.split()
        if atom_cols and not bond_cols:
            if len(toks) >= 15 and mode == "header" and toks[0] == cc.comp_id.upper():
                pass
        if not toks:
            continue
        # decide section by which column map is active from header context
        # (heuristic: atom rows begin with comp_id and have >=15 fields)
        if len(toks) >= 15 and toks[0].upper() == cc.comp_id.upper():
            atom_id = toks[1]
            elem = toks[3]
            charge = int(float(toks[4]))
            try:
                x, y, z = float(toks[12]), float(toks[13]), float(toks[14])
            except (ValueError, IndexError):
                x = y = z = 0.0
            arom = (toks[6] if len(toks) > 6 else "?") == "Y"
            cc.atoms.append(CCDAtom(atom_id, elem, charge, x, y, z, arom))
        elif len(toks) >= 7 and toks[0].upper() == cc.comp_id.upper() and atom_cols:
            # bond rows: comp_id atom1 atom2 value_order aromatic stereo ordinal
            try:
                order_s = toks[3].upper()
                # CCD value_order codes are upper-case: SING/DOUB/TRIP/AROM/DELO
                order = {"SING": 1, "SINGLE": 1, "DOUB": 2, "DOUBLE": 2,
                         "TRIP": 3, "TRIPLE": 3, "AROM": 1, "AROMATIC": 1,
                         "DELO": 1}.get(order_s, 1)
                arom = toks[4] == "Y" or order_s in ("AROM", "AROMATIC")
                cc.bonds.append(CCDBond(toks[1], toks[2], order, arom))
            except Exception:
                pass
    # Fallback: precise two-phase scan if heuristic missed rows
    if len(cc.atoms) < 2 or len(cc.bonds) < 2:
        return _parse_ccd_phase2(lines, cc.comp_id)
    return cc


def _parse_ccd_phase2(lines: List[str], comp_id: str) -> CCD:
    """Deterministic phase-2 parser driven by explicit column headers."""
    cc = CCD(comp_id=comp_id)
    sections = []  # list of dict(name, cols, rows)
    cur = None
    for line in lines:
        ls = line.strip()
        if ls.startswith("_chem_comp_atom."):
            if cur and cur["kind"] == "atom":
                cur = None
            cur = {"kind": "atom", "cols": [], "rows": []}
            sections.append(cur)
            cur["cols"].append(ls.split(".")[1])
        elif ls.startswith("_chem_comp_bond."):
            if cur and cur["kind"] == "bond":
                cur = None
            cur = {"kind": "bond", "cols": [], "rows": []}
            sections.append(cur)
            cur["cols"].append(ls.split(".")[1])
        elif ls.startswith("_") or ls.startswith("loop_"):
            continue
        elif cur is not None and ls and not ls.startswith("#"):
            cur["rows"].append(ls.split())
        elif ls.startswith("#"):
            cur = None
    for sec in sections:
        if sec["kind"] == "atom":
            want = ["comp_id", "atom_id", "type_symbol", "charge", "model_Cartn_x",
                    "model_Cartn_y", "model_Cartn_z", "pdbx_aromatic_flag"]
            idx = {c: sec["cols"].index(c) for c in want if c in sec["cols"]}
            for row in sec["rows"]:
                if not row:
                    continue
                try:
                    cc.atoms.append(CCDAtom(
                        row[idx["atom_id"]],
                        row[idx["type_symbol"]],
                        int(float(row[idx["charge"]])),
                        float(row[idx["model_Cartn_x"]]),
                        float(row[idx["model_Cartn_y"]]),
                        float(row[idx["model_Cartn_z"]]),
                        row[idx["pdbx_aromatic_flag"]] == "Y",
                    ))
                except (KeyError, ValueError, IndexError):
                    continue
        else:
            want = ["comp_id", "atom_id_1", "atom_id_2", "value_order", "pdbx_aromatic_flag"]
            idx = {c: sec["cols"].index(c) for c in want if c in sec["cols"]}
            for row in sec["rows"]:
                if not row:
                    continue
                try:
                    os_ = row[idx["value_order"]].upper()
                    order = {"SING": 1, "SINGLE": 1, "DOUB": 2, "DOUBLE": 2,
                             "TRIP": 3, "TRIPLE": 3}.get(os_, 1)
                    arom = row[idx["pdbx_aromatic_flag"]] == "Y" or os_ in ("AROM", "AROMATIC")
                    cc.bonds.append(CCDBond(row[idx["atom_id_1"]], row[idx["atom_id_2"]], order, arom))
                except (KeyError, ValueError, IndexError):
                    continue
    return cc


# ---------------------------------------------------------------------------
# SDF construction
# ---------------------------------------------------------------------------
def ccd_to_sdf_local_coords(
    pdb_path: str,
    resname: str,
    ccd: CCD,
    out_sdf: str,
    chain: str = "",
    resseq: Optional[int] = None,
) -> dict:
    """Build the crystal-reference SDF: CCD connectivity, local PDB coordinates.

    Returns a stats dict (n_atoms, centroid checks, sanitize status).
    """
    local = [r for r in st.parse_structure(pdb_path)
             if r.record == "HETATM" and r.resname == resname
             and (not chain or r.chain == chain)
             and (resseq is None or r.resseq == resseq)]
    if not local:
        # fall back to any chain instance
        local = [r for r in st.parse_structure(pdb_path)
                 if r.record == "HETATM" and r.resname == resname]
    if not local:
        raise ValueError(f"residue {resname} not found in {pdb_path}")
    local_heavy = [r for r in local if r.element != "H"]
    local_centroid = st.heavy_centroid(local_heavy)
    # NOTE: CCD *model* coordinates and the local PDB coordinates are in two
    # unrelated frames, so a centroid offset is NOT a quality metric.  It is kept
    # here only as an informational field.  The meaningful quality checks are the
    # heavy-atom coverage of the CCD connectivity by the local PDB atoms and
    # successful sanitisation of the assembled molecule.
    ccd_centroid = ccd.centroid_model()
    dev = utils.distance(local_centroid, ccd_centroid)
    ccd_heavy_atoms = [a for a in ccd.atoms if a.element not in ("H", "D")]
    local_heavy_names = {a.name for a in local_heavy}
    ccd_heavy_names = {a.atom_id for a in ccd_heavy_atoms}
    n_common = len(ccd_heavy_names & local_heavy_names)
    coverage = n_common / max(len(ccd_heavy_names), 1)

    # map CCD atom_id -> local atom
    coord = {}
    for a in local_heavy:
        coord.setdefault(a.name, a.coord())
    # Build RWMol from CCD (only atoms present locally)
    rw = Chem.RWMol()
    idx_of = {}
    for a in ccd.atoms:
        if a.atom_id not in coord:
            continue
        at = Chem.Atom(a.element)
        if a.element == "H":
            at.SetNoImplicit(True)
            at.SetNumExplicitHs(0)
        if a.charge:
            at.SetFormalCharge(a.charge)
        if a.element == "D":
            continue
        idx_of[a.atom_id] = rw.AddAtom(at)
    for b in ccd.bonds:
        if b.a1 in idx_of and b.a2 in idx_of:
            rw.AddBond(idx_of[b.a1], idx_of[b.a2],
                       Chem.BondType.AROMATIC if b.aromatic else
                       Chem.BondType.SINGLE if b.order == 1 else
                       Chem.BondType.DOUBLE if b.order == 2 else Chem.BondType.TRIPLE)

    mol = rw.GetMol()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for aid, cid in idx_of.items():
        x, y, z = coord[aid]
        conf.SetAtomPosition(cid, Chem.rdGeometry.Point3D(x, y, z))
    mol.AddConformer(conf)

    stats = {"n_atoms": mol.GetNumAtoms(), "n_heavy": mol.GetNumHeavyAtoms(),
             "n_ccd_heavy": len(ccd_heavy_atoms),
             "n_local_heavy": len(local_heavy),
             "heavy_coverage_local_vs_ccd": round(coverage, 3),
             "ccd_model_centroid_offset_A": round(dev, 3),
             "sanitized": False,
             "fallback": ""}
    try:
        Chem.SanitizeMol(mol)
        stats["sanitized"] = True
    except Exception as e:
        # keep unsanitized connectivity but still write - meeko/RDKit may cope
        stats["fallback"] = f"CCD sanitize failed ({e}); connectivity kept as-is"
    try:
        Chem.MolToMolFile(mol, out_sdf)
    except Exception as e:
        stats["fallback"] += f" ; write failed ({e})"
        # fallback: PDB-distance based mol
        mol2 = _mol_from_pdb_heuristic(local, ccd)
        if mol2 is not None:
            Chem.MolToMolFile(mol2, out_sdf)
            stats["fallback"] += " ; used PDB distance-based connectivity"
    return stats


def _mol_from_pdb_heuristic(recs, ccd: CCD) -> Optional[Chem.Mol]:
    """Ultimate fallback: RDKit distance-based connectivity from PDB lines."""
    lines = []
    for i, a in enumerate(recs, 1):
        lines.append(st.format_pdb_atom(i, a))
    block = "".join(lines)
    try:
        mol = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=False)
        if mol is None:
            return None
        Chem.SanitizeMol(mol)
        return mol
    except Exception:
        return None


def write_cocrystal_ligand(
    pdb_path: str,
    resname: str,
    ccd_path: str,
    out_sdf: str,
    chain: str = "",
    resseq: Optional[int] = None,
    name: str = "",
) -> dict:
    ccd = parse_ccd(ccd_path)
    os.makedirs(os.path.dirname(out_sdf) or ".", exist_ok=True)
    stats = ccd_to_sdf_local_coords(pdb_path, resname, ccd, out_sdf, chain=chain, resseq=resseq)
    stats["resname"] = resname
    stats["ccd_file"] = os.path.basename(ccd_path)
    stats["chain"] = chain
    stats["resseq"] = resseq
    stats["name"] = name or resname
    return stats
