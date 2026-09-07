"""Ligand preparation (sections 4.2 and 5.2).

Each molecule of a library SDF is handled independently: sanitize -> analyse
-> protonate (documented simple rules) -> explicit 3D hydrogens -> a single
``mk_prepare_ligand.py`` call per molecule. This keeps the ligand-id to PDBQT
mapping deterministic and lets one broken molecule fail without killing the
library.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import List, Optional

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

from . import utils
from .utils import run_cmd, which

_CARBOXYL = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
_BASIC_N = Chem.MolFromSmarts(
    "[N;X3;!$(N-C=O);!$(N-c);!$(N=*);!$(N~[$([#6]=[#6]),$([#6]#[#6]),$([#7]=[#6]),$([#8]=[#6])]);!$(N-[OX2H1])]"
)
_BASIC_N_RING = Chem.MolFromSmarts("[N;X3;R;!$(N-C=O);!$(N-c);!$(N=*)]")


@dataclass
class LigandEntry:
    source: str
    index: int
    name: str
    formula: str = ""
    mw: float = 0.0
    charge: int = 0
    rotatable: int = 0
    heavy_atoms: int = 0
    fragments: int = 1
    has_3d: bool = False
    sanitized: bool = True
    warnings: List[str] = field(default_factory=list)
    pdbqt: str = ""
    prep_sdf: str = ""
    pdbqt_ok: bool = False
    pdbqt_note: str = ""
    ph_notes: str = ""
    smiles: str = ""          # source SMILES when the molecule came from SMILES/CSV
    mol_source: str = "sdf"   # "sdf" | "smiles"

    def to_row(self) -> dict:
        return {
            "ligand_id": self.name, "source_file": os.path.basename(self.source),
            "index_in_sdf": self.index, "formula": self.formula, "MW": round(self.mw, 2),
            "formal_charge": self.charge, "rotatable_bonds": self.rotatable,
            "heavy_atoms": self.heavy_atoms, "fragments": self.fragments,
            "has_3D": self.has_3d, "sanitized": self.sanitized,
            "ph_notes": self.ph_notes, "smiles": self.smiles,
            "pdbqt": os.path.basename(self.pdbqt) if self.pdbqt else "",
        }


def _mol_name(mol: Chem.Mol, fallback: str) -> str:
    n = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else ""
    return utils.safe_name(n or fallback)


def iter_molecules(sdf_path: str):
    suppl = Chem.SDMolSupplier(sdf_path, sanitize=True, removeHs=False)
    for i, mol in enumerate(suppl, start=1):
        yield i, mol


def analyze_molecule(mol: Chem.Mol, name: str, source: str, index: int) -> LigandEntry:
    entry = LigandEntry(source=source, index=index, name=name)
    try:
        frags = list(Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=True))
        entry.fragments = len(frags)
        if entry.fragments > 1:
            entry.warnings.append("multiple covalent fragments (not a single ligand)")
        entry.mw = Descriptors.MolWt(mol)
        entry.charge = Chem.GetFormalCharge(mol)
        entry.rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
        entry.heavy_atoms = mol.GetNumHeavyAtoms()
        entry.formula = rdMolDescriptors.CalcMolFormula(mol)
        conf = mol.GetConformer() if mol.GetNumConformers() else None
        entry.has_3d = conf is not None
        if entry.has_3d:
            ok = n = 0
            for b in mol.GetBonds():
                i1 = b.GetBeginAtomIdx(); i2 = b.GetEndAtomIdx()
                if mol.GetAtomWithIdx(i1).GetSymbol() == "H" or mol.GetAtomWithIdx(i2).GetSymbol() == "H":
                    continue
                p1 = conf.GetAtomPosition(i1); p2 = conf.GetAtomPosition(i2)
                d = p1.Distance(p2)
                n += 1
                if 0.8 < d < 2.2:
                    ok += 1
            if n and ok / n < 0.5:
                entry.has_3d = False
                entry.warnings.append("bond lengths inconsistent with a real 3D geometry")
    except Exception as e:  # pragma: no cover - defensive
        entry.sanitized = False
        entry.warnings.append(f"RDKit analysis failed: {e}")
    return entry


# Simple-rule protonation cut-offs (v1.1). These are deliberately crude and
# documented - they are NOT a pKa predictor:
#   carboxyl acids  (R-COOH, pKa ~4.5)  -> deprotonated when pH  >= 5.0
#   basic N (aliphatic amine, conj-acid pKa ~9-10) -> protonated when pH <= 8.5
# Within the physiological range used for docking this reproduces the dominant
# protomer; at pH values near a group's pKa the true population is a mixture and
# this single-protomer approximation is disclosed in every report.
_CARBOXYL_DEPROTONATE_PH = 5.0
_BASIC_N_PROTONATE_PH = 8.5


def apply_ph_rules(mol: Chem.Mol, ph: float) -> Chem.Mol:
    """Documented simple pH rules (see report; no pKa predictor used)."""
    mol = Chem.Mol(mol)
    Chem.SanitizeMol(mol)
    notes = []
    if ph >= _CARBOXYL_DEPROTONATE_PH:
        for m in mol.GetSubstructMatches(_CARBOXYL):
            ox = mol.GetAtomWithIdx(m[1])
            h = [a for a in ox.GetNeighbors() if a.GetSymbol() == "H"]
            for a in h:
                mol.RemoveAtom(a.GetIdx())
            if h:
                notes.append(f"carboxyl deprotonated (COO-) [pH {ph} >= {_CARBOXYL_DEPROTONATE_PH:g}]")
    mol.UpdatePropertyCache()
    seen = set()
    if ph <= _BASIC_N_PROTONATE_PH:
        for sm in (_BASIC_N, _BASIC_N_RING):
            for m in mol.GetSubstructMatches(sm):
                if m[0] in seen:
                    continue
                seen.add(m[0])
                a = mol.GetAtomWithIdx(m[0])
                a.SetNoImplicit(False)
                a.SetNumExplicitHs(a.GetTotalNumHs() + 1)
                a.SetFormalCharge(1)
            if seen and sm is _BASIC_N:
                # Only try ring amines if no aliphatic basic N was handled.
                break
    if seen:
        notes.append(f"basic amine protonated (+N) [pH {ph} <= {_BASIC_N_PROTONATE_PH:g}]")
    Chem.SanitizeMol(mol)
    mol.SetProp("_DockStudio_ph", f"{ph}")
    mol.SetProp("_DockStudio_ph_notes", "; ".join(notes) or "no simple-rule protonation applied")
    return mol


def _single_meeko(script: str, sdf: str, out_pdbqt: str) -> str:
    cmd = [script, "-i", sdf, "-o", out_pdbqt, "--charge_model", "gasteiger"]
    rc, out, err = run_cmd(cmd)
    if rc != 0:
        raise RuntimeError(err[-2000:] or out[-2000:])
    return out_pdbqt


def _find_meeko_ligand_script() -> str:
    script = which("mk_prepare_ligand.py", "mk_prepare_ligand")
    if script:
        return script
    try:
        import meeko
        for cand in (os.path.join(os.path.dirname(meeko.__file__), "cli", "mk_prepare_ligand.py"),
                     os.path.join(os.path.dirname(meeko.__file__), "mk_prepare_ligand.py")):
            if os.path.exists(cand):
                return cand
    except Exception:
        pass
    raise RuntimeError("mk_prepare_ligand.py not found (meeko not installed correctly)")


def prepare_ligand_pdbqt_single(sdf_path: str, out_pdbqt: str) -> str:
    """Prepare a single-molecule SDF into a PDBQT via meeko."""
    script = _find_meeko_ligand_script()
    os.makedirs(os.path.dirname(out_pdbqt) or ".", exist_ok=True)
    _single_meeko(script, sdf_path, out_pdbqt)
    return out_pdbqt


def prepare_cocrystal_pdbqt(sdf_path: str, out_pdbqt: str, ph: float = 7.4) -> str:
    """Prepare a cocrystal/reference SDF that may lack hydrogens.

    Adds explicit H (heavy-atom coordinates preserved), applies the simple pH
    rules, then runs meeko. Used for self-docking inputs where the crystal
    coordinates are the RMSD reference.
    """
    script = _find_meeko_ligand_script()
    mol = Chem.MolFromMolFile(sdf_path, sanitize=True, removeHs=False)
    if mol is None:
        raise ValueError(f"cannot parse reference SDF {sdf_path}")
    if mol.GetNumConformers() == 0:
        raise ValueError(f"reference SDF {sdf_path} has no coordinates")
    mol = Chem.AddHs(mol, addCoords=True)
    mol = apply_ph_rules(mol, ph)
    tmp_sdf = os.path.join(os.path.dirname(out_pdbqt) or ".", os.path.basename(out_pdbqt) + ".prep.sdf")
    Chem.MolToMolFile(mol, tmp_sdf)
    os.makedirs(os.path.dirname(out_pdbqt) or ".", exist_ok=True)
    _single_meeko(script, tmp_sdf, out_pdbqt)
    return out_pdbqt


def prepare_ligand_library(
    sdf_path: str,
    out_dir: str,
    ph: float = 7.4,
    library_name: str = "",
) -> List[LigandEntry]:
    """Prepare all molecules in an SDF to individual PDBQTs.

    Returns one LigandEntry per SDF record in the same order; entries for
    unparseable / multi-fragment molecules have ``pdbqt_ok=False``.
    """
    os.makedirs(out_dir, exist_ok=True)
    script = _find_meeko_ligand_script()
    entries: List[LigandEntry] = []
    base = library_name or os.path.splitext(os.path.basename(sdf_path))[0]
    for i, mol in iter_molecules(sdf_path):
        name = f"{utils.safe_name(base)}_{i}"
        if mol is None:
            entries.append(LigandEntry(source=sdf_path, index=i, name=name,
                                       sanitized=False, warnings=["unparseable SDF record"]))
            continue
        entry = analyze_molecule(mol, name, sdf_path, i)
        if not entry.sanitized:
            entries.append(entry)
            continue
        if entry.fragments > 1:
            entry.warnings.append("multi-fragment molecule not exported for docking")
            entries.append(entry)
            continue
        try:
            m3 = Chem.AddHs(Chem.Mol(mol), addCoords=True)
            if not entry.has_3d:
                m3 = Chem.AddHs(mol)
                AllChem.EmbedMolecule(m3, randomSeed=0xC0FFEE)
                entry.warnings.append("no reliable 3D coordinates - embedded a conformer")
            m4 = apply_ph_rules(m3, ph)
            m4.SetProp("_Name", name)
            if not entry.has_3d:
                try:
                    AllChem.MMFFOptimizeMolecule(m4, maxIters=200)
                except Exception:
                    pass
            prep_sdf = os.path.join(out_dir, f"{name}.prep.sdf")
            Chem.MolToMolFile(m4, prep_sdf)
            entry.prep_sdf = prep_sdf
            entry.ph_notes = m4.GetProp("_DockStudio_ph_notes")
            out_pdbqt = os.path.join(out_dir, f"{name}.pdbqt")
            _single_meeko(script, prep_sdf, out_pdbqt)
            entry.pdbqt = out_pdbqt
            entry.pdbqt_ok = True
        except Exception as e:
            entry.warnings.append(f"preparation failed: {e}")
            entry.pdbqt_note = str(e)
        entries.append(entry)
    return entries


def torsdof_of(pdbqt_path: str) -> int:
    try:
        for line in open(pdbqt_path, errors="replace"):
            if "TORSDOF" in line:
                m = re.search(r"TORSDOF\s+(\d+)", line)
                if m:
                    return int(m.group(1))
    except Exception:
        pass
    return -1


# ---------------------------------------------------------------------------
# v2.0 - SMILES / CSV ligand input
# ---------------------------------------------------------------------------
def is_smiles_source(path: str) -> bool:
    ext = os.path.splitext(str(path))[1].lower()
    return ext in (".smi", ".smiles", ".csv", ".txt")


def read_smiles_records(path: str) -> List[dict]:
    """Parse a SMILES/CSV ligand source into records.

    Each record: ``{"index": int, "smiles": str, "name": str}``.
    Supported:
      * ``.smi``/``.txt`` - one SMILES per line; optional name after whitespace;
        ``#`` comments and blank lines are ignored.
      * ``.csv`` - header must contain a ``SMILES`` column (case-insensitive);
        an optional ``ID``/``NAME``/``TITLE`` column is used as the name hint.
    Raises ``ValueError`` on files with no usable SMILES.
    """
    path = str(path)
    ext = os.path.splitext(path)[1].lower()
    records: List[dict] = []
    if ext == ".csv":
        import csv as _csv
        with open(path, encoding="utf-8-sig", newline="") as fh:
            rows = list(_csv.reader(fh))
        if not rows:
            raise ValueError(f"SMILES CSV 为空: {path}")
        header = [h.strip().lower() for h in rows[0]]
        si = next((i for i, h in enumerate(header)
                   if h in ("smiles", "smi", "canonical_smiles", "isosmiles")), None)
        if si is None:
            raise ValueError(f"CSV 需包含 SMILES 列(表头之一:SMILES/SMI/…): {path}")
        ni = next((i for i, h in enumerate(header)
                   if h in ("id", "name", "title", "ligand_id", "compound", "mol_name")), None)
        idx = 0
        for row in rows[1:]:
            if len(row) <= si:
                continue
            smi = row[si].strip()
            if not smi or smi.startswith("#"):
                continue
            idx += 1
            nm = ""
            if ni is not None and len(row) > ni and row[ni].strip():
                nm = row[ni].strip()
            records.append({"index": idx, "smiles": smi, "name": nm})
    else:
        with open(path, encoding="utf-8-sig") as fh:
            idx = 0
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                toks = line.split()
                smi = toks[0]
                nm = " ".join(toks[1:]).strip()
                idx += 1
                records.append({"index": idx, "smiles": smi, "name": nm})
    if not records:
        raise ValueError(f"未在文件中解析到任何 SMILES: {path}")
    return records


def _entry_name(base: str, rec: dict) -> str:
    hint = utils.safe_name(rec.get("name") or "")
    if hint:
        return hint
    return f"{base}_{rec['index']}"


def prepare_smiles_library(
    source_path: str,
    out_dir: str,
    ph: float = 7.4,
    library_name: str = "",
) -> List[LigandEntry]:
    """Build individual PDBQTs from a SMILES/CSV ligand source (v2.0).

    Each SMILES is converted to a 3D conformer with RDKit ETKDG (then MMFF
    minimised), protonated with the documented simple pH rules, written to a
    ``prep.sdf`` and passed to meeko exactly like an SDF record. Failures are
    reported per molecule (``pdbqt_ok=False`` with a warning), never fatal for
    the whole library.
    """
    os.makedirs(out_dir, exist_ok=True)
    script = _find_meeko_ligand_script()
    base = utils.safe_name(library_name or os.path.splitext(os.path.basename(source_path))[0])
    entries: List[LigandEntry] = []
    for rec in read_smiles_records(source_path):
        name = _entry_name(base, rec)
        entry = LigandEntry(source=source_path, index=rec["index"], name=name,
                            mol_source="smiles", smiles=rec["smiles"])
        mol = None
        try:
            mol = Chem.MolFromSmiles(rec["smiles"])
        except Exception:
            mol = None
        if mol is None:
            entry.sanitized = False
            entry.warnings.append("unparseable SMILES")
            entries.append(entry)
            continue
        entry = analyze_molecule(mol, name, source_path, rec["index"])
        entry.mol_source = "smiles"
        entry.smiles = rec["smiles"]
        if not entry.sanitized:
            entries.append(entry)
            continue
        if entry.fragments > 1:
            entry.warnings.append("multi-fragment molecule not exported for docking")
            entries.append(entry)
            continue
        try:
            m3 = Chem.AddHs(mol)
            ok = AllChem.EmbedMolecule(m3, randomSeed=0xC0FFEE)
            if ok != 0:
                # retry once with a different seed after clearing the conformer
                m3.RemoveAllConformers()
                ok = AllChem.EmbedMolecule(m3, randomSeed=0xC0FFEE + 1)
            if ok != 0:
                raise ValueError("3D embedding failed for SMILES")
            try:
                AllChem.MMFFOptimizeMolecule(m3, maxIters=200)
            except Exception:
                pass
            entry.has_3d = True
            entry.warnings.append("no input 3D coordinates - embedded from SMILES (ETKDG+MMFF)")
            m4 = apply_ph_rules(m3, ph)
            m4.SetProp("_Name", name)
            prep_sdf = os.path.join(out_dir, f"{name}.prep.sdf")
            Chem.MolToMolFile(m4, prep_sdf)
            entry.prep_sdf = prep_sdf
            entry.ph_notes = m4.GetProp("_DockStudio_ph_notes")
            out_pdbqt = os.path.join(out_dir, f"{name}.pdbqt")
            _single_meeko(script, prep_sdf, out_pdbqt)
            entry.pdbqt = out_pdbqt
            entry.pdbqt_ok = True
        except Exception as e:
            entry.warnings.append(f"preparation failed: {e}")
            entry.pdbqt_note = str(e)
        entries.append(entry)
    return entries


def prepare_ligand_file(
    path: str,
    out_dir: str,
    ph: float = 7.4,
    library_name: str = "",
) -> List[LigandEntry]:
    """Prepare every molecule of an SDF/MOL *or* SMILES/CSV ligand source.

    Dispatches on file extension; returns ``List[LigandEntry]`` with the same
    semantics as :func:`prepare_ligand_library`.
    """
    if is_smiles_source(path):
        return prepare_smiles_library(path, out_dir, ph=ph, library_name=library_name)
    return prepare_ligand_library(path, out_dir, ph=ph, library_name=library_name)
