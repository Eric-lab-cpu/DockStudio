"""Configuration and result data-models for DockStudio runs.

All objects in this module are plain dataclasses with ``to_dict``/``from_dict``
so a complete run can be persisted to ``config.json`` and re-opened for resume.
"""

from __future__ import annotations

import dataclasses
import json
import os
from dataclasses import dataclass, field
from typing import Any, Optional

from .._version import (
    DEFAULT_CPU,
    DEFAULT_ENERGY_RANGE,
    DEFAULT_EXHAUSTIVENESS,
    DEFAULT_N_POSES,
    DEFAULT_PH,
    DEFAULT_REFINE_EXHAUSTIVENESS,
    DEFAULT_TOP_K,
    DEFAULT_TOP_REFINE,
)

#: standard L-amino acids (3-letter). Modified residues are mapped below.
STANDARD_AAS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

#: modified residue name -> standard residue it should be converted to.
MODIFIED_TO_STANDARD = {
    "TPO": "THR",  # phosphothreonine
    "SEP": "SER",  # phosphoserine
    "PTR": "TYR",  # phosphotyrosine
    "HYP": "PRO",
    "MSE": "MET",
    "CSO": "CYS",
    "SEC": "CYS",
    "KCX": "LYS",
    "MLY": "LYS",
    "LLP": "LYS",
    "PCA": "GLU",
    "HIC": "HIS",
    "HID": "HIS",
    "HIP": "HIS",
}

#: atom names to delete when converting a phospho-residue to its parent amino acid
PHOSPHO_DROP_ATOMS = {"P", "OP1", "OP2", "OP3", "O1P", "O2P", "O3P"}

#: common ion / solvent residue names removed from receptors by default
IONS_SOLVENTS = {
    "HOH", "WAT", "DOD", "NA", "K", "CL", "CA", "MG", "ZN", "MN", "FE",
    "CO", "NI", "CU", "CD", "SO4", "PO4", "GOL", "EDO", "ACT", "PEG", "DMS",
    "FMT", "NO3", "NH4", "CS", "BR", "IOD", "LI", "SR", "BA", "HG", "PT",
}

#: standard DNA/RNA residues (deleted from docking receptors unless the target
#: is a nucleic-acid binding protein whose interaction surface is wanted).
NUCLEIC_RES = {
    "A", "C", "G", "U", "I", "DA", "DC", "DG", "DT", "DU",
    "ADE", "CYT", "GUA", "URI", "THY",
}


@dataclass
class BoxDef:
    """Search-box definition (Angstrom)."""

    method: str = "from_ligand"  # from_ligand | explicit | from_residues | blind
    center: list = field(default_factory=lambda: [0.0, 0.0, 0.0])
    size: list = field(default_factory=lambda: [22.0, 22.0, 22.0])
    basis: str = ""          # human readable rationale
    ligand_ref: Optional[str] = None  # cocrystal ligand id used as center (from_ligand)
    residues: list = field(default_factory=list)  # residue ids used for from_residues

    def to_dict(self) -> dict:
        return dataclasses.asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "BoxDef":
        return cls(**{k: d.get(k) for k in ("method", "center", "size", "basis", "ligand_ref", "residues")})


@dataclass
class InputReceptor:
    """One receptor input file (PDB)."""

    path: str
    name: str = ""            # stem derived from path
    chains: list = field(default_factory=list)  # chains chosen for docking

    def to_dict(self) -> dict:
        return {"path": self.path, "name": self.name, "chains": list(self.chains)}

    @classmethod
    def from_dict(cls, d: dict) -> "InputReceptor":
        return cls(path=d["path"], name=d.get("name", ""), chains=list(d.get("chains", [])))


@dataclass
class InputLigand:
    """One ligand library input file (SDF)."""

    path: str
    name: str = ""            # library name (stem)

    def to_dict(self) -> dict:
        return {"path": self.path, "name": self.name}

    @classmethod
    def from_dict(cls, d: dict) -> "InputLigand":
        return cls(path=d["path"], name=d.get("name", ""))


@dataclass
class RunConfig:
    """Everything needed to (re)run a docking project."""

    receptors: list = field(default_factory=list)   # list[InputReceptor] dicts or objects
    ligands: list = field(default_factory=list)     # list[InputLigand] dicts or objects
    out_dir: str = ""
    title: str = ""

    # boxes: receptor-name -> BoxDef
    boxes: dict = field(default_factory=dict)
    # chain overrides: receptor name -> chains list (also stored per InputReceptor)

    # docking parameters
    exhaustiveness: int = DEFAULT_EXHAUSTIVENESS
    refine_exhaustiveness: int = DEFAULT_REFINE_EXHAUSTIVENESS
    n_poses: int = DEFAULT_N_POSES
    energy_range: float = DEFAULT_ENERGY_RANGE
    top_refine: int = DEFAULT_TOP_REFINE
    top_k: int = DEFAULT_TOP_K
    ph: float = DEFAULT_PH
    cpu: int = DEFAULT_CPU
    time_slice_s: Optional[int] = None   # optional per-call wall clock budget

    # pipeline switches
    run_selfdock: bool = True
    run_plip: bool = True
    run_visualization: bool = True
    run_refine: bool = True
    write_pse: bool = True
    overwrite: bool = False

    # validation helpers
    delete_bad_res: bool = True
    max_box_axis: float = 32.0
    min_box_axis: float = 22.0
    ligand_pad_ang: float = 12.0
    selfdock_pass_rmsd: float = 2.0

    # internal accounting (filled as the run proceeds)
    created_at: str = ""
    updated_at: str = ""

    # -- helpers ---------------------------------------------------------
    @property
    def receptor_names(self) -> list:
        return [r["name"] for r in self.receptors]

    def receptor_by_name(self, name: str) -> InputReceptor:
        for r in self.receptors:
            if r["name"] == name:
                return InputReceptor.from_dict(r)
        raise KeyError(name)

    def to_dict(self) -> dict:
        d = dataclasses.asdict(self)
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "RunConfig":
        known = {f.name for f in dataclasses.fields(cls)}
        kwargs = {k: v for k, v in d.items() if k in known}
        return cls(**kwargs)

    def save(self, path: str) -> None:
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh, ensure_ascii=False, indent=2)

    @classmethod
    def load(cls, path: str) -> "RunConfig":
        with open(path, "r", encoding="utf-8") as fh:
            return cls.from_dict(json.load(fh))


# ---------------------------------------------------------------------------
# Convenience directory helpers
# ---------------------------------------------------------------------------
def project_subdirs(out_dir: str) -> dict:
    return {
        "prepared_receptors": os.path.join(out_dir, "prepared", "receptors"),
        "prepared_ligands": os.path.join(out_dir, "prepared", "ligands"),
        "cocrystal": os.path.join(out_dir, "prepared", "cocrystal"),
        "screening": os.path.join(out_dir, "docking"),
        "refine": os.path.join(out_dir, "refine"),
        "results": os.path.join(out_dir, "results"),
        "plots3d": os.path.join(out_dir, "results", "3D_poses"),
        "composite": os.path.join(out_dir, "results", "composite"),
        "topk_panels": os.path.join(out_dir, "results", "topK_panels"),
        "pml": os.path.join(out_dir, "results", "pymol_scripts"),
        "pse": os.path.join(out_dir, "results", "pse"),
        "lig2d": os.path.join(out_dir, "results", "ligand_2D"),
        "plip": os.path.join(out_dir, "results", "plip"),
        "reports": os.path.join(out_dir, "reports"),
    }


def ensure_project_layout(out_dir: str) -> dict:
    subs = project_subdirs(out_dir)
    for p in subs.values():
        os.makedirs(p, exist_ok=True)
    return subs
