#!/usr/bin/env python3
"""Build the demo ligand library SDF for the DockStudio smoke test.

Receptor        : 4DFR (E. coli dihydrofolate reductase, chain A)
Native ligand   : MTX (co-crystal; coordinates taken from the PDB, connectivity
                  from the RCSB Chemical Component Dictionary file MTX.cif)
Other ligands   : trimethoprim (TMP) and pyrimethamine (PYR) - SMILES-derived
                  3D structures, DHFR inhibitors known to bind the same pocket.
"""

from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))

from rdkit import Chem
from rdkit.Chem import AllChem

from dockstudio.core import cofactor

PDB = os.path.join(HERE, "4DFR.pdb")
CCD = os.path.join(HERE, "MTX.cif")
OUT = os.path.join(HERE, "ligand_library.sdf")


def build(out_path: str = OUT) -> str:
    writer = Chem.SDWriter(out_path)
    # 1. native ligand from crystal coordinates (CCD connectivity)
    crystal_sdf = os.path.join(HERE, "_tmp_mtx.sdf")
    stats = cofactor.write_cocrystal_ligand(PDB, "MTX", CCD, crystal_sdf,
                                            chain="A", resseq=161, name="MTX")
    mtx = Chem.MolFromMolFile(crystal_sdf, sanitize=True, removeHs=False)
    mtx.SetProp("_Name", "MTX")
    writer.write(mtx)
    # 2/3. TMP and PYR by SMILES
    for nm, smi in [
        ("TMP", "COc1cc(cc(OC)c1OC)Cc2cnc(nc2N)N"),
        ("PYR", "CCC1=C(C(=NC(=N1)N)N)C2=CC=C(C=C2)Cl"),
    ]:
        m = Chem.AddHs(Chem.MolFromSmiles(smi))
        AllChem.EmbedMolecule(m, randomSeed=7)
        m.SetProp("_Name", nm)
        writer.write(m)
    writer.close()
    if os.path.exists(crystal_sdf):
        os.remove(crystal_sdf)
    return out_path


if __name__ == "__main__":
    print(build())
