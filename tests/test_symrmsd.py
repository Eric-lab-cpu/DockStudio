"""Unit tests for the symmetry-aware RMSD (v2.0 feature #1)."""

import os

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from dockstudio.core import symrmsd


def _mol_from_smiles(smiles: str):
    m = Chem.MolFromSmiles(smiles)
    m = Chem.AddHs(m)
    AllChem.EmbedMolecule(m, randomSeed=0xC0FFEE)
    return m


def _ref_from_smiles(smiles: str):
    m = _mol_from_smiles(smiles)
    els, xyz, orbits = symrmsd.reference_from_mol(m)
    return m, els, xyz, orbits


def test_automorphism_orbits_detects_carboxylate_pair():
    """Two oxygens of a carboxylate must be in the same symmetry orbit."""
    m = Chem.MolFromSmiles("CC(=O)[O-]")
    orbits = symrmsd.automorphism_orbits(m)
    assert orbits is not None
    o_idx = [i for i, a in enumerate(Chem.RemoveHs(m).GetAtoms())
             if a.GetSymbol() == "O"]
    merged = [set(o) for o in orbits]
    # the two O atoms appear together in one orbit
    assert any(set(o_idx) <= g for g in merged), orbits


def test_symmetry_rmsd_improves_on_greedy_for_ambiguous_orbit():
    """On an ambiguous same-element orbit, Hungarian re-assignment must lower RMSD."""
    ref_els = ["O", "O", "O"]
    ref_xyz = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [10.0, 0.0, 0.0]]
    orbits = [[0, 1, 2]]
    pose = [("O", [3.0, 0.0, 0.0]),
            ("O", [4.5, 0.0, 0.0]),
            ("O", [9.5, 0.0, 0.0])]
    res = symrmsd.symmetry_aware_rmsd(ref_els, ref_xyz, orbits, pose)
    assert res["rmsd_greedy"] is not None
    assert res["rmsd"] < res["rmsd_greedy"]
    assert res["symmetry_used"] is True


def test_symmetry_aware_identical_pose_is_zero():
    """A pose identical to the reference has ~0 symmetry-aware RMSD."""
    m, els, xyz, orbits = _ref_from_smiles("CC(=O)[O-]")
    pose = [(els[i], xyz[i]) for i in range(len(els))]
    res = symrmsd.symmetry_aware_rmsd(els, xyz, orbits, pose)
    assert res["rmsd"] is not None and res["rmsd"] < 1e-6


def test_phosphate_swap_does_not_report_fail():
    """Swapping two chemically equivalent phosphate oxygens is not a real error."""
    # reference: tetrahedral P with 4 equivalent-ish O's (methyl phosphate minus H)
    mol = Chem.MolFromSmiles("COP(=O)([O-])[O-]")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    els, xyz, orbits = symrmsd.reference_from_mol(mol)
    # build a pose that is exactly the reference but with the two O atoms of the
    # phosphate swapped in the atom list (coordinate multiset unchanged)
    o_pos = [i for i, e in enumerate(els) if e == "O"][:2]
    pose = [(els[i], xyz[i]) for i in range(len(els))]
    i0, i1 = o_pos[0], o_pos[1]
    pose[i0], pose[i1] = pose[i1], pose[i0]
    res = symrmsd.symmetry_aware_rmsd(els, xyz, orbits, pose)
    assert res["rmsd"] is not None
    # below the 2.0 A PASS threshold => no false FAIL
    assert res["rmsd"] <= 2.0
