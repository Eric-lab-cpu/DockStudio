"""Unit tests for SMILES/CSV ligand parsing (v2.0 feature #2)."""

import os
import shutil

import pytest
from rdkit import Chem

from dockstudio.core import ligand


def test_read_smiles_records_plain(tmp_path):
    p = tmp_path / "lig.smi"
    p.write_text("# comment\n"
                 "CC(=O)O\tacetic\n"
                 "c1ccccc1 benzene\n"
                 "\n"
                 "CCN  ethylamine\n", encoding="utf-8")
    recs = ligand.read_smiles_records(str(p))
    assert len(recs) == 3
    assert recs[0]["name"] == "acetic"
    assert recs[0]["smiles"] == "CC(=O)O"
    assert recs[2]["name"] == "ethylamine"


def test_read_smiles_records_csv(tmp_path):
    p = tmp_path / "lig.csv"
    p.write_text("ID,SMILES\nmol1,CC(=O)O\nmol2,c1ccccc1\n", encoding="utf-8")
    recs = ligand.read_smiles_records(str(p))
    assert len(recs) == 2
    assert recs[0]["name"] == "mol1"
    assert recs[1]["smiles"] == "c1ccccc1"


def test_is_smiles_source():
    assert ligand.is_smiles_source("a.smi")
    assert ligand.is_smiles_source("a.csv")
    assert not ligand.is_smiles_source("a.sdf")


def test_prepare_smiles_library_ok_or_skips(tmp_path):
    if not (shutil.which("mk_prepare_ligand.py") or shutil.which("mk_prepare_ligand")):
        pytest.skip("meeko mk_prepare_ligand not available")
    p = tmp_path / "l.smi"
    p.write_text("CC(=O)O acetic\nc1ccccc1 benzene\n", encoding="utf-8")
    outdir = tmp_path / "prep"
    entries = ligand.prepare_smiles_library(str(p), str(outdir), ph=7.4,
                                            library_name="lib")
    ok = [e for e in entries if e.pdbqt_ok]
    assert len(ok) >= 1
    for e in ok:
        assert e.mol_source == "smiles"
        assert e.smiles
        assert os.path.exists(e.pdbqt)
