"""mmCIF(.cif) receptor-input tests (v1.1 feature)."""

import os

import pytest

from dockstudio.core import inventory, structure as st

DATA = os.path.join(os.path.dirname(__file__), "..", "examples", "demo", "4DFR.pdb")

gemmi = pytest.importorskip("gemmi")


@pytest.fixture(scope="module")
def cif_path(tmp_path_factory):
    p = tmp_path_factory.mktemp("cif") / "4DFR.cif"
    doc = gemmi.read_structure(DATA)
    doc.make_mmcif_document().write_file(str(p))
    return str(p)


def test_parse_structure_cif_matches_pdb(cif_path):
    pdb_recs = st.parse_structure(DATA)
    cif_recs = st.parse_structure(cif_path)
    assert len(pdb_recs) == len(cif_recs) > 1000
    # same chain assignment
    pdb_inv = st.inventory_pdb(DATA, name="p")
    cif_inv = st.inventory_pdb(cif_path, name="c")
    assert sorted(pdb_inv.chains.keys()) == sorted(cif_inv.chains.keys())
    for c, ci in cif_inv.chains.items():
        assert ci.polymer_res > 100
        ligs = [l["resname"] for l in ci.het_ligands]
        assert "MTX" in ligs


def test_inventory_receptor_cif(cif_path):
    rows = inventory.inventory_receptor(cif_path, "4DFR", ["A"])
    assert rows and rows[0]["ligand"] == "MTX"
    assert rows[0]["selected_for_docking"] == "yes"


def test_coords_of_ligand_cif(cif_path):
    # coordinates must be identical in the same frame as the PDB input
    pdb_mtx = st.coords_of_ligand(DATA, "MTX", chain="A")
    cif_mtx = st.coords_of_ligand(cif_path, "MTX", chain="A")
    assert len(pdb_mtx) == len(cif_mtx) > 0
    for a, b in zip(pdb_mtx, cif_mtx):
        assert a.x == pytest.approx(b.x, abs=1e-3)
        assert a.y == pytest.approx(b.y, abs=1e-3)
        assert a.z == pytest.approx(b.z, abs=1e-3)


def test_clean_receptor_cif_matches_pdb(cif_path, tmp_path):
    p_out = st.clean_receptor(DATA, str(tmp_path / "p.clean.pdb"), keep_chains=["A"])
    c_out = st.clean_receptor(cif_path, str(tmp_path / "c.clean.pdb"), keep_chains=["A"])
    assert p_out.kept == c_out.kept > 1000
    assert p_out.dropped == c_out.dropped
