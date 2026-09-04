import os

import pytest

from dockstudio.core import structure as st

DATA = os.path.join(os.path.dirname(__file__), "..", "examples", "demo", "4DFR.pdb")


@pytest.fixture(scope="module")
def pdb(tmp_path_factory):
    out = tmp_path_factory.mktemp("clean")
    act = st.clean_receptor(DATA, str(out / "4dfr.clean.pdb"), keep_chains=["A"])
    return out / "4dfr.clean.pdb", act


def test_clean_writes_valid_pdb(pdb):
    path, act = pdb
    n = sum(1 for l in open(path) if l.startswith("ATOM"))
    assert n == act.kept
    assert n > 1000


def test_clean_drops_ligand_water(pdb):
    path, act = pdb
    assert "MTX" in act.dropped
    assert act.dropped.get("HOH", 0) > 0


def test_parse_records(pdb):
    path, _ = pdb
    recs = st.parse_pdb(str(path))
    assert recs
    # all kept heavy atoms have correct element column
    for r in recs[:200]:
        assert r.element in ("C", "N", "O", "S")


def test_centroid_and_extent(pdb):
    path, _ = pdb
    recs = st.parse_pdb(str(path))
    c = st.heavy_centroid(recs)
    lo, hi = st.coord_extent(recs)
    assert lo[0] <= c[0] <= hi[0]
