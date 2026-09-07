"""Unit tests for the Chemical Component Dictionary (CCD) parser (cofactor)."""

import os

import pytest

from dockstudio.core import cofactor

DEMO_CCD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                        "examples", "demo", "MTX.cif")


def test_parse_ccd_bond_orders_are_kept():
    """CCD bond order codes (DOUB/AROM/...) must not collapse to single bonds."""
    if not os.path.exists(DEMO_CCD):
        pytest.skip("demo MTX.cif not present")
    ccd = cofactor.parse_ccd(DEMO_CCD)
    assert len(ccd.bonds) > 20
    orders = {b.order for b in ccd.bonds}
    assert 2 in orders, "expected at least one double bond in MTX (carbonyls)"
    assert any(b.aromatic for b in ccd.bonds), "expected aromatic bonds in MTX"


def test_parse_ccd_elements_present():
    if not os.path.exists(DEMO_CCD):
        pytest.skip("demo MTX.cif not present")
    ccd = cofactor.parse_ccd(DEMO_CCD)
    elements = {a.element for a in ccd.atoms}
    assert {"C", "N", "O"} <= elements
