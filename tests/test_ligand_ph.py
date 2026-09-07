"""Unit tests for the documented simple pH protonation rules (v1.1+/v2.0).

These rules are deliberately crude (not a pKa predictor) but they must actually
do what the reports claim: carboxyl acids deprotonate at pH >= 5, basic amines
protonate at pH <= 8.5.
"""

from rdkit import Chem

from dockstudio.core import ligand


def _apply(smiles: str, ph: float = 7.4):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    return ligand.apply_ph_rules(mol, ph)


def _smiles_of(mol) -> str:
    return Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(mol)))


def test_carboxyl_deprotonated_at_physiological_ph():
    m = _apply("CC(=O)O")
    assert "[O-]" in _smiles_of(m)
    assert "carboxyl deprotonated" in m.GetProp("_DockStudio_ph_notes")


def test_two_carboxyls_both_deprotonated():
    m = _apply("OC(=O)C(=O)O")
    s = _smiles_of(m)
    assert s.count("[O-]") == 2


def test_no_carboxyl_deprotonation_at_acidic_ph():
    # pH 4 (< 5 cutoff) -> acetic acid stays neutral
    m = _apply("CC(=O)O", ph=4.0)
    assert "[O-]" not in _smiles_of(m)


def test_basic_amine_protonated_at_physiological_ph():
    m = _apply("CCN")
    assert "+" in _smiles_of(m)          # ammonium form
    assert "basic amine protonated" in m.GetProp("_DockStudio_ph_notes")


def test_zwitterionic_residue():
    m = _apply("NCCCCC(N)C(=O)O")  # lysine-like
    s = _smiles_of(m)
    assert s.count("[NH3+]") >= 1
    assert "[O-]" in s


def test_imidazole_not_treated_as_basic_amine():
    # aromatic ring N is not an sp3 basic amine under the documented rules
    m = _apply("c1cnc[nH]1")
    assert "+" not in _smiles_of(m)
