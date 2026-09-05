"""PDBQT pose-parsing tests, incl. AutoDock type->element mapping (v1.1 fix)."""

import os

from dockstudio.core import docking


def _write(txt, tmp_path):
    p = os.path.join(str(tmp_path), "pose.pdbqt")
    with open(p, "w") as fh:
        fh.write(txt)
    return p


def test_element_mapping_autodock_types(tmp_path):
    txt = (
        "MODEL     1\n"
        "REMARK VINA RESULT:   -8.1  0.000  0.000\n"
        "ATOM      1  C1  LIG A   1       1.000   2.000   3.000  1.00  0.00     0    A\n"
        "ATOM      2  N2  LIG A   1       1.000   2.000   4.000  1.00  0.00     0    NA\n"
        "ATOM      3  O3  LIG A   1       1.000   2.000   5.000  1.00  0.00     0    OA\n"
        "ATOM      4  CL1 LIG A   1       1.000   2.000   6.000  1.00  0.00     0   Cl\n"
        "ATOM      5  H4  LIG A   1       1.000   2.000   7.000  1.00  0.00     0    HD\n"
        "ENDMDL\n"
    )
    poses = docking.parse_pdbqt_poses(_write(txt, tmp_path))
    assert len(poses) == 1
    elems = [(a.name, a.element) for a in poses[0].atoms]
    assert elems == [("C1", "C"), ("N2", "N"), ("O3", "O"),
                     ("CL1", "Cl"), ("H4", "H")]
    assert poses[0].affinity == -8.1


def test_no_model_returns_empty(tmp_path):
    assert docking.parse_pdbqt_poses(_write("ATOM ... no model\n", tmp_path)) == []
