"""Accuracy-assessment report tests (v1.1 feature)."""

import os
import shutil

from dockstudio.core import accuracy, models, utils

PASS_ROW = {
    "poses_available": 9,
    "mode1_rmsd": 1.1,
    "min_rmsd": 0.9,
    "min_rmsd_pose": 3,
    "energy_min_rmsd": -8.0,
    "energy_delta_min_rmsd": 0.4,
    "n_crystal_atoms": 22,
    "pass": True,
    "poses": [{"pose": 1, "affinity": -8.4, "rmsd": 1.1},
              {"pose": 3, "affinity": -8.0, "rmsd": 0.9}],
    "receptor": "rec1",
    "cocrystal_ligand": "rec1_LIG",
    "cocrystal_resname": "LIG",
}


def _cfg(out_dir):
    return models.RunConfig(
        receptors=[{"path": "/tmp/x.pdb", "name": "rec1", "chains": ["A"]},
                   {"path": "/tmp/y.cif", "name": "rec2", "chains": ["A"]}],
        ligands=[{"path": "/tmp/lib.sdf", "name": "demo"}],
        out_dir=out_dir,
        title="accuracy unit test",
        boxes={"rec1": {"method": "from_ligand", "center": [0, 0, 0],
                        "size": [22, 22, 22], "ligand_ref": "LIG"}},
        exhaustiveness=8, refine_exhaustiveness=10,
    )


def test_accuracy_report_honest_no_data(tmp_path):
    out = str(tmp_path / "out")
    os.makedirs(os.path.join(out, "results"), exist_ok=True)
    cfg = _cfg(out)
    summ = accuracy.build_accuracy_report(cfg, out)
    md = os.path.join(out, "reports", "03_accuracy_assessment.md")
    csvp = os.path.join(out, "reports", "accuracy_assessment.csv")
    assert os.path.exists(md) and os.path.exists(csvp)
    assert summ["n_assessable"] == 0
    txt = open(md, encoding="utf-8").read()
    assert "不可评估" in txt or "not assessable" in txt or "无可评估" in txt


def test_accuracy_report_uses_selfdock(tmp_path):
    out = str(tmp_path / "out")
    results_dir = os.path.join(out, "results")
    os.makedirs(results_dir, exist_ok=True)
    cfg = _cfg(out)
    utils.write_json({"selfdock": [PASS_ROW]}, os.path.join(results_dir,
                                                            "selfdock_summary.json"))
    summ = accuracy.build_accuracy_report(cfg, out)
    assert summ["n_assessable"] == 1
    assert summ["n_pass_2A"] == 1
    md = open(os.path.join(out, "reports", "03_accuracy_assessment.md"),
              encoding="utf-8").read()
    assert "rec1" in md
    assert "可接受" in md or "高精度" in md
    # rec2 (no cocrystal) is honestly marked not assessable
    rows = open(os.path.join(out, "reports", "accuracy_assessment.csv"),
                encoding="utf-8-sig").read()
    assert "rec2" in rows
    assert "not_assessable" in rows
