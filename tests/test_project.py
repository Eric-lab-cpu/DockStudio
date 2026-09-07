"""Unit tests for .dsproj save/load and zip export (v2.0 features #18)."""

import os
import zipfile

from dockstudio.core import models, project


def _cfg(out_dir: str) -> models.RunConfig:
    return models.RunConfig(
        receptors=[{"path": "/x/1.pdb", "name": "R1", "chains": ["A"]}],
        ligands=[{"path": "/x/l.sdf", "name": "lib"}],
        out_dir=out_dir,
        title="t",
        top_k=3,
        n_workers=2,
        run_enrichment=True,
        actives_path="/x/act.smi",
        decoys_path="/x/dec.smi",
    )


def test_dsproj_roundtrip(tmp_path):
    cfg = _cfg(str(tmp_path / "out"))
    p = tmp_path / "proj.dsproj"
    project.save_project(cfg, str(p))
    cfg2 = project.load_run_config(str(p))
    assert cfg2.title == cfg.title
    assert cfg2.top_k == 3
    assert cfg2.n_workers == 2
    assert cfg2.run_enrichment is True
    assert cfg2.receptor_names == ["R1"]


def test_export_zip_packs_real_files(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "config.json").write_text("{}", encoding="utf-8")
    (out / "results").mkdir()
    (out / "results" / "final_top5.csv").write_text("receptor,ligand\nR1,L1\n",
                                                    encoding="utf-8")
    dest = tmp_path / "export.zip"
    created = project.export_run_zip(str(out), str(dest))
    assert os.path.exists(created)
    with zipfile.ZipFile(created) as zf:
        names = zf.namelist()
    assert any(n.endswith("config.json") for n in names)
    assert any(n.endswith("final_top5.csv") for n in names)
    # the archive should NOT contain the zip itself (no recursion)
    assert not any(n.endswith("export.zip") for n in names)


def test_load_run_config_reads_legacy_config(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    cfg = _cfg(str(out))
    cfg.save(str(out / "config.json"))
    cfg2 = project.load_run_config(str(out / "config.json"))
    assert cfg2.top_k == 3
