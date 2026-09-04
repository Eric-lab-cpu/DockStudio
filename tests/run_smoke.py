#!/usr/bin/env python3
"""End-to-end smoke test for DockStudio on the 4DFR demo (headless engine).

Usage:
    python tests/run_smoke.py [time_slice_s] [--out DIR] [--phases p1,p2,...]

Run repeatedly; docking stages stop after ``time_slice_s`` seconds and resume
from where they stopped (``done.flag`` files). When the run returns
``interrupted=false`` every requested phase is complete.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, ROOT)

from dockstudio.core import models
from dockstudio.core.pipeline import PHASES, run_project

DEMO = os.path.join(ROOT, "examples", "demo")


def make_cfg(out_dir: str) -> models.RunConfig:
    lib = os.path.join(DEMO, "ligand_library.sdf")
    if not os.path.exists(lib):
        sys.path.insert(0, os.path.join(DEMO))
        import make_library
        make_library.build(lib)
    cfg = models.RunConfig(
        receptors=[{"path": os.path.join(DEMO, "4DFR.pdb"), "name": "4DFR",
                    "chains": ["A"]}],
        ligands=[{"path": lib, "name": "demo"}],
        out_dir=out_dir,
        title="4DFR DHFR smoke test",
        exhaustiveness=8,        # reduced for smoke test (default 16)
        refine_exhaustiveness=10,
        n_poses=9,
        energy_range=4.0,
        top_refine=3,
        top_k=3,
        cpu=1,
        run_selfdock=True,
        run_plip=True,
        run_visualization=True,
        run_refine=True,
        write_pse=False,
        ph=7.4,
        delete_bad_res=True,
    )
    return cfg


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("slice_s", nargs="?", type=float, default=150.0)
    ap.add_argument("--out", default=os.path.join(ROOT, "examples", "demo_out"))
    ap.add_argument("--phases", default=",".join(PHASES))
    args = ap.parse_args(argv)
    cfg = make_cfg(args.out)
    phases = [p for p in args.phases.split(",") if p]

    def log(m):
        print(m, flush=True)

    summary = run_project(cfg, on_log=log, phases=phases, time_slice_s=args.slice_s)
    print("\n=== SMOKE SUMMARY ===")
    print(json.dumps(summary, ensure_ascii=False, indent=1))
    return 0 if not summary.get("interrupted") else 2


if __name__ == "__main__":
    raise SystemExit(main())
