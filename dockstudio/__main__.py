"""Entry point for DockStudio.

Usage:
    python -m dockstudio                 # launch the graphical application (GUI)
    python -m dockstudio --version       # print version
    python -m dockstudio run <project> [options]
                                         # headless CLI (v2.0 escape hatch)

The ``run`` subcommand reads the exact same project file the GUI produces
(``.dsproj`` or the legacy ``config.json`` inside an output directory) and
executes the engine pipeline headlessly. DockStudio stays GUI-first; this CLI is
an optional automation / cluster escape hatch and never changes GUI behaviour.
"""

from __future__ import annotations

import argparse
import multiprocessing
import sys

from dockstudio._version import __version__


def _run_headless(argv) -> int:
    from dockstudio.core import pipeline, project

    p = argparse.ArgumentParser(prog="dockstudio run", description=__doc__)
    p.add_argument("project", help="path to a .dsproj project file or a legacy "
                                   "config.json inside an output directory")
    p.add_argument("--out", default=None,
                   help="override the output directory (default: stored in project)")
    p.add_argument("--phases", default=None,
                   help="comma-separated phase subset, e.g. "
                        "prepare,screen,refine,selfdock,analyze,enrich,visualize,report,html")
    p.add_argument("--time", type=float, default=None,
                   help="optional wall-clock time-slice in seconds per run")
    args = p.parse_args(argv)

    try:
        cfg = project.load_run_config(args.project)
    except Exception as e:
        print(f"[dockstudio] cannot load project {args.project}: {e}", file=sys.stderr)
        return 2
    if args.out:
        cfg.out_dir = args.out
    if not cfg.out_dir:
        print("[dockstudio] no output directory set (pass --out).", file=sys.stderr)
        return 2

    phases = None
    if args.phases:
        phases = [s.strip() for s in args.phases.split(",") if s.strip()]

    def _log(m):
        print(m, flush=True)

    print(f"[dockstudio] DockStudio v{__version__} headless run")
    print(f"[dockstudio] project: {args.project}")
    print(f"[dockstudio] out_dir : {cfg.out_dir}")
    if phases:
        print(f"[dockstudio] phases  : {phases}")
    try:
        summary = pipeline.run_project(cfg, on_log=_log, phases=phases,
                                       time_slice_s=args.time)
    except Exception as e:
        print(f"[dockstudio] pipeline FAILED: {e}", file=sys.stderr)
        return 1
    inter = bool(summary.get("interrupted"))
    print(f"[dockstudio] finished. interrupted={inter} phases={summary.get('phases')}")
    return 0 if not inter else 3


def main(argv=None) -> int:
    multiprocessing.freeze_support()
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        from dockstudio.gui.app import main as gui_main
        return gui_main(sys.argv)
    if argv[0] in ("--version", "-V"):
        from dockstudio._version import VERSION_LINE
        print(VERSION_LINE)
        return 0
    if argv[0] == "run":
        return _run_headless(argv[1:])
    if argv[0] in ("--help", "-h", "help"):
        print(__doc__)
        return 0
    # fall back to GUI (legacy behaviour)
    from dockstudio.gui.app import main as gui_main
    return gui_main(sys.argv)


if __name__ == "__main__":
    raise SystemExit(main())
