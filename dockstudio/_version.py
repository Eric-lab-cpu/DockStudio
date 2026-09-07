"""Version, branding and application identity for DockStudio."""

# ---------------------------------------------------------------------------
# Branding (Eric Studio)
# ---------------------------------------------------------------------------
BRAND = "Eric Studio"
BRAND_CN = "Eric Studio"
COPYRIGHT = "Copyright (c) 2026 Eric Studio. All rights reserved."
COPYRIGHT_CN = "版权所有 © 2026 Eric Studio。保留所有权利。"

APP_NAME = "DockStudio 分子自动对接平台"
APP_NAME_ASCII = "DockStudio"
APP_TITLE = f"{APP_NAME} · {BRAND}"
__version__ = "2.0.0"
VERSION_LINE = f"v{__version__}"
RELEASE_NOTES = "v2.0.0"

# ---------------------------------------------------------------------------
# Scientific defaults required by the protocol specification
# ---------------------------------------------------------------------------
DEFAULT_EXHAUSTIVENESS = 16
DEFAULT_REFINE_EXHAUSTIVENESS = 32
DEFAULT_N_POSES = 9
DEFAULT_ENERGY_RANGE = 4.0
DEFAULT_TOP_REFINE = 12
DEFAULT_TOP_K = 5
DEFAULT_PH = 7.4
DEFAULT_CPU = 2
SELFDOCK_PASS_RMSD = 2.0  # angstrom: pose is "accurate / PASS" below this RMSD

# Docking-accuracy assessment (self-docking redocking of the native ligand)
DEFAULT_RUN_ACCURACY_REPORT = True
# When reporting RMSD of a docked pose vs the crystal reference, hydrogens are
# never counted; both structures are treated as the same Cartesian frame (no
# re-superposition), which is the standard "redocking RMSD" definition.
ACCURACY_REPORT_BASELINE_STRICT = 1.0   # A: high-confidence docking accuracy
ACCURACY_REPORT_BASELINE_MEDIUM = 2.0  # A: acceptable / typical PASS threshold

# ---------------------------------------------------------------------------
# v2.0.0 — desktop-level, publication-oriented virtual-screening workbench
# ---------------------------------------------------------------------------
# Enrichment (ROC/AUC) validation workflow.  OFF by default: it only runs when
# the user supplies experimental known-actives and decoys (never auto-generated).
DEFAULT_RUN_ENRICHMENT = False
ENRICHMENT_EF_PERCENTILES = (1.0, 5.0)  # EF1% and EF5% reported
# Symmetry-aware RMSD is the primary redocking metric since v2.0; the old
# same-element greedy nearest-neighbour value is always kept as a comparison
# column. When symmetry perception is impossible the symmetry-aware value
# degrades to the greedy one and the report says so.
DEFAULT_USE_SYMMETRY_RMSD = True
# Throughput: process-pool workers for the docking stages.  1 = historical
# serial behaviour (each vina invocation may use `cpu` threads).  When > 1 the
# engine docks `n_workers` ligands in parallel, one Vina thread each, and
# records this in the reports (keeps single-run reproducibility).
DEFAULT_N_WORKERS = 1
# Interactive HTML report + 3D viewer are exported by default in v2.0; they are
# sibling files of the markdown reports and never replace them.
DEFAULT_RUN_HTML_REPORT = True

