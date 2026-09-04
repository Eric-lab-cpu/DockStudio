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
__version__ = "1.0.0"

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
SELFDOCK_PASS_RMSD = 2.0  # angstrom
