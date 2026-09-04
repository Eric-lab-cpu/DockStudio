"""Search-box definition helpers (section 5.4 decision tree)."""

from __future__ import annotations

from typing import List, Optional

from . import structure as st
from .models import BoxDef
from .utils import distance


def _pad_size(extent_axis: float, pad: float, min_axis: float, max_axis: float) -> float:
    v = max(min_axis, extent_axis + pad)
    return min(v, max_axis)


def box_from_ligand(recs: List[st.AtomRec],
                    min_axis: float = 22.0,
                    pad: float = 12.0,
                    max_axis: float = 32.0,
                    ligand_tag: str = "") -> BoxDef:
    center = st.heavy_centroid(recs)
    lo, hi = st.coord_extent(recs)
    size = []
    for i in range(3):
        size.append(_pad_size(hi[i] - lo[i], pad, min_axis, max_axis))
    return BoxDef(method="from_ligand", center=[round(c, 3) for c in center],
                  size=[round(s, 3) for s in size],
                  ligand_ref=ligand_tag or None,
                  basis=f"centered on co-crystal ligand '{ligand_tag}' heavy-atom centroid; "
                        f"axis = max({min_axis}, ligand extent + {pad}), capped at {max_axis}")


def box_from_center(center: List[float], size: List[float], basis: str = "user-specified") -> BoxDef:
    return BoxDef(method="explicit", center=[round(c, 3) for c in center],
                  size=[round(s, 3) for s in size], basis=basis)


def box_from_residues(recs: List[st.AtomRec],
                      residues: List[str],
                      pad: float = 8.0,
                      min_axis: float = 22.0,
                      max_axis: float = 32.0) -> BoxDef:
    """Box centered on selected residue heavy atoms, sized from their spread."""
    sel = [r for r in recs
           if f"{r.resname}{r.resseq}" in residues or f"{r.chain}:{r.resname}{r.resseq}" in residues
           or f"{r.resseq}" in residues]
    if not sel:
        raise ValueError("no atoms matched the supplied residue list")
    center = st.heavy_centroid(sel)
    lo, hi = st.coord_extent(sel, pad=pad)
    size = [max(min_axis, hi[i] - lo[i]) for i in range(3)]
    size = [min(max_axis, s) for s in size]
    return BoxDef(method="from_residues", center=[round(c, 3) for c in center],
                  size=[round(s, 3) for s in size], residues=list(residues),
                  basis=f"centered on residues {','.join(residues)} ({len(sel)} heavy atoms); "
                        f"padding {pad} A per axis, capped at {max_axis}")


def box_blind(recs: List[st.AtomRec], min_axis: float = 22.0, max_axis: float = 48.0) -> BoxDef:
    center = st.heavy_centroid(recs)
    lo, hi = st.coord_extent(recs, pad=6.0)
    size = [max(min_axis, min(max_axis, hi[i] - lo[i])) for i in range(3)]
    return BoxDef(method="blind", center=[round(c, 3) for c in center],
                  size=[round(s, 3) for s in size],
                  basis="whole-protein blind docking box (EXPLORATORY, high false-positive risk)")


def box_auto_decision(recs: List[st.AtomRec],
                      cocrystal: Optional[str] = None,
                      cocrystal_records: Optional[List[st.AtomRec]] = None,
                      user_box: Optional[BoxDef] = None,
                      user_residues: Optional[List[str]] = None) -> BoxDef:
    """Protocol 5.4 decision tree.

    Priority: user box > co-crystal ligand in chosen chain > residue list >
    blind (exploratory). This helper returns the box; scientific caveats must
    still be surfaced by the caller.
    """
    if user_box is not None and user_box.method in ("explicit", "from_ligand"):
        return user_box
    if user_residues:
        return box_from_residues(recs, user_residues)
    if cocrystal and cocrystal_records:
        return box_from_ligand(cocrystal_records, ligand_tag=cocrystal)
    return box_blind(recs)
