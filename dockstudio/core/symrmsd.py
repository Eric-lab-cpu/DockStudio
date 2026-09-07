"""Symmetry-aware RMSD for redocking validation (v2.0).

The historical DockStudio self-docking RMSD used a same-element greedy
nearest-neighbour match without any symmetry correction.  For ligands that
contain chemically equivalent atoms (carboxylate oxygens, phosphate oxygens,
the two halves of a symmetric aromatic ring, ...) an arbitrary "swap" of two
equivalent atoms between the crystal reference and the docked pose is not a
real geometric error, yet the plain greedy value would count it as one and can
turn a correct pose into a false FAIL.

This module computes a symmetry-aware RMSD:

* equivalence classes (orbits) of reference heavy atoms are derived from a
  **Weisfeiler-Lehman colour-refinement** partition of the reference graph
  (bond orders and charges ignored - so resonance-equivalent atoms such as
  carboxylate / phosphate oxygens are treated as one class; hydrogens removed
  before the computation);
* an element-preserving greedy match seeds the reference<->pose correspondence;
* within every symmetry orbit the pose atoms assigned to that orbit are then
  re-permuted with the Hungarian (linear-sum-assignment) algorithm to minimise
  the RMSD - i.e. swapping two chemically equivalent atoms costs nothing;
* the plain greedy RMSD is always returned alongside as a comparison column.

Honesty / method disclosure
---------------------------
* Atoms from different orbits are never swapped, so the correction cannot
  fabricate an artificially low RMSD by re-mapping chemically distinct atoms.
* If symmetry perception is not possible (molecule too large) the module
  degrades to the greedy value and reports ``symmetry_used = False`` so reports
  stay truthful.
* The method is intended for *redocking* (same molecule, same Cartesian
  frame). It is not a superposition/alignment RMSD.

No GUI/tkinter dependency: safe for headless testing.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

try:  # scipy is a declared runtime dependency (packaging/requirements-*.txt)
    from scipy.optimize import linear_sum_assignment  # type: ignore
    _HAVE_SCIPY = True
except Exception:  # pragma: no cover - defensive
    _HAVE_SCIPY = False

from rdkit import Chem

#: maximum heavy-atom count for automorphism enumeration (safety bound)
_MAX_HEAVY_ORBITS = 150
#: maximum number of self-substructure matches used for orbit detection
_MAX_AUTOMORPHISMS = 2000


def automorphism_orbits(mol: Chem.Mol) -> Optional[List[List[int]]]:
    """Return symmetry-equivalence orbits of the *heavy* atoms of ``mol``.

    Each orbit is a list of heavy-atom indices (indices refer to a copy of
    ``mol`` with hydrogens removed, i.e. the same ordering as the arrays built
    by :func:`reference_from_mol`).  Returns ``None`` when the molecule is too
    large for a cheap computation (safety bound).

    Equivalence is computed with the **Weisfeiler-Lehman colour-refinement**
    algorithm on the *topological graph* (bond orders and formal charges are
    ignored).  This deliberately treats resonance-equivalent atoms (the two
    oxygens of a carboxylate, the four oxygens of a phosphate) as one symmetry
    class even when the reference file depicts one localised single/double bond
    - which is the correct behaviour for redocking RMSD symmetry correction.
    Colour refinement is a recognised, deterministic approximation of the
    automorphism orbits used for molecular symmetry perception.
    """
    try:
        hm = Chem.RemoveHs(mol)
    except Exception:
        return None
    if hm.GetNumHeavyAtoms() > _MAX_HEAVY_ORBITS or hm.GetNumAtoms() > 220:
        return None
    heavy = [a.GetIdx() for a in hm.GetAtoms() if a.GetAtomicNum() > 1]
    k = len(heavy)
    if k == 0:
        return []
    pos_of = {aidx: p for p, aidx in enumerate(heavy)}
    adj: List[List[int]] = [[] for _ in range(k)]
    for b in hm.GetBonds():
        i = pos_of.get(b.GetBeginAtomIdx())
        j = pos_of.get(b.GetEndAtomIdx())
        if i is not None and j is not None and i != j:
            adj[i].append(j)
            adj[j].append(i)
    # initial colour = element symbol
    colours = [hm.GetAtomWithIdx(aidx).GetSymbol() for aidx in heavy]
    n_iter = 0
    max_iter = k + 2
    while n_iter < max_iter:
        n_iter += 1
        # signature = (own colour, sorted neighbour colours)
        sig_map = {}
        next_colours = []
        for p in range(k):
            sig = (colours[p], tuple(sorted(colours[q] for q in adj[p])))
            if sig not in sig_map:
                sig_map[sig] = len(sig_map)
            next_colours.append(sig_map[sig])
        if next_colours == colours:
            break
        colours = next_colours
    groups: Dict[int, List[int]] = {}
    for p in range(k):
        groups.setdefault(colours[p], []).append(p)
    orbits = sorted((sorted(v) for v in groups.values()), key=lambda v: v[0])
    return orbits


def reference_from_mol(mol: Chem.Mol):
    """Extract heavy-atom element/coordinate arrays + automorphism orbits.

    Returns ``(elements, xyz, orbits)`` where ``elements[i]`` is the element
    symbol, ``xyz[i]`` the [x,y,z] list and ``orbits`` indexes those arrays.
    ``orbits`` is ``None`` when symmetry perception is not available.
    """
    hm = Chem.RemoveHs(Chem.Mol(mol))
    if hm.GetNumConformers() == 0:
        raise ValueError("reference molecule has no coordinates")
    conf = hm.GetConformer()
    elements: List[str] = []
    xyz: List[List[float]] = []
    for a in hm.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        elements.append(a.GetSymbol())
        xyz.append([p.x, p.y, p.z])
    orbits = automorphism_orbits(hm)
    return elements, xyz, orbits


def reference_from_sdf(sdf_path: str):
    """Load a single-molecule SDF and return ``(elements, xyz, orbits)``."""
    mol = Chem.MolFromMolFile(sdf_path, sanitize=False, removeHs=True)
    if mol is None:
        return [], [], None
    return reference_from_mol(mol)


def greedy_rmsd_and_pairs(
    ref_elements: Sequence[str],
    ref_xyz: Sequence[Sequence[float]],
    pose_atoms: Sequence[Tuple[str, Sequence[float]]],
) -> Tuple[Optional[float], List[Tuple[int, int]]]:
    """Same-element greedy nearest-neighbour matching (historical method).

    ``pose_atoms`` is a list of ``(element, [x,y,z])``. Returns
    ``(rmsd, pairs)`` with pairs ``(ref_index, pose_index)``; ``rmsd`` is
    ``None`` when no pair could be formed.
    """
    n_ref = len(ref_elements)
    used = [False] * n_ref
    pairs: List[Tuple[int, int]] = []
    total = 0.0
    for pi, (pel, pxyz) in enumerate(pose_atoms):
        best_i = -1
        best_d2 = float("inf")
        for ri in range(n_ref):
            if used[ri] or ref_elements[ri] != pel:
                continue
            d2 = sum((a - b) ** 2 for a, b in zip(ref_xyz[ri], pxyz))
            if d2 < best_d2:
                best_d2 = d2
                best_i = ri
        if best_i < 0:
            continue
        used[best_i] = True
        pairs.append((best_i, pi))
        total += best_d2
    if not pairs:
        return None, []
    return float(np.sqrt(total / len(pairs))), pairs


def _rmsd_of_pairs(ref_xyz, pose_xyz, pairs) -> float:
    if not pairs:
        return float("nan")
    s = 0.0
    for ri, pi in pairs:
        s += sum((a - b) ** 2 for a, b in zip(ref_xyz[ri], pose_xyz[pi]))
    return float(np.sqrt(s / len(pairs)))


def _refine_orbits(ref_xyz, ref_elements, orbits, pose_xyz, pairs):
    """Re-permute pose atoms inside each symmetry orbit (Hungarian)."""
    if not _HAVE_SCIPY or not orbits:
        return pairs
    # pose index currently paired to each reference index
    pair_of_ref: Dict[int, int] = {ri: pi for ri, pi in pairs}
    new_pairs = list(pairs)
    for orbit in orbits:
        members = [ri for ri in orbit if ri in pair_of_ref]
        if len(members) < 2:
            continue
        k = len(members)
        pos_pairs = [pair_of_ref[ri] for ri in members]
        # cost matrix: members (ref) x those pose atoms
        cost = np.zeros((k, k), dtype=float)
        for a, ri in enumerate(members):
            for b, pi in enumerate(pos_pairs):
                cost[a, b] = sum(
                    (x - y) ** 2 for x, y in zip(ref_xyz[ri], pose_xyz[pi]))
        row, col = linear_sum_assignment(cost)
        # col[b] = pose position assigned to ref row b (i.e. members[row_idx])
        for r_idx, c_idx in zip(row, col):
            # members[row[r]] gets pose atom pos_pairs[col[c]]
            assigned_pi = pos_pairs[c_idx]
            ri = members[r_idx]
            # replace old pair (ri -> old pi) by (ri -> assigned_pi)
            for idx, (old_ri, old_pi) in enumerate(new_pairs):
                if old_ri == ri:
                    new_pairs[idx] = (ri, assigned_pi)
                    break
    return new_pairs


def symmetry_aware_rmsd(
    ref_elements: Sequence[str],
    ref_xyz: Sequence[Sequence[float]],
    orbits: Optional[List[List[int]]],
    pose_atoms: Sequence[Tuple[str, Sequence[float]]],
) -> dict:
    """Compute greedy + symmetry-aware RMSD between a reference and a pose.

    ``ref_elements`` / ``ref_xyz`` come from :func:`reference_from_sdf`;
    ``orbits`` are the symmetry-equivalence orbits (or ``None``).  ``pose_atoms``
    is the heavy-atom list ``(element, [x,y,z])`` of one docked pose.

    Returns a dict with keys: ``rmsd`` (symmetry-aware, primary), ``rmsd_greedy``,
    ``n_matched``, ``n_ref``, ``n_pose``, ``symmetry_used`` and ``method``.
    When no match is possible the RMSD values are ``None``.
    """
    pose_elements = [p[0] for p in pose_atoms]
    pose_xyz = [p[1] for p in pose_atoms]
    greedy_rmsd, pairs = greedy_rmsd_and_pairs(ref_elements, ref_xyz, pose_atoms)
    if greedy_rmsd is None:
        return {"rmsd": None, "rmsd_greedy": None, "n_matched": 0,
                "n_ref": len(ref_elements), "n_pose": len(pose_atoms),
                "symmetry_used": False, "method": "no match"}
    use_orbits = bool(orbits) and any(len(o) > 1 for o in orbits)
    if use_orbits:
        pairs = _refine_orbits(ref_xyz, ref_elements, orbits, pose_xyz, pairs)
    rmsd = _rmsd_of_pairs(ref_xyz, pose_xyz, pairs)
    if use_orbits:
        method = ("symmetry-aware (graph colour-refinement equivalence classes + "
                  "Hungarian re-assignment)")
    elif orbits:
        method = ("same-element greedy nearest-neighbour "
                  "(symmetry orbits detected but all trivial)")
    else:
        method = ("same-element greedy nearest-neighbour "
                  "(symmetry perception not applied)")
    return {
        "rmsd": float(rmsd),
        "rmsd_greedy": float(greedy_rmsd),
        "n_matched": len(pairs),
        "n_ref": len(ref_elements),
        "n_pose": len(pose_atoms),
        "symmetry_used": use_orbits,
        "method": method,
    }
