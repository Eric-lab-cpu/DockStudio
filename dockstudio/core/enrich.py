"""Enrichment (ROC / AUC / EF) validation workflow (v2.0, feature #16).

Goal
----
Turn a virtual screen into a *validation* experiment: when the user supplies a
known-active set and a (putative) decoy / inactive set that were part of the
docked ligand library, DockStudio measures how well its ranking separates the
two classes. This is the evidence reviewers of a docking paper most often ask
for, and it is computed from *real* docking scores only.

Method (reported honestly)
--------------------------
* Labels are attached to docked molecules by **canonical isomeric SMILES of the
  neutralised molecule** (actives / decoys files may be SDF, SMILES or CSV).
  Neutralisation makes the matching robust to DockStudio's pH protonation
  (R-COO- / R-NH3+ vs the neutral label), otherwise every ionizable compound
  would be falsely counted as label attrition.  Molecules that docked but were
  not labelled are ignored (they carry no ground truth); labels that did not
  produce a docking score are counted and reported as attrition.
* The per-molecule score is the **screening mode-1 Vina affinity** - one single,
  consistent docking stage across the whole library.  Refined shortlist energies
  are deliberately not mixed in (refinement re-docks only the best-scoring
  subset, so mixing would inflate the enrichment statistics).
  Scores are ranked best-first (more negative affinity = better).
* ROC curve / AUC: standard rank-based (Mann-Whitney) AUC with average ranks
  for ties. EF1%/EF5% = (fraction of actives found in the top 1%/5% of the
  ranked docked list) / (0.01 / 0.05).
* Nothing is auto-generated: if the user does not provide actives and decoys
  this module does nothing.

No matplotlib dependency: the ROC figure is emitted as an inline SVG and the
ROC points are also stored as CSV/JSON so the numbers can be re-plotted.
"""

from __future__ import annotations

import csv
import os
from typing import Dict, List, Optional, Sequence

import numpy as np

try:
    from scipy.stats import rankdata  # type: ignore
    _HAVE_SCIPY = True
except Exception:  # pragma: no cover - defensive
    _HAVE_SCIPY = False

from rdkit import Chem

from . import batch, models, utils
from .ligand import is_smiles_source, read_smiles_records
from .models import RunConfig, project_subdirs

EF_PERCENTILES = (1.0, 5.0)  # EF1%, EF5%


# ---------------------------------------------------------------------------
# molecule identity (canonical SMILES)
# ---------------------------------------------------------------------------
def _make_uncharger():
    """Return an RDKit Uncharger across rdkit builds/module layouts."""
    try:
        from rdkit.Chem import rdMolStandardize
        return rdMolStandardize.Uncharger()
    except Exception:
        from rdkit.Chem import MolStandardize
        return MolStandardize.rdMolStandardize.Uncharger()


_UNCHARGER = None


def canonical_smiles(mol: Chem.Mol) -> str:
    """Canonical isomeric SMILES of the *neutralised* molecule.

    DockStudio protonates ligands before docking (simple pH rules: carboxylates
    become R-COO-, basic amines R-NH3+).  A docked molecule must therefore be
    matched to its label on the neutralised heavy-atom graph, otherwise every
    ionizable compound would be counted as "not docked / label lost" even though
    it clearly is present.  Neutralisation keeps stereochemistry (isomeric
    SMILES) but removes formal charges before the canonical string is built.
    """
    try:
        m = Chem.RemoveHs(Chem.Mol(mol))
        global _UNCHARGER
        if _UNCHARGER is None:
            _UNCHARGER = _make_uncharger()
        try:
            m = _UNCHARGER.uncharge(m)
        except Exception:
            pass
        return Chem.MolToSmiles(m, isomericSmiles=True)
    except Exception:
        try:
            return Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(mol)), isomericSmiles=True)
        except Exception:
            return ""


def mols_from_sdf(path: str):
    suppl = Chem.SDMolSupplier(path, sanitize=True, removeHs=True)
    for mol in suppl:
        if mol is not None:
            yield mol


def read_label_smiles(path: str) -> List[str]:
    """Return canonical SMILES of every molecule in an SDF / SMILES / CSV set."""
    out: List[str] = []
    if is_smiles_source(path):
        for rec in read_smiles_records(path):
            try:
                m = Chem.MolFromSmiles(rec["smiles"])
            except Exception:
                m = None
            if m is not None:
                s = canonical_smiles(m)
                if s:
                    out.append(s)
    else:
        for m in mols_from_sdf(path):
            s = canonical_smiles(m)
            if s:
                out.append(s)
    return out


# ---------------------------------------------------------------------------
# scoring: mode-1 Vina affinity per receptor x docked molecule
# ---------------------------------------------------------------------------
def _affinity_of(res: dict) -> Optional[float]:
    v = res.get("affinity_best")
    if v is None:
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def smiles_of_ligand(meta: dict, cache: Dict[str, str]) -> str:
    """Neutralised canonical SMILES of a prepared ligand.

    Source order: cached -> source SMILES (canonicalised) -> prepared SDF
    (canonicalised).  Returning the raw input SMILES string would break matching
    against the (canonicalised) label set, so the source string is always
    canonicalised here.
    """
    key = meta.get("name", "")
    if key in cache:
        return cache[key]
    cs = ""
    s = meta.get("smiles") or ""
    if s:
        try:
            m = Chem.MolFromSmiles(s)
        except Exception:
            m = None
        if m is not None:
            cs = canonical_smiles(m)
    if not cs:
        prep = meta.get("prep_sdf")
        if prep and os.path.exists(prep):
            try:
                m = Chem.MolFromMolFile(prep, sanitize=True, removeHs=True)
                if m is not None:
                    cs = canonical_smiles(m)
            except Exception:
                cs = ""
    if cs and key:
        cache[key] = cs
    return cs


def docked_scores_by_receptor(
    cfg: RunConfig,
    out_dir: str,
    lig_meta: Dict[str, dict],
) -> Dict[str, Dict[str, dict]]:
    """Map each receptor to ``{canonical_smiles: {"score": float, "source": str}}``.

    ``score`` is the **screening** mode-1 Vina affinity (a single, consistent
    docking stage across the whole library).  Refined shortlist energies are
    deliberately NOT mixed in: refinement only re-docks the best-scoring subset,
    so using "refine if available else screening" would systematically favour the
    best-scoring compounds (usually the actives) and inflate the enrichment
    statistics.
    """
    sub = project_subdirs(out_dir)
    screen_best: Dict[tuple, float] = {}
    for res in batch.iter_results(sub["screening"]):
        a = _affinity_of(res)
        if a is not None:
            screen_best[(res.get("receptor"), res.get("ligand"))] = a

    cache: Dict[str, str] = {}
    per_rec: Dict[str, Dict[str, dict]] = {}
    for (rec, lig), score in screen_best.items():
        if rec not in per_rec:
            per_rec[rec] = {}
        meta = lig_meta.get(lig) or {}
        smi = smiles_of_ligand(meta, cache) if lig else ""
        if not smi:
            continue
        existing = per_rec[rec].get(smi)
        if existing is None or score < existing["score"]:
            per_rec[rec][smi] = {"score": score, "source": "screening"}
    return per_rec


# ---------------------------------------------------------------------------
# ROC / AUC / EF
# ---------------------------------------------------------------------------
def _auc_from_ranks(pos_ranks: Sequence[float], n_pos: int, n_neg: int) -> float:
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    return float((sum(pos_ranks) - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def roc_auc_ef(active_scores: Sequence[float], decoy_scores: Sequence[float]):
    """Return a dict with auc/ef1/ef5 and ROC operating points.

    Inputs are Vina affinities (kcal/mol): more negative = better. Internally the
    score is inverted (``score = -affinity``) so that higher = better, then:

    * AUC is the Mann-Whitney statistic with average ranks for ties
      (``P(score_active > score_decoy)``);
    * ROC points step from the best-scoring molecule to the worst;
    * EF1%/EF5% count actives found in the top 1%/5% of the ranked list divided
      by the expected fraction under random picking.
    """
    n_pos = len(active_scores)
    n_neg = len(decoy_scores)
    auc = float("nan")
    ef = {f"{x:g}": float("nan") for x in EF_PERCENTILES}
    roc_points: List[tuple] = []
    if n_pos and n_neg:
        all_scores = list(active_scores) + list(decoy_scores)
        labels = [1] * n_pos + [0] * n_neg
        # higher = better score
        score = [-float(s) for s in all_scores]
        if _HAVE_SCIPY:
            ranks = rankdata(score)  # ascending; best score gets high rank
            pos_ranks = [r for r, lab in zip(ranks, labels) if lab == 1]
            auc = _auc_from_ranks(pos_ranks, n_pos, n_neg)
        else:  # pragma: no cover - scipy is a declared dependency
            auc = float("nan")
        # best -> worst order (descending score); ties share one threshold
        order = np.argsort(score, kind="mergesort")[::-1]
        sorted_score = np.asarray(score)[order]
        sorted_labels = np.asarray(labels)[order]
        roc_points = [(0.0, 0.0)]
        tp = 0
        fp = 0
        i = 0
        n_tot = n_pos + n_neg
        while i < n_tot:
            s = sorted_score[i]
            j = i
            while j < n_tot and sorted_score[j] == s:
                if sorted_labels[j] == 1:
                    tp += 1
                else:
                    fp += 1
                j += 1
            roc_points.append((fp / n_neg, tp / n_pos))
            i = j
        if roc_points[-1] != (1.0, 1.0):
            roc_points.append((1.0, 1.0))
        # EF at top x%
        for x in EF_PERCENTILES:
            k = int(np.ceil(n_tot * x / 100.0))
            act_in_top = int(sorted_labels[:k].sum()) if k > 0 else 0
            ef[f"{x:g}"] = (act_in_top / n_pos) / (x / 100.0) if n_pos else float("nan")
    return {"auc": auc, "n_active": n_pos, "n_decoy": n_neg,
            "ef": ef, "roc_points": roc_points}


# ---------------------------------------------------------------------------
# report writing
# ---------------------------------------------------------------------------
def _svg_roc(roc_points: Sequence[tuple], auc: float, receptor: str,
             n_act: int, n_dec: int) -> str:
    W, H, PAD = 480, 360, 44
    plot_w, plot_h = W - 2 * PAD, H - 2 * PAD
    pts = list(roc_points)
    if not pts:
        pts = [(0.0, 0.0), (1.0, 1.0)]
    xvals = [p[0] for p in pts]
    yvals = [p[1] for p in pts]

    def _px(x): return PAD + x * plot_w

    def _py(y): return H - PAD - y * plot_h
    path = " ".join(f"{_px(x):.1f},{_py(y):.1f}" for x, y in pts)
    auctxt = f"n/a" if auc != auc else f"{auc:.4f}"
    grid = ""
    for g in (0.0, 0.25, 0.5, 0.75, 1.0):
        grid += (f'<line x1="{_px(g):.1f}" y1="{_py(0):.1f}" '
                 f'x2="{_px(g):.1f}" y2="{_py(1):.1f}" stroke="#ececec" '
                 f'stroke-width="1"/>')
        grid += (f'<line x1="{_px(0):.1f}" y1="{_py(g):.1f}" '
                 f'x2="{_px(1):.1f}" y2="{_py(g):.1f}" stroke="#ececec" '
                 f'stroke-width="1"/>')
    return f"""<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}"
 viewBox="0 0 {W} {H}" font-family="Arial, sans-serif">
<rect width="{W}" height="{H}" fill="white"/>
{grid}
<line x1="{_px(0):.1f}" y1="{_py(0):.1f}" x2="{_px(1):.1f}" y2="{_py(1):.1f}"
 stroke="#bbbbbb" stroke-width="1" stroke-dasharray="4,3"/>
<polyline points="{path}" fill="none" stroke="#1f77b4" stroke-width="2"/>
<text x="{PAD}" y="{H - PAD + 22}" font-size="11" fill="#333">False positive rate</text>
<text x="{PAD - 8}" y="{_py(0):.1f}" text-anchor="end" font-size="11" fill="#333">TPR</text>
<text x="{_px(0):.1f}" y="{_py(1.02):.1f}" font-size="12" fill="#333">ROC</text>
<text x="{W - PAD}" y="{_py(1.02):.1f}" text-anchor="end" font-size="12"
 fill="#1f77b4" font-weight="bold">AUC = {auctxt}</text>
</svg>"""


def build_enrichment_report(cfg: RunConfig, out_dir: str, lig_meta: Dict[str, dict],
                            log=None) -> dict:
    """Run the enrichment workflow when configured; write files; return summary.

    Returns ``{"enabled": False, "note": ...}`` when no enrichment requested, and
    per-receptor summaries under ``{"receptors": {...}}`` otherwise.
    """
    def _log(m):
        if log:
            try:
                log(m)
            except Exception:
                pass

    if not cfg.run_enrichment:
        return {"enabled": False, "note": "run_enrichment is False"}

    missing = [k for k, p in (("actives", cfg.actives_path),
                              ("decoys", cfg.decoys_path)) if not p or not os.path.exists(p)]
    if missing:
        _log("[enrich] SKIPPED: missing files -> " + ", ".join(missing))
        return {"enabled": True, "note": "missing actives/decoys files; enrichment not computed",
                "error": missing}

    actives = set(read_label_smiles(cfg.actives_path))
    decoys = set(read_label_smiles(cfg.decoys_path))
    overlap = actives & decoys
    if overlap:
        # a molecule in both lists cannot be a clean positive and a clean
        # negative; treat it as active and disclose the removal.
        decoys = decoys - actives
        _log(f"[enrich] WARNING: {len(overlap)} molecule(s) present in BOTH actives and "
             f"decoys; counted as active (removed from decoys) and disclosed in the summary.")
    _log(f"[enrich] labels: {len(actives)} actives, {len(decoys)} decoys "
         f"(overlap removed: {len(overlap)})")

    sub = project_subdirs(out_dir)
    os.makedirs(sub["enrichment"], exist_ok=True)
    per_rec = docked_scores_by_receptor(cfg, out_dir, lig_meta)
    rec_summaries = {}
    all_rows_csv = []
    label_attrition = {"active": 0, "decoy": 0}

    for rec, smi_scores in per_rec.items():
        # scores -> labels
        pos_scores = []
        neg_scores = []
        not_scored = []
        for smi, info in smi_scores.items():
            sc = info["score"]
            if smi in actives:
                pos_scores.append(sc)
            elif smi in decoys:
                neg_scores.append(sc)
            else:
                not_scored.append(smi)
        pos_missing = len([a for a in actives if a not in smi_scores])
        neg_missing = len([d for d in decoys if d not in smi_scores])
        label_attrition["active"] += pos_missing
        label_attrition["decoy"] += neg_missing
        res = roc_auc_ef(pos_scores, neg_scores)
        auc = res["auc"]
        if auc != auc:
            _log(f"[enrich] {rec}: insufficient labelled docked data "
                 f"(actives_scored={len(pos_scores)}, decoys_scored={len(neg_scores)})")
        else:
            _log(f"[enrich] {rec}: AUC={auc:.4f}  EF1%={res['ef'].get('1', float('nan')):.2f}  "
                 f"EF5%={res['ef'].get('5', float('nan')):.2f}  "
                 f"(actives_scored={len(pos_scores)}, decoys_scored={len(neg_scores)})")

        svg = _svg_roc(res["roc_points"], auc, rec, len(pos_scores), len(neg_scores))
        svg_path = os.path.join(sub["enrichment"], f"roc_{rec}.svg")
        with open(svg_path, "w", encoding="utf-8") as fh:
            fh.write(svg)

        # ROC points csv
        roc_path = os.path.join(sub["enrichment"], f"roc_points_{rec}.csv")
        with open(roc_path, "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.writer(fh)
            w.writerow(["fpr", "tpr"])
            for fpr, tpr in res["roc_points"]:
                w.writerow([f"{fpr:.6f}", f"{tpr:.6f}"])

        # per-molecule decision csv
        dec_rows = []
        seen_smi = set()
        for smi in list(actives) + list(decoys):
            if smi in seen_smi:
                continue
            seen_smi.add(smi)
            lab = "active" if smi in actives else "decoy"
            info = smi_scores.get(smi)
            dec_rows.append({"receptor": rec, "label": lab,
                             "canonical_smiles": smi,
                             "score_kcal_mol": round(info["score"], 2) if info else "",
                             "source": info["source"] if info else "not-docked",
                             "in_top": ""})
        dec_path = os.path.join(sub["enrichment"], f"molecules_{rec}.csv")
        with open(dec_path, "w", newline="", encoding="utf-8-sig") as fh:
            cols = ["receptor", "label", "canonical_smiles", "score_kcal_mol",
                    "source", "in_top"]
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            w.writerows(dec_rows)

        summary = {"receptor": rec,
                   "n_active_scored": len(pos_scores),
                   "n_decoy_scored": len(neg_scores),
                   "n_active_labels": len(actives),
                   "n_decoy_labels": len(decoys),
                   "auc": round(auc, 4) if auc == auc else None,
                   "ef1": round(res["ef"].get("1", float("nan")), 3),
                   "ef5": round(res["ef"].get("5", float("nan")), 3),
                   "roc_svg": os.path.basename(svg_path),
                   "roc_csv": os.path.basename(roc_path),
                   "molecules_csv": os.path.basename(dec_path),
                   "scoring": "screening mode-1 Vina affinity (single consistent stage)"}
        rec_summaries[rec] = summary
        all_rows_csv.append(summary)
        for row in dec_rows:
            if row["source"] == "not-docked":
                pass

    # aggregate summary file
    summary_path = os.path.join(sub["enrichment"], "enrichment_summary.json")
    utils.write_json({"enabled": True, "receptors": rec_summaries,
                      "label_attrition": label_attrition,
                      "overlap_active_decoy": len(overlap),
                      "actives_file": cfg.actives_path,
                      "decoys_file": cfg.decoys_path,
                      "method": "ROC/AUC(Mann-Whitney, avg-rank ties); EF1%/5%; "
                                "labels matched by canonical isomeric SMILES of the "
                                "NEUTRALISED molecule (pH protomers are matched on the "
                                "charge-free graph); score = screening mode-1 Vina "
                                "affinity (single consistent stage, not mixed with "
                                "refined shortlist energies)"},
                      summary_path)
    if all_rows_csv:
        csv_path = os.path.join(sub["enrichment"], "enrichment_summary.csv")
        with open(csv_path, "w", newline="", encoding="utf-8-sig") as fh:
            cols = list(all_rows_csv[0].keys())
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            w.writerows(all_rows_csv)

    return {"enabled": True, "receptors": rec_summaries,
            "label_attrition": label_attrition,
            "overlap_active_decoy": len(overlap),
            "summary_json": summary_path,
            "auc": next((v["auc"] for v in rec_summaries.values() if v["auc"] is not None), None)}
