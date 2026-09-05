"""Docking-accuracy assessment report (self-docking redocking validation).

v1.1 addition. For every receptor that carries a co-crystal ligand, DockStudio
already performs a "self-docking" run during :mod:`dockstudio.core.selfdock`
(redock the native ligand into its own crystal pocket with the same box used
for screening).  This module turns those raw RMSD numbers into a formal,
publication-oriented **docking-accuracy assessment report**:

    reports/03_accuracy_assessment.md
    reports/accuracy_assessment.csv      (one row per receptor)
    reports/accuracy_poses.csv           (one row per pose, when available)

Honesty rules followed here (do not fabricate):
  * only numbers that exist on disk under ``out_dir`` are used;
  * receptors without a co-crystal ligand are reported as *not assessable*
    (blind / exploratory docking cannot be validated by redocking);
  * the report always restates method parameters actually used.

Definition of "docking accuracy" used by this report
------------------------------------------------------
A pose is a spatial prediction of how the ligand binds.  The *redocking RMSD*
compares the docked pose against the crystal pose of the same ligand in the
*identical* Cartesian frame (no re-superposition), which is the standard
self-docking protocol definition.  Interpretation bands (documented, not
over-claimed):

    RMSD <= 1.0 A                 high accuracy  (near-crystallographic)
    1.0 A < RMSD <= 2.0 A         acceptable
    RMSD > 2.0 A                  poor / FAIL    (pose not reproduced)

These bands are literature conventions for self-docking redocking experiments
and are only meaningful *for targets that have a co-crystal reference*.
"""

from __future__ import annotations

import csv
import os
from typing import Dict, List, Optional

from .._version import (
    ACCURACY_REPORT_BASELINE_MEDIUM,
    ACCURACY_REPORT_BASELINE_STRICT,
    BRAND,
    COPYRIGHT_CN,
    SELFDOCK_PASS_RMSD,
)
from . import models, utils
from .models import RunConfig, project_subdirs


def _verdict(rmsd: float) -> str:
    if rmsd <= ACCURACY_REPORT_BASELINE_STRICT:
        return "high"
    if rmsd <= ACCURACY_REPORT_BASELINE_MEDIUM:
        return "acceptable"
    return "poor"


def _collect_receptor_rows(cfg: RunConfig, out_dir: str) -> Dict[str, dict]:
    """Read selfdock_summary.json -> per-receptor evaluation dicts."""
    summary_path = os.path.join(out_dir, "results", "selfdock_summary.json")
    if not os.path.exists(summary_path):
        return {}
    try:
        data = utils.read_json(summary_path).get("selfdock", [])
    except Exception:
        return {}
    rows: Dict[str, dict] = {}
    for r in data:
        rows.setdefault(r.get("receptor"), {}).update(r)
    return rows


def _rmsd_of(r: dict, key: str) -> Optional[float]:
    v = r.get(key)
    if v is None or v == "":
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def build_accuracy_report(cfg: RunConfig, out_dir: str, log=None) -> dict:
    """Build the docking-accuracy assessment report files.

    Returns a summary dict for the pipeline / GUI (never raises on missing
    data -- it writes an honest "not assessable" report instead).
    """
    sub = project_subdirs(out_dir)
    reports_dir = sub["reports"]
    os.makedirs(reports_dir, exist_ok=True)
    selfdock_rows = _collect_receptor_rows(cfg, out_dir)

    def _log(m: str):
        if log:
            try:
                log(m)
            except Exception:
                pass

    receptor_rows = []
    pose_rows = []
    n_assessable = 0
    n_pass_2a = 0
    n_pass_1a = 0
    rmsd_values = []

    for rec in cfg.receptors:
        name = rec["name"]
        box = cfg.boxes.get(name, {})
        lig_ref = box.get("ligand_ref")
        sd = selfdock_rows.get(name, {})
        mode1 = _rmsd_of(sd, "mode1_rmsd")
        min_rmsd = _rmsd_of(sd, "min_rmsd")
        # energy bookkeeping from selfdock evaluation (stored in selfdock_summary):
        energy_gap = _rmsd_of(sd, "energy_delta_min_rmsd")
        for p in sd.get("poses", []) or []:
            pose_rows.append({
                "receptor": name,
                "pose": p.get("pose", ""),
                "affinity_kcal_mol": p.get("affinity", ""),
                "rmsd_A": p.get("rmsd", ""),
            })

        if mode1 is None:
            reason = (
                "no co-crystal reference in docking chains (blind/exploratory docking; "
                "redocking accuracy not assessable)"
                if not lig_ref else "self-docking did not produce a mode-1 RMSD "
                                     "(check selfdock_summary.json / run log)"
            )
            receptor_rows.append({
                "receptor": name,
                "cocrystal_ligand": lig_ref or "",
                "box_method": box.get("method", ""),
                "exhaustiveness": cfg.refine_exhaustiveness if cfg.run_refine else cfg.exhaustiveness,
                "mode1_rmsd_A": "",
                "min_rmsd_A": "",
                "min_rmsd_pose": "",
                "energy_delta_mode1_min_A": "",
                "pass_2A": "",
                "verdict": "not_assessable",
                "note": reason,
            })
            continue

        n_assessable += 1
        rmsd_values.append(mode1)
        pass_2a = mode1 <= SELFDOCK_PASS_RMSD
        if pass_2a:
            n_pass_2a += 1
        if mode1 <= ACCURACY_REPORT_BASELINE_STRICT:
            n_pass_1a += 1
        min_pose = sd.get("min_rmsd_pose")
        # gap = energy(lowest-RMSD pose) - energy(mode 1); a large positive gap
        # with a low RMSD pose elsewhere means the scorer did not rank the best
        # geometry first (scoring deficiency), a small gap means it did.
        gap = energy_gap
        receptor_rows.append({
            "receptor": name,
            "cocrystal_ligand": lig_ref or sd.get("cocrystal_resname", ""),
            "box_method": box.get("method", ""),
            "exhaustiveness": cfg.refine_exhaustiveness if cfg.run_refine else cfg.exhaustiveness,
            "mode1_rmsd_A": round(mode1, 3),
            "min_rmsd_A": round(min_rmsd, 3) if min_rmsd is not None else "",
            "min_rmsd_pose": min_pose if min_pose is not None else "",
            "energy_delta_mode1_min_A": round(gap, 2) if gap is not None else "",
            "pass_2A": "yes" if pass_2a else "no",
            "verdict": _verdict(mode1),
            "note": (
                f"lowest-energy pose (mode 1) RMSD = {mode1:.2f} A; "
                f"best-RMSD pose #{min_pose} = "
                f"{min_rmsd:.2f} A" if min_rmsd is not None
                else f"lowest-energy pose (mode 1) RMSD = {mode1:.2f} A"
            ),
        })

    # ---- write per-receptor CSV -------------------------------------------------
    rec_cols = ["receptor", "cocrystal_ligand", "box_method", "exhaustiveness",
                "mode1_rmsd_A", "min_rmsd_A", "min_rmsd_pose",
                "energy_delta_mode1_min_A", "pass_2A", "verdict", "note"]
    csv_path = os.path.join(reports_dir, "accuracy_assessment.csv")
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=rec_cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(receptor_rows)

    pose_path = os.path.join(reports_dir, "accuracy_poses.csv")
    if pose_rows:
        with open(pose_path, "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=["receptor", "pose", "affinity_kcal_mol", "rmsd_A"])
            w.writeheader()
            w.writerows(pose_rows)

    # ---- write markdown report --------------------------------------------------
    L = []
    L.append("# 分子对接准确性评估报告(Docking Accuracy Assessment)")
    L.append("")
    L.append(f"- 生成时间:{utils.now_str()}")
    L.append(f"- 项目标题:{cfg.title or '(未命名)'}")
    L.append(f"- 评估方式:**共晶配体回贴(self-docking redocking)**,同一笛卡尔坐标系下"
             f"对接姿态 vs 晶体姿态的重原子 RMSD(无重新叠加)。")
    L.append(f"- 判定阈值:<={ACCURACY_REPORT_BASELINE_STRICT:.1f} A = 高精度;"
             f"<= {ACCURACY_REPORT_BASELINE_MEDIUM:.1f} A = 可接受;"
             f"> {SELFDOCK_PASS_RMSD:.1f} A = FAIL / 姿态未重现。")
    L.append("")
    L.append("## 1. 总览")
    L.append("")
    L.append(f"- 受体总数:{len(cfg.receptors)}")
    L.append(f"- 可评估(含共晶配体):{n_assessable}")
    L.append(f"- 通过(<= {SELFDOCK_PASS_RMSD:.1f} A):{n_pass_2a}")
    L.append(f"- 高精度(<= {ACCURACY_REPORT_BASELINE_STRICT:.1f} A):{n_pass_1a}")
    if rmsd_values:
        mean = sum(rmsd_values) / len(rmsd_values)
        med = sorted(rmsd_values)[len(rmsd_values) // 2]
        L.append(f"- mode-1 RMSD 均值:{mean:.2f} A;中位数:{med:.2f} A")
    else:
        L.append("- mode-1 RMSD:无可评估受体(无共晶配体或自对接未运行)。")
    L.append("")
    L.append("## 2. 逐受体评估")
    L.append("")
    if receptor_rows:
        L.append("| 受体 | 共晶配体 | 盒子方法 | exhaustiveness | mode-1 RMSD (A) | "
                 "min RMSD (A) | 判定 |")
        L.append("| --- | --- | --- | --- | --- | --- | --- |")
        for r in receptor_rows:
            ver = r["verdict"]
            if ver == "high":
                label = "高精度"
            elif ver == "acceptable":
                label = "可接受"
            elif ver == "poor":
                label = "FAIL"
            else:
                label = "不可评估"
            L.append(f"| {r['receptor']} | {r['cocrystal_ligand'] or '—'} | "
                     f"{r['box_method'] or '—'} | {r['exhaustiveness']} | "
                     f"{r['mode1_rmsd_A'] if r['mode1_rmsd_A'] != '' else '—'} | "
                     f"{r['min_rmsd_A'] if r['min_rmsd_A'] != '' else '—'} | {label} |")
    else:
        L.append("_无逐受体记录。_")
    L.append("")
    L.append("## 3. 姿态级数据")
    L.append("")
    if pose_rows:
        L.append("每个可评估受体的姿态(RMSD vs 晶体,按 pose 序号):详见 "
                 "`accuracy_poses.csv`。mode-1 为 Vina 打分最低能姿态。")
        L.append("")
        L.append("**解释:** 若最低能姿态(mode-1)RMSD 小(<=2 A),说明打分函数能把"
                 "接近晶体的姿态排在第一位(打分与几何一致);若 mode-1 RMSD 大但存在"
                 "低 RMSD 姿态且能量差很小,说明几何搜索覆盖了正确姿态但打分排序不够"
                 "锐利——属于打分而非采样问题。")
    else:
        L.append("_无姿态级数据。_")
    L.append("")
    L.append("## 4. 局限与诚实声明")
    L.append("")
    L.append("- 本报告**只反映自对接回贴精度**,不评估:(i) 无共晶配体靶点的对接"
             "(探索性,不可验证);(ii) 打分/富集精度(需要活性/非活性实验数据);"
             "(iii) 蛋白柔性/诱导契合。")
    L.append("- RMSD 使用同元素原子贪心最近邻匹配,未做对称性修正;对称配体(如苯环"
             "翻转、甲基旋转)的 RMSD 可能被高估。")
    L.append(f"- 质子化采用文档化简单规则(目标 pH {cfg.ph};羧酸 pH≥5 去质子化、"
             f"碱性胺 pH≤8.5 质子化),未使用 pKa 预测;手性/互变异构需人工复核后再下结论。")
    L.append(f"- 运行环境与完整参数见 `01_methods_report.md`;自动一致性检查见 "
             f"`02_quality_assessment.md`。")
    L.append("")
    L.append(f"---")
    L.append("")
    L.append(f"{BRAND} · {COPYRIGHT_CN}")
    md_path = os.path.join(reports_dir, "03_accuracy_assessment.md")
    with open(md_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")

    summary = {
        "report_md": md_path,
        "report_csv": csv_path,
        "pose_csv": pose_path if pose_rows else "",
        "n_receptors": len(cfg.receptors),
        "n_assessable": n_assessable,
        "n_pass_2A": n_pass_2a,
        "n_pass_strict_1A": n_pass_1a,
        "mean_mode1_rmsd_A": round(sum(rmsd_values) / len(rmsd_values), 3) if rmsd_values else None,
    }
    _log(f"[accuracy] report written: {os.path.relpath(md_path, out_dir)} "
         f"(assessable={n_assessable}, pass2A={n_pass_2a})")
    return summary
