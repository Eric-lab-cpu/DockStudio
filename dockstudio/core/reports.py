"""Methods report, quality assessment and README writers (section 12)."""

from __future__ import annotations

import os
from typing import Dict, List, Optional

from .._version import BRAND, COPYRIGHT_CN
from . import batch, envinfo, qc as qc_mod, utils
from .models import project_subdirs


def _fmt_rows_csv(path: str) -> str:
    import csv
    if not os.path.exists(path):
        return "_missing_"
    with open(path, encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        return "_empty_"
    return "\n".join("| " + " | ".join(str(r.get(c, "")) for c in rows[0]) + " |" for r in rows[:20])


def directory_tree(root: str, maxdepth: int = 4) -> str:
    lines = []
    root = os.path.abspath(root)
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
        depth = dirpath[len(root):].count(os.sep)
        if depth > maxdepth:
            dirnames[:] = []
            continue
        indent = "  " * depth
        lines.append(f"{indent}{os.path.basename(dirpath) or root}/")
        for f in sorted(filenames):
            if f.endswith((".pyc", ".tmp")):
                continue
            try:
                sz = os.path.getsize(os.path.join(dirpath, f))
            except OSError:
                sz = 0
            lines.append(f"{indent}  {f}  ({sz} B)")
    return "\n".join(lines)


def write_methods_report(cfg, out_dir: str, env: dict, inventory_summary: dict,
                         boxes: dict, prep_summary: dict) -> str:
    sub = project_subdirs(out_dir)
    results = batch.collect_results(sub["screening"])
    refine_res = batch.collect_results(sub["refine"])
    path = os.path.join(sub["reports"], "01_methods_report.md")
    L = []
    L.append("# 方法学与参数记录(Methods & Parameters)\n")
    L.append(f"- 生成时间:{utils.now_str()}")
    L.append(f"- 项目标题:{cfg.title or '(未命名)'}")
    L.append(f"- 输出目录:`{cfg.out_dir}`")
    L.append("\n## 1. 环境\n")
    L.append("| 工具 | 版本 |")
    L.append("| --- | --- |")
    for t, v in env.get("tools", {}).items():
        L.append(f"| {t} | {v} |")
    L.append("| python | " + env.get("python", {}).get("python", "?") + " |")
    L.append("\n## 2. 输入盘点\n")
    L.append("### 受体\n")
    L.append(_fmt_rows_csv(os.path.join(out_dir, "receptor_inventory.csv")))
    L.append("\n### 配体\n")
    L.append(_fmt_rows_csv(os.path.join(out_dir, "ligand_list.csv")))
    L.append("\n## 3. 搜索盒子\n")
    for rec, b in boxes.items():
        L.append(f"- `{rec}`: method={b.get('method')} center={b.get('center')} "
                 f"size={b.get('size')} basis={b.get('basis')}")
    L.append("\n## 4. 结构准备\n")
    for rec, p in (prep_summary or {}).items():
        L.append(f"- 受体 `{rec}`: atoms_in_pdbqt={p.get('n_atoms_pdbqt')} "
                 f"deleted_bad_residues={p.get('deleted_bad_residues')}")
        for d in (p.get("deleted_bad") or []):
            L.append(f"  - 删除不完整残基(距位点>10 A): {d}")
        for w in (p.get("clean_warnings") or []):
            L.append(f"  - 清洗提示(如实记录): {w}")
    L.append("\n## 5. 对接参数\n")
    L.append(f"- Vina scoring: `vina`; screening exhaustiveness={cfg.exhaustiveness}; "
             f"refinement exhaustiveness={cfg.refine_exhaustiveness}; n_poses={cfg.n_poses}; "
             f"energy_range={cfg.energy_range}; Top-K={cfg.top_k}")
    L.append(f"- 配体质子化: 文档化简单规则,目标 pH={cfg.ph}(羧酸 pH≥5 去质子化;"
             f"碱性胺 pH≤8.5 质子化)。非 pKa 预测,见质量评估报告局限。")
    L.append("- 随机种子: 每对 `crc32('<receptor>|<ligand>') & 0x7fffffff`(写入各 result.json)")
    L.append(f"- 批处理: 断点续跑(done.flag), 幂等续跑; cpu={cfg.cpu}")
    L.append("\n## 6. 完成的对接任务\n")
    L.append(f"- 初筛完成对子数:{len(results)}(含 pre-exist 幂等跳过)")
    L.append(f"- 精修完成对子数:{len(refine_res)}")
    L.append("\n### 初筛 mode-1 结合能矩阵\n")
    L.append(_fmt_rows_csv(os.path.join(out_dir, "results", "screening_matrix.csv")))
    L.append("\n### 最终 Top-K\n")
    L.append(_fmt_rows_csv(os.path.join(out_dir, "results",
                                        f"final_top{cfg.top_k}.csv")))
    L.append("\n## 7. 回贴验证(Self-docking)\n")
    L.append(_fmt_rows_csv(os.path.join(out_dir, "results", "selfdock_summary.csv")))
    L.append("\n## 8. 相互作用分析\n")
    L.append(_fmt_rows_csv(os.path.join(out_dir, "results", "plip", "plip_interactions_all.csv")))
    L.append("\n## 9. 已执行/未执行的说明\n")
    L.append("本报告只记录本目录中真实生成的文件与参数;任何未执行的步骤不在此列出。")
    L.append("\n## 10. 平台约束\n")
    L.append(f"- 平台:{env.get('python', {}).get('platform', '?')}")
    L.append("- 说明: Windows 下若使用中文路径,PyMOL .pml 脚本已设计为纯 ASCII 相对路径,"
             "规避 GBK 解码问题。")
    L.append(f"\n---\n\n{BRAND} · {COPYRIGHT_CN}\n")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))
    return path


def _selfdock_pass_details(out_dir: str) -> List[str]:
    """Read selfdock_summary.json (if present) into human lines."""
    sd_path = os.path.join(out_dir, "results", "selfdock_summary.json")
    if not os.path.exists(sd_path):
        return []
    try:
        rows = utils.read_json(sd_path).get("selfdock", [])
    except Exception:
        return []
    out = []
    for r in rows:
        rmsd = r.get("mode1_rmsd")
        if rmsd is None:
            continue
        try:
            rmsd = float(rmsd)
            ok = rmsd <= 2.0
        except (TypeError, ValueError):
            ok = None
        state = "PASS" if ok else ("FAIL" if ok is False else "n/a")
        out.append(f"receptor {r.get('receptor')} ({r.get('cocrystal_ligand', '')}): "
                   f"self-dock mode-1 RMSD = {rmsd} A -> {state}")
    return out


def write_quality_report(cfg, out_dir: str, checks: List[dict],
                         selfdock_rows: Optional[list] = None) -> str:
    sub = project_subdirs(out_dir)
    path = os.path.join(sub["reports"], "02_quality_assessment.md")
    L = []
    L.append("# 准确性与质量评估\n")
    L.append("\n## 1. 自动一致性检查\n")
    for c in checks:
        L.append(f"- [{'x' if c['ok'] else ' '}] {c['check']}: {c['detail']}")
    sd_lines = _selfdock_pass_details(out_dir)
    if sd_lines:
        L.append("\n## 2. 共晶配体回贴(self-dock)一致性\n")
        L.extend("- " + x for x in sd_lines)
        L.append("\n> 阈值:mode-1 RMSD <= 2.0 A 视为通过。RMSD 为同一笛卡尔坐标系下"
                 "重原子匹配(无重新叠加),见 `03_accuracy_assessment.md`。")
    L.append("\n## 3. 已发现并修正的错误\n")
    L.append("- 以实际工作日志为准;本报告列出已知记录缺口,不做虚构。")
    L.append("\n## 4. 记录缺口(如存在)\n")
    gaps = []
    if not os.path.exists(os.path.join(out_dir, "results", "plip", "plip_interactions_all.csv")):
        gaps.append("PLIP 交互表缺失(可能 PLIP 不可用或阶段被关闭)。")
    if not checks:
        gaps.append("未执行自动一致性检查。")
    L.append("\n".join("- " + g for g in gaps) if gaps else "- 无。")
    L.append("\n## 5. 科学局限\n")
    L.append("- Vina 打分函数的已知噪声(~0.5-1 kcal/mol);名次差 <0.1 kcal/mol 的命中不应过度解读。")
    L.append(f"- 质子化为文档化简单规则(目标 pH {cfg.ph};羧酸 pH≥5 去质子、碱性胺 "
             f"pH≤8.5 质子化),未使用 pKa 预测工具;手性/互变异构/关键可电离基团需人工复核。")
    L.append("- 无共晶配体受体的对接为探索性,回贴验证不适用,假阳性风险高。")
    L.append("- 图片未经人工目视终审;发布前应检查 PyMOL 渲染与结合模式合理性。")
    L.append("\n## 6. 不宜过度解读的结果\n")
    L.append("- 被标记为 `screening-fallback` 的最终 Top-K 条目(精修失败时);")
    L.append("- 无共晶配体受体(exploratory/blind)的全部结果;")
    L.append("- 回贴 mode-1 RMSD > 2.0 A 受体上的所有排序(打分/姿态未被实验几何验证)。")
    L.append(f"\n---\n\n{BRAND} · {COPYRIGHT_CN}")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    return path


def write_receptor_protocol(cfg, out_dir: str, selfdock_rows: Optional[list] = None) -> str:
    """Write receptor_protocol.csv (protocol section 4.3): per receptor chain,
    cocrystal ligand, box, exhaustiveness and self-dock RMSD."""
    import csv
    sd = {}
    for row in (selfdock_rows or []):
        sd[row.get("receptor")] = row
    path = os.path.join(out_dir, "receptor_protocol.csv")
    cols = ["receptor", "chains", "cocrystal_ligand", "box_method", "box_center",
            "box_size", "exhaustiveness", "refine_exhaustiveness",
            "selfdock_mode1_rmsd_A", "selfdock_pass"]
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for rec in cfg.receptors:
            nm = rec["name"]
            box = cfg.boxes.get(nm, {})
            s = sd.get(nm, {})
            w.writerow({
                "receptor": nm,
                "chains": ";".join(rec.get("chains") or []),
                "cocrystal_ligand": box.get("ligand_ref", ""),
                "box_method": box.get("method", ""),
                "box_center": ";".join(str(round(x, 2)) for x in box.get("center", [])),
                "box_size": ";".join(str(round(x, 2)) for x in box.get("size", [])),
                "exhaustiveness": cfg.exhaustiveness,
                "refine_exhaustiveness": cfg.refine_exhaustiveness,
                "selfdock_mode1_rmsd_A": s.get("mode1_rmsd", ""),
                "selfdock_pass": s.get("pass", ""),
            })
    return path


def write_readme(out_dir: str, cfg, env: dict) -> str:
    path = os.path.join(out_dir, "README.md")
    L = [
        "# 运行结果说明",
        "",
        f"- 标题:{cfg.title or '(未命名)'}",
        f"- 生成时间:{utils.now_str()}",
        "",
        "## 目录结构",
        "```",
        directory_tree(out_dir, maxdepth=3),
        "```",
        "",
        "## 主要产物",
        "- `receptor_inventory.csv` / `ligand_list.csv`:输入盘点",
        "- `docking/*/result.json`:每对初筛结果(含种子、盒子、affinities)",
        "- `refine/*/result.json`:精修结果",
        "- `results/screening_matrix.csv`:受体 x 配体初筛结合能矩阵",
        "- " + "`results/final_top" + str(cfg.top_k) + ".csv` + `_final_top" + str(cfg.top_k) + ".json`:最终 Top-K",
        "- `results/selfdock_summary.csv`:共晶配体回贴 RMSD",
        "- `reports/accuracy_assessment.csv` + `reports/03_accuracy_assessment.md`:"
        "对接准确性评估报告(如启用)",
        "- `results/plip/plip_interactions_all.csv`:相互作用表",
        "- `results/3D_poses/*.png`,`results/composite/*.png`,`results/topK_panels/*.png`:图件",
        "- `results/pymol_data/*/scene.pml` + `complex.pdb`:可编辑 ASCII PyMOL 源",
        "- `results/pse/*.pse`:PyMOL 会话(如启用)",
        "- `reports/`:方法与质量报告",
        "",
        "详细方法见 `reports/01_methods_report.md`;质量评估见 `reports/02_quality_assessment.md`;",
        "对接准确性评估(启用时)见 `reports/03_accuracy_assessment.md`。",
        "",
    ]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))
    return path
