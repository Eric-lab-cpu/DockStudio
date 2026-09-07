"""Interactive HTML exports (v2.0, feature #17).

Two sibling products of the markdown reports are produced:

* a per-complex **interactive 3D viewer** (``results/html_viewers/<complex>.html``)
  that shows the receptor cartoon, the docked ligand, the pocket residues
  (auto-detected, receptor atoms within 6 A of the ligand) and the PLIP
  interaction annotation as text;
* an **interactive HTML total report** (``reports/html/index.html``) with a
  table of contents, the Top-K table, the docking-accuracy assessment, the
  enrichment (ROC/AUC/EF) section, links to every viewer and the environment /
  methods summary.

3D rendering uses **3Dmol.js**. The generated pages try to embed a locally
bundled copy of ``3Dmol-min.js`` (from ``dockstudio/resources``) so they are
fully offline; when that file is not present the page falls back to a CDN
``<script src>`` and the footer says so honestly.

Important honesty rule: the 3D view is an *auxiliary interactive inspection*
tool, not a scientific computation. Geometrically guessed "interaction lines"
are deliberately NOT drawn; the interaction annotation shown is the real PLIP
text (residue / type / distance). No GUI/tkinter dependency.
"""

from __future__ import annotations

import base64
import csv
import json
import os
from typing import Dict, List, Optional

from .._version import __version__ as VERSION
from . import batch, utils
from .models import project_subdirs

_3DMOL_CDN = "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.1/3Dmol-min.js"
_LIG_RESNAME = "LIG"


def _local_3dmol_snippet() -> str:
    """Embed bundled 3Dmol.js when present; else return the CDN tag."""
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # dockstudio/
    for cand in (os.path.join(here, "resources", "3Dmol-min.js"),
                 os.path.join(here, "..", "assets", "3Dmol-min.js"),
                 os.path.join(os.path.dirname(here), "assets", "3Dmol-min.js")):
        if os.path.exists(cand):
            try:
                with open(cand, encoding="utf-8", errors="replace") as fh:
                    text = fh.read()
                if text.strip():
                    return f"<script>\n{text}\n</script>"
            except Exception:
                break
    return (f'<script src="{_3DMOL_CDN}"></script>\n'
            f'<script>window._3DmolFromCDN=true;</script>')


def _pocket_residues(pdb_text: str, cutoff: float = 6.0) -> List[dict]:
    """Receptor residues (chain, resi) with any atom within ``cutoff`` A of LIG."""
    atoms = []
    lig = []
    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        except ValueError:
            continue
        resn = line[17:20].strip()
        chain = line[21].strip()
        try:
            resi = int(line[22:26].strip())
        except ValueError:
            continue
        rec = {"x": x, "y": y, "z": z, "resn": resn, "chain": chain, "resi": resi}
        atoms.append(rec)
        if resn == _LIG_RESNAME:
            lig.append(rec)
    if not lig:
        return []
    out = []
    for a in atoms:
        if a["resn"] == _LIG_RESNAME:
            continue
        hit = False
        for l in lig:
            d2 = (a["x"] - l["x"]) ** 2 + (a["y"] - l["y"]) ** 2 + (a["z"] - l["z"]) ** 2
            if d2 <= cutoff * cutoff:
                hit = True
                break
        if hit and {"chain": a["chain"], "resi": a["resi"]} not in out:
            out.append({"chain": a["chain"], "resi": a["resi"]})
    out.sort(key=lambda r: (r["chain"], r["resi"]))
    return out


def _interaction_lines(name: str, plip_summary: Optional[dict]) -> List[str]:
    """Real PLIP text annotations for one complex (never fabricated)."""
    if not plip_summary:
        return []
    entry = plip_summary.get(name) or {}
    if entry.get("error"):
        return [f"PLIP: {entry['error']}"]
    lines = []
    for it in entry.get("interactions") or []:
        t = it.get("type", "")
        res = it.get("residue", "")
        dist = it.get("distance_A", "")
        lines.append(f"{t} - {res}" + (f"  ({dist} A)" if dist else ""))
    return lines


# ---------------------------------------------------------------------------
# per-complex 3D viewer
# ---------------------------------------------------------------------------
def _viewer_html(name: str, pdb_text: str, pocket: List[dict],
                 interactions: List[str], used_cdn: bool,
                 mol_snippet: str = "") -> str:
    b64 = base64.b64encode(pdb_text.encode("utf-8")).decode("ascii")
    pocket_json = json.dumps(pocket)
    inter_json = json.dumps(interactions or [])
    offline_note = ("<p style='color:#a66;font-size:12px'>首次打开需联网从 CDN 加载 "
                    "3Dmol.js(本构建未内置本地副本);完整离线请将 3Dmol-min.js 放至 "
                    "software 的 dockstudio/resources 目录后重新生成。</p>" if used_cdn
                    else "<p style='color:#2a7;font-size:12px'>3Dmol.js 已内置,可离线查看。</p>")
    n_pocket = len(pocket)
    template = """<!DOCTYPE html>
<html lang="zh"><head><meta charset="utf-8">
<title>DockStudio 3D 查看器 - @@NAME@@</title>
<style>
body{font-family:'Segoe UI','Microsoft YaHei',sans-serif;margin:0;background:#f5f7fa;}
#viewer{position:absolute;top:0;bottom:0;left:0;right:0;}
.toolbar{position:absolute;top:8px;right:12px;z-index:5;}
.toolbar button{margin-left:6px;padding:4px 10px;border:1px solid #bcd;border-radius:4px;
background:#fff;cursor:pointer;font-size:12px;}
.info{position:absolute;left:10px;bottom:10px;z-index:5;background:rgba(255,255,255,.94);
border:1px solid #dde;border-radius:6px;padding:8px 12px;font-size:12px;max-width:360px;
box-shadow:0 2px 8px rgba(0,0,0,.06);}
h1{font-size:14px;margin:0 0 4px;} ul{margin:4px 0 0;padding-left:18px;}
</style></head><body>
@@3DMOL@@
<div id="viewer"></div>
<div class="toolbar">
 <button onclick="viewer.spin(true)">旋转</button>
 <button onclick="viewer.spin(false)">停止</button>
 <button onclick="viewer.zoomTo()">复位</button>
 <button onclick="togglePocket()">口袋残基开关</button>
</div>
<div class="info">
 <h1>@@NAME@@ · DockStudio v@@VERSION@@</h1>
 <div>口袋残基数:@@NPOCKET@@ | 配体原子按橙色显示,口袋残基按青色 stick 显示。</div>
 <div id="inter"><b>PLIP 相互作用注释(文本):</b><ul id="interList"></ul></div>
 @@OFFLINE@@
 <div style="color:#888;margin-top:4px">说明:查看器仅用于交互式查看,非科学计算依据;不绘制几何猜想的相互作用连线。</div>
</div>
<script>
var pdbTxt = decodeURIComponent(escape(atob("@@PDB@@")));
var pocket = @@POCKET@@;
var inter = @@INTER@@;
var viewer = null;
try {
  viewer = $3Dmol.createViewer("viewer",{backgroundColor:"white"});
  viewer.addModel(pdbTxt,"pdb");
  viewer.setStyle({},{cartoon:{color:"spectrum"}});
  viewer.setStyle({resn:"@@LIG@@"},{stick:{colorscheme:"orangeCarbon",radius:0.18},
                          sphere:{scale:0.22}});
  viewer.setStyle({resi:pocket.map(function(p){return p.resi;})},
                  {stick:{colorscheme:"cyanCarbon",radius:0.14}});
  viewer.zoomTo();
  viewer.render();
} catch(e) {
  document.getElementById("viewer").innerHTML=
    "<p style='padding:20px'>3D 查看器初始化失败("+e+")。请联网后刷新,或将内置 3Dmol.js 补齐。</p>";
}
var _pocketOn=true;
function togglePocket() {
  _pocketOn=!_pocketOn;
  try {
    viewer.setStyle({resi:pocket.map(function(p){return p.resi;})},{});
    viewer.setStyle({resn:"@@LIG@@"},{stick:{colorscheme:"orangeCarbon",radius:0.18}});
    if(_pocketOn) viewer.setStyle({resi:pocket.map(function(p){return p.resi;})},
                      {stick:{colorscheme:"cyanCarbon",radius:0.14}});
    viewer.render();
  } catch(e) {}
}
var il=document.getElementById("interList");
if(!inter.length){il.innerHTML="<li>无(PLIP 未运行或无相互作用记录)</li>";}
else {inter.forEach(function(s){var li=document.createElement("li");li.textContent=s;il.appendChild(li);});}
</script>
</body></html>
"""
    return (template.replace("@@NAME@@", name)
                    .replace("@@VERSION@@", VERSION)
                    .replace("@@NPOCKET@@", str(n_pocket))
                    .replace("@@PDB@@", b64)
                    .replace("@@POCKET@@", pocket_json)
                    .replace("@@INTER@@", inter_json)
                    .replace("@@OFFLINE@@", offline_note)
                    .replace("@@3DMOL@@", mol_snippet)
                    .replace("@@LIG@@", _LIG_RESNAME))


def _write_viewer(name: str, sub: dict, plip_summary: Optional[dict],
                  used_cdn: bool, mol_snippet: str = "") -> Optional[str]:
    cands = []
    cp = os.path.join(sub["results"], "complexes", f"{name}.pdb")
    pp = os.path.join(sub["results"], "pymol_data", name, "complex.pdb")
    if os.path.exists(cp):
        cands.append(cp)
    if os.path.exists(pp) and pp not in cands:
        cands.append(pp)
    if not cands:
        return None
    with open(cands[0], encoding="utf-8", errors="replace") as fh:
        pdb_text = fh.read()
    pocket = _pocket_residues(pdb_text)
    inter = _interaction_lines(name, plip_summary)
    html = _viewer_html(name, pdb_text, pocket, inter, used_cdn, mol_snippet)
    os.makedirs(sub["html_viewers"], exist_ok=True)
    out = os.path.join(sub["html_viewers"], f"{name}.html")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(html)
    return out


# ---------------------------------------------------------------------------
# total HTML report
# ---------------------------------------------------------------------------
def _top_rows(out_dir: str, top_k: int) -> List[dict]:
    sub = project_subdirs(out_dir)
    path = os.path.join(sub["results"], f"final_top{top_k}.csv")
    if not os.path.exists(path):
        return []
    rows = []
    with open(path, encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh))
    return rows


def _html_table(rows: List[dict]) -> str:
    if not rows:
        return "<p>_无数据_</p>"
    cols = list(rows[0].keys())
    head = "".join(f"<th>{c}</th>" for c in cols)
    body = ""
    for r in rows:
        body += "<tr>" + "".join(f"<td>{r.get(c, '')}</td>" for c in cols) + "</tr>"
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def _env_table(env: dict) -> str:
    rows = []
    for k, v in (env.get("tools") or {}).items():
        rows.append({"tool": k, "version": v})
    for k in ("python", "platform"):
        pass
    py = env.get("python") or {}
    rows.append({"tool": "python", "version": py.get("python", "?")})
    return _html_table(rows)


def _roc_svg_inline(path: str) -> str:
    if not os.path.exists(path):
        return "<p>ROC 图缺失</p>"
    try:
        with open(path, encoding="utf-8") as fh:
            return fh.read()
    except Exception:
        return "<p>ROC 图读取失败</p>"


def build_html_report(cfg, out_dir: str, log=None) -> dict:
    """Generate per-complex viewers + the HTML total report (v2.0)."""
    def _log(m):
        if log:
            try:
                log(m)
            except Exception:
                pass

    sub = project_subdirs(out_dir)
    viewers_dir = sub["html_viewers"]
    report_dir = sub["html_report"]
    os.makedirs(viewers_dir, exist_ok=True)
    os.makedirs(report_dir, exist_ok=True)
    snippet = _local_3dmol_snippet()
    used_cdn = "window._3DmolFromCDN" in snippet

    # data sources (real files only)
    top_rows = _top_rows(out_dir, max(1, int(getattr(cfg, "top_k", 5) or 5)))
    env = {}
    env_path = os.path.join(out_dir, "environment.json")
    if os.path.exists(env_path):
        try:
            env = utils.read_json(env_path)
        except Exception:
            env = {}
    plip_summary = None
    psp = os.path.join(sub["results"], "plip", "plip_summary.json")
    if os.path.exists(psp):
        try:
            plip_summary = utils.read_json(psp)
        except Exception:
            plip_summary = None

    viewer_links = []
    for row in top_rows:
        name = f"{row.get('receptor')}__{row.get('ligand')}"
        html_path = _write_viewer(name, sub, plip_summary, used_cdn, snippet)
        if html_path:
            rel = os.path.relpath(html_path, report_dir).replace(os.sep, "/")
            viewer_links.append(f'<li><a href="{rel}" target="_blank">{name}</a></li>')
    _log(f"[html] wrote {len(viewer_links)} per-complex viewer page(s)")

    # enrichment data
    enrich_summary = {}
    esj = os.path.join(sub["enrichment"], "enrichment_summary.json")
    if os.path.exists(esj):
        try:
            enrich_summary = utils.read_json(esj)
        except Exception:
            enrich_summary = {}

    # accuracy summary
    acc_rows = []
    acsv = os.path.join(sub["reports"], "accuracy_assessment.csv")
    if os.path.exists(acsv):
        with open(acsv, encoding="utf-8-sig", newline="") as fh:
            acc_rows = list(csv.DictReader(fh))

    # selfdock rows
    sd_rows = []
    sdcsv = os.path.join(sub["results"], "selfdock_summary.csv")
    if os.path.exists(sdcsv):
        with open(sdcsv, encoding="utf-8-sig", newline="") as fh:
            sd_rows = list(csv.DictReader(fh))

    # ---- assemble the report ---------------------------------------------
    enrich_html = ""
    recs = enrich_summary.get("receptors") or {}
    if recs:
        sec_rows = []
        for rec, s in recs.items():
            svg_rel = None
            svg_path = os.path.join(sub["enrichment"], s.get("roc_svg", ""))
            if os.path.exists(svg_path):
                svg_rel = os.path.relpath(svg_path, report_dir).replace(os.sep, "/")
            auc_txt = f'{s.get("auc"):.4f}' if isinstance(s.get("auc"), (int, float)) else "n/a"
            sec_rows.append(f"<tr><td>{rec}</td><td>{auc_txt}</td>"
                            f"<td>{s.get('ef1', 'n/a')}</td><td>{s.get('ef5', 'n/a')}</td>"
                            f"<td>{s.get('n_active_scored', 0)}/{s.get('n_active_labels', 0)}</td>"
                            f"<td>{s.get('n_decoy_scored', 0)}/{s.get('n_decoy_labels', 0)}</td></tr>")
        roc_tables = ""
        for rec, s in recs.items():
            svg_path = os.path.join(sub["enrichment"], s.get("roc_svg", ""))
            if os.path.exists(svg_path):
                roc_tables += (f"<h4>ROC - {rec}</h4>"
                               f"<div>{_roc_svg_inline(svg_path)}</div>")
        enrich_html = f"""
<h2 id="enrich">富集度验证(ROC/AUC/EF)</h2>
<p>仅当提供已知活性/诱饵数据时生成;数据来自真实对接得分(refine 优先,否则 screening),标签按规范 SMILES 匹配。方法学详见 md 报告。</p>
<table><thead><tr><th>受体</th><th>AUC</th><th>EF1%</th><th>EF5%</th>
<th>活性(计分/标记)</th><th>诱饵(计分/标记)</th></tr></thead><tbody>
{''.join(sec_rows)}</tbody></table>{roc_tables}"""

    acc_html = ""
    if acc_rows:
        acc_html = f"<h2 id='accuracy'>对接准确性评估</h2>{_html_table(acc_rows)}"
    sd_html = ""
    if sd_rows:
        sd_html = f"<h2 id='selfdock'>共晶配体回贴(Self-docking)</h2>{_html_table(sd_rows)}"

    viewer_block = ("<ul>" + "".join(viewer_links) + "</ul>" if viewer_links
                    else "<p>无 Top-K 复合物可查看(检查 analyze/visualize 是否完成)。</p>")

    topk_html = (f"<h2 id='topk'>最终 Top-{cfg.top_k}</h2>" + _html_table(top_rows)
                 if top_rows else "<h2 id='topk'>最终 Top-K</h2><p>无 final_top CSV。</p>")

    ts = utils.now_str()
    toc = f"""
<h1>DockStudio 交互式总报告</h1>
<p>版本 v{VERSION} · 生成于 {ts} · © 2026 Eric Studio</p>
<ul>
<li><a href="#env">环境</a></li>
<li><a href="#topk">最终 Top-K</a></li>
<li><a href="#selfdock">共晶配体回贴</a></li>
<li><a href="#accuracy">对接准确性评估</a></li>
<li><a href="#viewers">3D 查看器</a></li>
<li><a href="#enrich">富集度验证</a></li>
<li><a href="#qc">质量与文件</a></li>
</ul>"""

    # list report files (md siblings) as QC pointer
    rep_files = []
    if os.path.isdir(sub["reports"]):
        for f in sorted(os.listdir(sub["reports"])):
            if f.endswith((".md", ".csv", ".txt")):
                rep_files.append(f)
    qc_html = ("<h2 id='qc'>质量与文件</h2>"
               "<p>方法与参数:<a href='../01_methods_report.md'>01_methods_report.md</a> · "
               "质量评估:<a href='../02_quality_assessment.md'>02_quality_assessment.md</a></p>"
               "<ul>" + "".join(f"<li>{f}</li>" for f in rep_files) + "</ul>")

    css = """
body{font-family:'Segoe UI','Microsoft YaHei',sans-serif;margin:24px auto;
max-width:1100px;background:#fff;color:#222;line-height:1.5;}
h1{border-bottom:2px solid #2a5db0;padding-bottom:4px;}
h2{margin-top:28px;border-left:4px solid #2a5db0;padding-left:8px;}
table{border-collapse:collapse;font-size:13px;margin:10px 0;overflow-x:auto;display:block;}
th,td{border:1px solid #ddd;padding:4px 8px;text-align:left;}
th{background:#f0f4fa;} li{margin:2px 0;} a{color:#2a5db0;text-decoration:none;}
p{font-size:13px;color:#333;}
"""
    html = f"""<!DOCTYPE html>
<html lang="zh"><head><meta charset="utf-8">
<title>DockStudio 交互式总报告 v{VERSION}</title>
<style>{css}</style>
</head><body>
{toc}
<h2 id="env">环境</h2>
{_env_table(env)}
{topk_html}
{sd_html}
{acc_html}
<h2 id="viewers">3D 查看器(Top-K 复合物)</h2>
{viewer_block}
<p style="color:#888;font-size:12px">查看器用于交互式查看,不构成科学计算依据;相互作用为 PLIP 文本标注。</p>
{enrich_html}
{qc_html}
<hr><p style="color:#888;font-size:12px">本 HTML 报告为 md 报告的姊妹件;科学方法与参数以
<a href="../01_methods_report.md">01_methods_report.md</a> 为准。</p>
</body></html>
"""
    out_path = os.path.join(report_dir, "index.html")
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(html)
    _log(f"[html] total report written: {os.path.relpath(out_path, out_dir)}")
    return {"report_path": out_path, "n_viewers": len(viewer_links)}
