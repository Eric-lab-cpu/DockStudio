"""DockStudio 中文图形界面 (v1.1, modern ttkbootstrap design).

Entry points:
    python -m dockstudio
or the packaged DockStudio.exe.

Design notes (v1.1)
--------------------
* Modern, "top-tier desktop software" style: left navigation rail, header
  action bar, card-based step pages, progress/feedback and consistent spacing.
* Styling is provided by **ttkbootstrap** when installed (recommended); the UI
  degrades gracefully to the classic ``ttk`` look when it is not.
* The output directory is chosen with the native OS folder picker into a
  read-only field (no manual path typing required).
* New in v1.1:
    - receptor inputs accept PDB *and* mmCIF (``.pdb/.ent/.cif/.mmcif``)
    - optional "molecular-docking accuracy assessment report" export
    - Help menu / in-app documentation viewer

Workflow steps
  1) 输入   : receptors (PDB/mmCIF) + ligand library (SDF/MOL) + output folder
  2) 位点   : per-receptor search box (auto cocrystal / manual / residues / blind)
  3) 参数   : scoring & compute budget + accuracy-report switch
  4) 运行   : start / resume / live log
  5) 结果   : Top-K table, reports & accuracy assessment

The engine runs on a background thread; the UI refreshes via a queue.
"""

from __future__ import annotations

import os
import queue
import threading

# --- toolkit bootstrap -----------------------------------------------------
# ttkbootstrap needs tkinter at import time.  If either is unavailable the app
# still tries to fall back to stock tkinter/ttk (e.g. system without theme pkg).
try:
    import tkinter as tk
    from tkinter import filedialog, messagebox
    from tkinter.scrolledtext import ScrolledText
    _TK_OK = True
except Exception:  # pragma: no cover - non-GUI host
    _TK_OK = False

try:
    import ttkbootstrap as ttkb
    from ttkbootstrap import ttk
    HAVE_TTB = True
except Exception:  # pragma: no cover
    ttkb = None
    if _TK_OK:
        from tkinter import ttk
    HAVE_TTB = False

if _TK_OK:
    from dockstudio._version import (
        ACCURACY_REPORT_BASELINE_MEDIUM,
        ACCURACY_REPORT_BASELINE_STRICT,
        APP_NAME,
        APP_NAME_ASCII,
        APP_TITLE,
        BRAND,
        COPYRIGHT_CN,
        __version__,
    )
    from dockstudio.core import envinfo, inventory, models, pipeline, structure as st

OUT_FOLDER_MARKER = "未选择输出目录"
NAV_ITEMS = [
    ("输入", "受体 / 配体 / 输出目录"),
    ("结合位点", "定义搜索盒子"),
    ("参数", "打分与计算预算"),
    ("运行", "启动 / 续跑 / 日志"),
    ("结果", "Top-K 与报告"),
]


# --------------------------------------------------------------------------
# small helpers (kept toolkit-agnostic)
# --------------------------------------------------------------------------
def _stem(path: str) -> str:
    return os.path.splitext(os.path.basename(path))[0]


def _open_path(path: str):
    """Open a file/folder with the OS default application."""
    import subprocess
    import sys
    try:
        if sys.platform.startswith("win"):
            os.startfile(path)  # type: ignore[attr-defined]
        elif sys.platform == "darwin":
            subprocess.Popen(["open", path])
        else:
            subprocess.Popen(["xdg-open", path])
    except Exception:
        return False
    return True


def _scan_polymer_chains(path: str) -> list:
    try:
        inv = st.inventory_pdb(path)
        return [c for c, ci in inv.chains.items() if ci.polymer_res > 0]
    except Exception:
        return []


def _find_cocrystal(path: str, chains: list) -> list:
    try:
        rows = inventory.inventory_receptor(path, _stem(path), chains)
        found = []
        for r in rows:
            if r.get("ligand") and r["chain"] in chains:
                found.append(f"{r['ligand']} ({r['chain']}{r.get('ligand_resseq')}, "
                             f"{r.get('ligand_n_heavy')} heavy)")
        seen = set()
        out = []
        for x in found:
            if x not in seen:
                seen.add(x)
                out.append(x)
        return out
    except Exception:
        return []


def _safe_int(var: tk.Variable, default: int) -> int:
    try:
        return int(float(var.get()))
    except Exception:
        return default


def _safe_float(var: tk.Variable, default: float) -> float:
    try:
        return float(var.get())
    except Exception:
        return default


# --------------------------------------------------------------------------
# main application
# --------------------------------------------------------------------------
class DockStudioApp:
    def __init__(self, root):
        self.root = root
        root.title(f"{APP_TITLE}  v{__version__}")
        try:
            root.geometry("1240x840")
            root.minsize(1080, 720)
        except Exception:
            pass
        self._icon = None
        self._apply_icon()

        self._q: queue.Queue = queue.Queue()
        self._worker: threading.Thread | None = None
        self._stop_event = threading.Event()

        # ---- data model ---------------------------------------------------
        self.receptors = []      # list[dict]: path, name, chains(list)
        self.ligands = []        # list[dict]: path, name
        self.out_dir = ""
        self.boxes = {}          # rec name -> dict(method, center, size, basis, ligand_ref)

        self._build_menu()
        self._build_ui()
        self._build_footer()
        self.root.after(120, self._poll_queue)

        # environment check note (only block on required tools)
        self._post_env_check()

    # ================================================================== helpers
    def _resource_path(self, name: str) -> str:
        import sys
        bundle = getattr(sys, "_MEIPASS", None)
        if bundle:
            cand = os.path.join(bundle, "dockstudio", "resources", name)
            if os.path.exists(cand):
                return cand
        here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # dockstudio/
        cand = os.path.join(here, "resources", name)
        if os.path.exists(cand):
            return cand
        return os.path.join(os.path.dirname(os.path.dirname(here)), "assets", name)

    def _apply_icon(self):
        try:
            ico = self._resource_path("icon.ico")
            if os.name == "nt" and os.path.exists(ico):
                self.root.iconbitmap(default=ico)
            else:
                png = self._resource_path("icon.png")
                if os.path.exists(png):
                    from PIL import Image as _PILImage
                    from PIL import ImageTk as _ImageTk
                    img = _ImageTk.PhotoImage(_PILImage.open(png).resize((32, 32)))
                    self._icon = img
                    self.root.iconphoto(True, img)
        except Exception:
            pass

    def _post_env_check(self):
        try:
            env = envinfo.gather_environment()
            missing = env.get("missing_required") or []
            miss_opt = env.get("missing_optional") or []
            if missing:
                messagebox.showwarning(
                    "依赖提示",
                    "以下必需组件未检测到,相关功能将不可用:\n" + ", ".join(missing) +
                    "\n\n请按 帮助 → 用户手册 中的安装说明配置依赖。",
                    parent=self.root)
            elif miss_opt:
                note = "可选组件未检测到(仅降级,不影响主流程): " + ", ".join(miss_opt)
                self._q.put(("log", note))
        except Exception:
            pass

    # ================================================================== widgets
    def _btn(self, parent, text, command, bootstyle="secondary", width=None):
        kw = {"text": text, "command": command}
        if width:
            kw["width"] = width
        if HAVE_TTB:
            return ttk.Button(parent, bootstyle=bootstyle, **kw)
        return ttk.Button(parent, **kw)

    def _card(self, parent, title: str):
        """A 'card': a labeled container with themed background."""
        if HAVE_TTB:
            lf = ttk.Labelframe(parent, text=title, padding=(10, 8), bootstyle="primary")
        else:
            lf = ttk.LabelFrame(parent, text=title, padding=(10, 8))
        return lf

    # ------------------------------------------------------------- menu
    def _build_menu(self):
        menubar = tk.Menu(self.root)
        m_file = tk.Menu(menubar, tearoff=0)
        m_file.add_command(label="打开输出目录", command=self._open_out)
        m_file.add_separator()
        m_file.add_command(label="退出", command=self.root.destroy)
        menubar.add_cascade(label="文件", menu=m_file)

        m_run = tk.Menu(menubar, tearoff=0)
        m_run.add_command(label="开始 / 断点续跑", command=self._start_run)
        m_run.add_command(label="停止(阶段结束后)", command=self._request_stop)
        menubar.add_cascade(label="运行", menu=m_run)

        m_help = tk.Menu(menubar, tearoff=0)
        m_help.add_command(label="帮助中心", command=self._open_help)
        m_help.add_command(label="关于 " + APP_NAME, command=self._about)
        menubar.add_cascade(label="帮助", menu=m_help)
        self.root.config(menu=menubar)

    def _about(self):
        top = tk.Toplevel(self.root)
        top.title("关于 DockStudio")
        top.geometry("520x330")
        top.transient(self.root)
        top.resizable(False, False)
        box = ttk.Frame(top, padding=20)
        box.pack(fill="both", expand=True)
        ttk.Label(box, text=APP_NAME, font=("Microsoft YaHei UI", 16, "bold")).pack(anchor="center")
        ttk.Label(box, text=f"版本 {__version__} · {APP_NAME_ASCII}", font=("", 11)).pack(pady=(6, 2))
        ttk.Separator(box).pack(fill="x", pady=10)
        ttk.Label(box, text="一站式自动化批量分子对接平台\n"
                            "AutoDock Vina · Meeko · RDKit · PLIP · PyMOL",
                  justify="center").pack(pady=4)
        ttk.Label(box, text=f"本产品由 {BRAND} 提供",
                  font=("", 11, "bold")).pack(pady=(14, 2))
        ttk.Label(box, text=COPYRIGHT_CN, foreground="#777").pack()
        if HAVE_TTB:
            ttk.Button(top, text="关闭", bootstyle="secondary", command=top.destroy).pack(pady=12)
        else:
            ttk.Button(top, text="关闭", command=top.destroy).pack(pady=12)

    # ------------------------------------------------------------- shell
    def _build_ui(self):
        """Header + body (left rail nav + notebook content) + footer."""
        self._build_header()

        body = ttk.Frame(self.root, padding=(6, 2))
        body.pack(fill="both", expand=True)

        # left navigation rail
        nav = ttk.Frame(body, width=206)
        nav.pack(side="left", fill="y", padx=(0, 6))
        nav.pack_propagate(False)
        ttk.Label(nav, text="流程向导", font=("Microsoft YaHei UI", 9, "bold"),
                  foreground="#666").pack(anchor="w", padx=10, pady=(4, 4))
        self._nav_btns = []
        for i, (title, sub) in enumerate(NAV_ITEMS):
            b = ttk.Button(nav, text=f"{i + 1}.  {title}", command=lambda k=i: self._goto(k))
            b.pack(fill="x", padx=6, pady=2)
            self._nav_btns.append(b)

        # content notebook (tabs kept for page management, visually hidden)
        self.nb = ttk.Notebook(body)
        self.nb.pack(side="left", fill="both", expand=True)
        self.tab_inputs = ttk.Frame(self.nb, padding=8)
        self.tab_site = ttk.Frame(self.nb, padding=8)
        self.tab_params = ttk.Frame(self.nb, padding=8)
        self.tab_run = ttk.Frame(self.nb, padding=8)
        self.tab_results = ttk.Frame(self.nb, padding=8)
        for pg, name in ((self.tab_inputs, "inputs"), (self.tab_site, "site"),
                         (self.tab_params, "params"), (self.tab_run, "run"),
                         (self.tab_results, "results")):
            self.nb.add(pg, text="")
        # hide tab strip (rail navigation instead)
        try:
            self.nb.configure(tabposition="nw")
        except Exception:
            pass
        # hide the notebook tab strip (navigation is via the left rail)
        try:
            if HAVE_TTB:
                _style = ttkb.Style()
            else:
                _style = ttk.Style(self.root)
            _style.layout("TNotebook.Tab", [])  # make tabs invisible
        except Exception:
            pass

        self._build_inputs_tab()
        self._build_site_tab()
        self._build_params_tab()
        self._build_run_tab()
        self._build_results_tab()

        self._goto(0)
        self.nb.bind("<<NotebookTabChanged>>", self._on_tab_changed)

    def _build_header(self):
        hdr = ttk.Frame(self.root, padding=(12, 8))
        hdr.pack(fill="x")
        ttk.Label(hdr, text=APP_NAME, font=("Microsoft YaHei UI", 15, "bold")).pack(side="left")
        ttk.Label(hdr, text=f"v{__version__}  ·  严谨可复现的批量分子对接",
                  foreground="#888").pack(side="left", padx=(10, 0), pady=(4, 0))
        if HAVE_TTB:
            ttk.Button(hdr, text="开始运行", bootstyle="success",
                       command=self._start_run).pack(side="right", padx=4)
            ttk.Button(hdr, text="打开输出目录", bootstyle="secondary-outline",
                       command=self._open_out).pack(side="right", padx=4)
            ttk.Button(hdr, text="帮助", bootstyle="info-outline",
                       command=self._open_help).pack(side="right", padx=4)
        else:
            ttk.Button(hdr, text="开始运行", command=self._start_run).pack(side="right", padx=4)
            ttk.Button(hdr, text="打开输出目录", command=self._open_out).pack(side="right", padx=4)
            ttk.Button(hdr, text="帮助", command=self._open_help).pack(side="right", padx=4)

    def _build_footer(self):
        bar = ttk.Frame(self.root)
        bar.pack(fill="x", padx=10, pady=(0, 4))
        ttk.Label(bar, text=f"{BRAND} · 仅供科研使用 · 结果请按发表规范复核",
                  foreground="#777").pack(side="left")
        self._foot_status = tk.StringVar(value="就绪")
        ttk.Label(bar, textvariable=self._foot_status, foreground="#888").pack(side="right")

    def _goto(self, idx: int):
        try:
            self.nb.select(idx)
        except Exception:
            return
        for i, b in enumerate(self._nav_btns):
            try:
                if i == idx and HAVE_TTB:
                    b.configure(bootstyle="primary")
                elif HAVE_TTB:
                    b.configure(bootstyle="secondary")
            except Exception:
                pass
        self._on_tab_changed()

    def _on_tab_changed(self, event=None):
        try:
            idx = self.nb.index(self.nb.select())
        except Exception:
            return
        if idx == 1:
            self._refresh_site_lists()

    # ------------------------------------------------------------- tab 1
    def _build_inputs_tab(self):
        f = self.tab_inputs

        # Output directory card (requirement: native folder picker, no typing)
        od = self._card(f, "输出目录")
        od.pack(fill="x", pady=(0, 8))
        orow = ttk.Frame(od)
        orow.pack(fill="x")
        self.out_var = tk.StringVar(value=OUT_FOLDER_MARKER)
        ent = ttk.Entry(orow, textvariable=self.out_var, state="readonly")
        ent.pack(side="left", fill="x", expand=True, padx=(0, 6))
        self._btn(orow, "浏览文件夹…", self._choose_out, bootstyle="primary").pack(side="left")
        self._btn(orow, "示例: 4DFR 演示", self._load_demo, bootstyle="secondary-outline").pack(side="left", padx=6)
        ttk.Label(od, text="所有结果、报告与图件都将写入此文件夹。",
                  foreground="#888").pack(anchor="w", pady=(4, 0))

        mid = ttk.Frame(f)
        mid.pack(fill="both", expand=True)
        mid.columnconfigure(0, weight=1)
        mid.columnconfigure(1, weight=1)
        mid.rowconfigure(0, weight=1)

        # --- receptor card -------------------------------------------------
        left = self._card(mid, "受体(大分子)— PDB / mmCIF,可添加多个")
        left.grid(row=0, column=0, sticky="nsew", padx=(0, 4))
        left.columnconfigure(0, weight=1)
        left.rowconfigure(1, weight=1)
        self.rec_list = ttk.Treeview(left, columns=("chains",), height=9)
        self.rec_list.heading("#0", text="受体文件")
        self.rec_list.heading("chains", text="对接链")
        self.rec_list.column("#0", width=280)
        self.rec_list.column("chains", width=110)
        self.rec_list.grid(row=1, column=0, sticky="nsew", padx=2, pady=4)
        bf = ttk.Frame(left)
        bf.grid(row=2, column=0, sticky="w")
        self._btn(bf, "添加受体…", self._add_receptor, bootstyle="primary").pack(side="left")
        self._btn(bf, "删除选中", self._del_receptor).pack(side="left", padx=4)
        self._btn(bf, "上移", lambda: self._move_receptor(-1)).pack(side="left")
        self._btn(bf, "下移", lambda: self._move_receptor(1)).pack(side="left", padx=4)
        self._btn(bf, "链 / 共晶…", self._edit_receptor).pack(side="left", padx=4)

        # --- ligand card ---------------------------------------------------
        right = self._card(mid, "配体库(SDF / MOL / SMILES / CSV)— 可添加多个文件")
        right.grid(row=0, column=1, sticky="nsew", padx=(4, 0))
        right.columnconfigure(0, weight=1)
        right.rowconfigure(1, weight=1)
        self.lig_list = ttk.Treeview(right, columns=("name",), height=9)
        self.lig_list.heading("#0", text="配体库文件")
        self.lig_list.heading("name", text="库名")
        self.lig_list.column("#0", width=280)
        self.lig_list.column("name", width=110)
        self.lig_list.grid(row=1, column=0, sticky="nsew", padx=2, pady=4)
        lf = ttk.Frame(right)
        lf.grid(row=2, column=0, sticky="w")
        self._btn(lf, "添加 SDF/SMILES…", self._add_ligand, bootstyle="primary").pack(side="left")
        self._btn(lf, "删除选中", self._del_ligand).pack(side="left", padx=4)

    def _add_receptor(self):
        paths = filedialog.askopenfilenames(
            parent=self.root, title="选择受体大分子文件",
            filetypes=[("结构文件", "*.pdb *.ent *.cif *.mmcif"),
                       ("PDB", "*.pdb *.ent"), ("mmCIF", "*.cif *.mmcif"),
                       ("所有文件", "*.*")])
        for p in paths:
            if any(r["path"] == p for r in self.receptors):
                continue
            chains = _scan_polymer_chains(p)
            self.receptors.append({"path": p, "name": _stem(p), "chains": chains})
        self._refresh_rec_list()

    def _del_receptor(self):
        sel = self.rec_list.selection()
        if not sel:
            return
        idx = int(sel[0])
        name = self.receptors[idx]["name"]
        self.receptors.pop(idx)
        self.boxes.pop(name, None)
        self._refresh_rec_list()

    def _move_receptor(self, d: int):
        sel = self.rec_list.selection()
        if not sel:
            return
        i = int(sel[0])
        j = i + d
        if 0 <= j < len(self.receptors):
            self.receptors[i], self.receptors[j] = self.receptors[j], self.receptors[i]
            self._refresh_rec_list()
            self.rec_list.selection_set(str(j))

    def _edit_receptor(self):
        sel = self.rec_list.selection()
        if not sel:
            return
        i = int(sel[0])
        rec = self.receptors[i]
        dlg = tk.Toplevel(self.root)
        dlg.title(f"受体设置 - {rec['name']}")
        dlg.geometry("480x380")
        dlg.transient(self.root)
        ttk.Label(dlg, text="对接链(可多选, Ctrl/Shift):").pack(anchor="w", padx=10, pady=4)
        lb = tk.Listbox(dlg, selectmode="extended", height=11, exportselection=False)
        all_chains = _scan_polymer_chains(rec["path"]) or (rec["chains"] or ["A"])
        for c in all_chains:
            lb.insert("end", c)
        for j, c in enumerate(all_chains):
            if c in (rec.get("chains") or []):
                lb.selection_set(j)
        lb.pack(fill="both", expand=True, padx=10)
        coc = _find_cocrystal(rec["path"], all_chains)
        info = ("检测到共晶配体:\n" + "\n".join("  - " + c for c in coc)) if coc else \
            "未在所选链中检测到共晶配体。"
        ttk.Label(dlg, text=info, justify="left").pack(anchor="w", padx=10, pady=4)

        def apply():
            rec["chains"] = [all_chains[x] for x in lb.curselection()] or all_chains
            if rec["name"] in self.boxes:
                self.boxes[rec["name"]]["basis"] = "auto (edited chains)"
            self._refresh_rec_list()
            dlg.destroy()

        self._btn(dlg, "确定", apply, bootstyle="primary").pack(pady=8)

    def _add_ligand(self):
        paths = filedialog.askopenfilenames(
            parent=self.root, title="选择配体库",
            filetypes=[("配体文件", "*.sdf *.mol *.sd *.smi *.smiles *.csv *.txt"),
                       ("SDF/MOL", "*.sdf *.mol *.sd"),
                       ("SMILES/CSV", "*.smi *.smiles *.csv *.txt"),
                       ("所有文件", "*.*")])
        for p in paths:
            if any(l["path"] == p for l in self.ligands):
                continue
            self.ligands.append({"path": p, "name": _stem(p)})
        self._refresh_lig_list()

    def _del_ligand(self):
        sel = self.lig_list.selection()
        if sel:
            self.ligands.pop(int(sel[0]))
            self._refresh_lig_list()

    def _pick_enrich_file(self, kind: str):
        p = filedialog.askopenfilename(
            parent=self.root, title=("选择已知活性文件" if kind == "actives" else "选择诱饵文件"),
            filetypes=[("配体文件", "*.sdf *.mol *.sd *.smi *.smiles *.csv *.txt"),
                       ("SDF/MOL", "*.sdf *.mol *.sd"),
                       ("SMILES/CSV", "*.smi *.smiles *.csv *.txt"),
                       ("所有文件", "*.*")])
        if p:
            if kind == "actives":
                self.act_var.set(p)
            else:
                self.dec_var.set(p)
            self._q.put(("log", f"富集度 {kind} 文件: {p}"))

    def _choose_out(self):
        p = filedialog.askdirectory(parent=self.root, title="选择输出目录")
        if p:
            self.out_dir = p
            self.out_var.set(p)
            self._q.put(("log", f"输出目录: {p}"))

    def _load_demo(self):
        demo = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__)))), "examples", "demo")
        pdb = os.path.join(demo, "4DFR.pdb")
        lib = os.path.join(demo, "ligand_library.sdf")
        if not os.path.exists(lib):
            try:
                import sys
                sys.path.insert(0, demo)
                import make_library
                make_library.build(lib)
            except Exception as e:
                messagebox.showerror("示例", f"构建示例配体库失败:{e}", parent=self.root)
                return
        self.receptors = [{"path": pdb, "name": "4DFR", "chains": ["A"]}]
        self.ligands = [{"path": lib, "name": "demo"}]
        self.out_dir = os.path.join(os.path.dirname(demo), "demo_out")
        self.out_var.set(self.out_dir)
        self.boxes.pop("4DFR", None)
        self._refresh_rec_list()
        self._refresh_lig_list()
        messagebox.showinfo("示例", "已载入 4DFR 演示数据。请在 ②结合位点 查看自动盒子,然后 ④运行。",
                            parent=self.root)

    def _refresh_rec_list(self):
        self.rec_list.delete(*self.rec_list.get_children())
        for i, r in enumerate(self.receptors):
            self.rec_list.insert("", "end", iid=str(i), text=os.path.basename(r["path"]),
                                 values=(",".join(r.get("chains") or []),))

    def _refresh_lig_list(self):
        self.lig_list.delete(*self.lig_list.get_children())
        for i, l in enumerate(self.ligands):
            self.lig_list.insert("", "end", iid=str(i), text=os.path.basename(l["path"]),
                                 values=(l["name"],))

    # ------------------------------------------------------------- tab 2
    def _build_site_tab(self):
        f = self.tab_site
        hint = ("每个受体需定义搜索盒子。选择受体后按共晶配体 / 残基自动或手动定义;"
                "无共晶配体时务必选择盲对接(探索性)或手动定义。")
        ttk.Label(f, text=hint, wraplength=900, foreground="#666").pack(anchor="w")

        top = ttk.Frame(f)
        top.pack(fill="x", pady=6)
        ttk.Label(top, text="受体:").pack(side="left")
        self.site_rec_var = tk.StringVar()
        self.site_combo = ttk.Combobox(top, textvariable=self.site_rec_var, width=32,
                                       state="readonly")
        self.site_combo.pack(side="left", padx=4)
        self.site_combo.bind("<<ComboboxSelected>>", lambda e: self._load_site_for_rec())
        self._btn(top, "刷新(重读 PDB/CIF)", self._refresh_site_lists).pack(side="left", padx=4)
        self._btn(top, "预览盒子", self._preview_box, bootstyle="info").pack(side="left", padx=4)

        body = ttk.Frame(f)
        body.pack(fill="both", expand=True)
        body.columnconfigure(0, weight=1)
        body.rowconfigure(1, weight=1)

        self.site_info = tk.Text(body, height=6, state="disabled",
                                 bg="#f6f8fa", relief="flat", font=("Consolas", 9))
        self.site_info.grid(row=0, column=0, sticky="ew", pady=(6, 4))

        stylef = self._card(body, "盒子定义方法")
        stylef.grid(row=1, column=0, sticky="nsew", pady=4)
        self.box_mode = tk.StringVar(value="auto")
        opts = [
            ("auto", "自动: 共晶配体质心(推荐,仅当所选链含共晶配体)"),
            ("explicit", "手动: 中心 + 尺寸"),
            ("residues", "残基中心(如 TRP168 或 168,可逗号分隔)"),
            ("blind", "盲对接(全蛋白,探索性,高假阳性)"),
        ]
        for i, (val, txt) in enumerate(opts):
            ttk.Radiobutton(stylef, text=txt, value=val, variable=self.box_mode,
                            command=self._on_box_mode).pack(anchor="w", pady=2)

        gridf = ttk.Frame(f)
        gridf.pack(fill="x", pady=4)
        self.cx = tk.DoubleVar(value=0.0)
        self.cy = tk.DoubleVar(value=0.0)
        self.cz = tk.DoubleVar(value=0.0)
        self.sx = tk.DoubleVar(value=24.0)
        self.sy = tk.DoubleVar(value=24.0)
        self.sz = tk.DoubleVar(value=24.0)
        for i, (lab, var) in enumerate([("中心 X", self.cx), ("中心 Y", self.cy),
                                        ("中心 Z", self.cz)]):
            ttk.Label(gridf, text=lab).grid(row=0, column=i * 2, sticky="e", padx=(2, 2))
            ttk.Entry(gridf, textvariable=var, width=9).grid(row=0, column=i * 2 + 1, padx=2)
        for i, (lab, var) in enumerate([("尺寸 X", self.sx), ("尺寸 Y", self.sy),
                                        ("尺寸 Z", self.sz)]):
            ttk.Label(gridf, text=lab).grid(row=1, column=i * 2, sticky="e", padx=(2, 2))
            ttk.Entry(gridf, textvariable=var, width=9).grid(row=1, column=i * 2 + 1, padx=2)
        ttk.Label(gridf, text="尺寸下限 8 Å,建议 >= 22 Å(默认 24)。").grid(
            row=2, column=0, columnspan=6, sticky="w", pady=(4, 0))
        self.site_res_var = tk.StringVar(value="")
        rr = ttk.Frame(f)
        rr.pack(fill="x", pady=4)
        ttk.Label(rr, text="残基列表:").pack(side="left")
        ttk.Entry(rr, textvariable=self.site_res_var, width=46).pack(side="left", padx=4)

    def _on_box_mode(self):
        """Enable/disable manual inputs depending on the chosen box method."""
        mode = self.box_mode.get()
        # manual center/size only relevant for 'explicit'; residue list for residues
        # (kept simple: always visible; validation happens at compute time)

    def _refresh_site_lists(self):
        names = [r["name"] for r in self.receptors]
        self.site_combo["values"] = names
        if names:
            cur = self.site_rec_var.get()
            if cur not in names:
                self.site_rec_var.set(names[0])
            self._load_site_for_rec()

    def _load_site_for_rec(self):
        name = self.site_rec_var.get()
        if not name:
            return
        rec = next((r for r in self.receptors if r["name"] == name), None)
        if not rec:
            return
        box = self.boxes.get(name)
        chains = rec.get("chains") or []
        lines = [f"受体: {name}   对接链: {','.join(chains) or '(全部)'}"]
        coc = _find_cocrystal(rec["path"], chains)
        lines.append("共晶配体: " + ("; ".join(coc) if coc else "无(请手动定义,或选盲对接)"))
        if box:
            lines.append(f"现有盒子: method={box.get('method')} center={box.get('center')} "
                         f"size={box.get('size')}")
        self.site_info.config(state="normal")
        self.site_info.delete("1.0", "end")
        self.site_info.insert("1.0", "\n".join(lines))
        self.site_info.config(state="disabled")

    def _preview_box(self):
        name = self.site_rec_var.get()
        if not name:
            messagebox.showinfo("提示", "请先在 ①输入 添加受体。", parent=self.root)
            return
        rec = next((r for r in self.receptors if r["name"] == name), None)
        if not rec:
            return
        box = self._compute_box_for(rec, self.box_mode.get(), silent=False)
        if box is None:
            return
        self.boxes[name] = box
        self._load_site_for_rec()
        messagebox.showinfo("盒子",
                            f"{name}\n方法: {box['method']}\n"
                            f"中心: {[round(x, 2) for x in box['center']]}\n"
                            f"尺寸: {[round(x, 2) for x in box['size']]}\n{box.get('basis', '')}",
                            parent=self.root)

    def _compute_box_for(self, rec, mode, silent=True):
        from dockstudio.core import box as boxmod
        pdb = rec["path"]
        chains = rec.get("chains") or []
        if mode == "auto":
            rows = inventory.inventory_receptor(pdb, rec["name"], chains)
            cand = None
            for r in rows:
                if r.get("ligand") and r["chain"] in chains:
                    if cand is None or (r.get("ligand_n_heavy") or 0) > \
                            (cand.get("ligand_n_heavy") or 0):
                        cand = r
            if cand:
                ligr = st.coords_of_ligand(pdb, cand["ligand"], cand["chain"],
                                           cand.get("ligand_resseq"))
                b = boxmod.box_from_ligand(ligr, ligand_tag=cand["ligand"])
                return b.to_dict()
            if not silent:
                messagebox.showwarning("位点", f"{rec['name']}: 所选链中无共晶配体,"
                                       "请手动定义或选择盲对接。", parent=self.root)
            return None
        if mode == "explicit":
            center = [self.cx.get(), self.cy.get(), self.cz.get()]
            size = [self.sx.get(), self.sy.get(), self.sz.get()]
            if min(size) < 8:
                messagebox.showwarning("参数", "盒子尺寸过小(<8 Å)", parent=self.root)
                return None
            return {"method": "explicit", "center": center, "size": size,
                    "basis": "user-specified center+size"}
        if mode == "residues":
            txt = self.site_res_var.get().strip()
            if not txt:
                messagebox.showwarning("位点", "请填写残基列表", parent=self.root)
                return None
            tokens = [t.strip() for t in txt.replace(";", ",").split(",")]
            recs = st.parse_structure(pdb, keep_chains=chains or None)
            try:
                b = boxmod.box_from_residues(recs, tokens)
                return b.to_dict()
            except Exception as e:
                messagebox.showerror("残基", str(e), parent=self.root)
                return None
        if mode == "blind":
            recs = st.parse_structure(pdb, keep_chains=chains or None)
            b = boxmod.box_blind(recs)
            return b.to_dict()
        return None

    def _collect_boxes(self):
        """Ensure every receptor has a box; compute auto/default where absent."""
        from dockstudio.core import box as boxmod
        missing = []
        for rec in self.receptors:
            name = rec["name"]
            if name in self.boxes:
                continue
            rows = inventory.inventory_receptor(rec["path"], name, rec.get("chains") or [])
            cand = None
            for r in rows:
                if r.get("ligand") and r["chain"] in (rec.get("chains") or []):
                    if cand is None or (r.get("ligand_n_heavy") or 0) > \
                            (cand.get("ligand_n_heavy") or 0):
                        cand = r
            if cand:
                ligr = st.coords_of_ligand(rec["path"], cand["ligand"], cand["chain"],
                                           cand.get("ligand_resseq"))
                self.boxes[name] = boxmod.box_from_ligand(
                    ligr, ligand_tag=cand["ligand"]).to_dict()
            else:
                missing.append(name)
        return missing

    # ------------------------------------------------------------- tab 3
    def _build_params_tab(self):
        f = self.tab_params
        p = self._card(f, "计算参数(默认与协议一致;可下调以省时)")
        p.pack(fill="x")
        grid = ttk.Frame(p)
        grid.pack(fill="x", padx=4, pady=4)
        self.exh = tk.IntVar(value=16)
        self.exh_refine = tk.IntVar(value=32)
        self.nposes = tk.IntVar(value=9)
        self.erange = tk.DoubleVar(value=4.0)
        self.top_ref = tk.IntVar(value=12)
        self.topk = tk.IntVar(value=5)
        self.cpu = tk.IntVar(value=2)
        self.ph = tk.DoubleVar(value=7.4)
        fields = [
            ("初筛 exhaustiveness", self.exh),
            ("精修 exhaustiveness", self.exh_refine),
            ("每配体 pose 数", self.nposes),
            ("energy_range (kcal/mol)", self.erange),
            ("精修短名单", self.top_ref),
            ("最终 Top-K", self.topk),
            ("并行 CPU", self.cpu),
            ("配体目标 pH", self.ph),
        ]
        for i, (lab, var) in enumerate(fields):
            r, c = divmod(i, 2)
            ttk.Label(grid, text=lab).grid(row=r, column=c * 2, sticky="e", padx=(2, 3), pady=3)
            ttk.Entry(grid, textvariable=var, width=10).grid(row=r, column=c * 2 + 1, padx=3)

        toggles = self._card(f, "开关")
        toggles.pack(fill="x", pady=8)
        self.chk_self = tk.BooleanVar(value=True)
        self.chk_plip = tk.BooleanVar(value=True)
        self.chk_viz = tk.BooleanVar(value=True)
        self.chk_refine = tk.BooleanVar(value=True)
        self.chk_pse = tk.BooleanVar(value=True)
        self.chk_acc = tk.BooleanVar(value=True)
        row1 = ttk.Frame(toggles)
        row1.pack(fill="x", padx=4, pady=2)
        for var, txt in [(self.chk_self, "回贴验证(self-docking)"),
                         (self.chk_refine, "精修(Top-K 前)"),
                         (self.chk_acc, "导出对接准确性评估报告")]:
            ttk.Checkbutton(row1, text=txt, variable=var).pack(side="left", padx=8)
        row2 = ttk.Frame(toggles)
        row2.pack(fill="x", padx=4, pady=2)
        for var, txt in [(self.chk_plip, "相互作用分析 PLIP"),
                         (self.chk_viz, "PyMOL 图件 + .pml"),
                         (self.chk_pse, "输出 .pse 会话")]:
            ttk.Checkbutton(row2, text=txt, variable=var).pack(side="left", padx=8)

        # --- v2.0: throughput / enrichment validation card ------------------
        v2 = self._card(f, "v2.0: 高通量并行 + 富集度验证")
        v2.pack(fill="x", pady=8)
        row0 = ttk.Frame(v2)
        row0.pack(fill="x", padx=6, pady=3)
        self.workers = tk.IntVar(value=1)
        ttk.Label(row0, text="并行任务数(进程池; >1 时每个 dock 用“并行 CPU”线程):").pack(side="left")
        ttk.Entry(row0, textvariable=self.workers, width=8).pack(side="left", padx=4)
        self.chk_enrich = tk.BooleanVar(value=False)
        ttk.Checkbutton(v2, text="启用“富集度验证模式”(ROC/AUC/EF;需提供已知活性 + 诱饵数据)",
                        variable=self.chk_enrich).pack(anchor="w", padx=6, pady=2)
        self.act_var = tk.StringVar(value="")
        self.dec_var = tk.StringVar(value="")
        actrow = ttk.Frame(v2)
        actrow.pack(fill="x", padx=6, pady=1)
        ttk.Label(actrow, text="已知活性文件(SDF/SMILES/CSV):").pack(side="left")
        ttk.Entry(actrow, textvariable=self.act_var, width=30).pack(side="left", padx=4)
        self._btn(actrow, "浏览…", lambda: self._pick_enrich_file("actives"),
                  bootstyle="secondary").pack(side="left")
        decrow = ttk.Frame(v2)
        decrow.pack(fill="x", padx=6, pady=1)
        ttk.Label(decrow, text="诱饵文件(SDF/SMILES/CSV):").pack(side="left")
        ttk.Entry(decrow, textvariable=self.dec_var, width=30).pack(side="left", padx=4)
        self._btn(decrow, "浏览…", lambda: self._pick_enrich_file("decoys"),
                  bootstyle="secondary").pack(side="left")
        ttk.Label(v2, justify="left", wraplength=940, foreground="#666",
                  text="说明:富集度模式要求 active/decoy 分子也包含在你的配体库中并成功对接;"
                       "引擎按规范 SMILES 匹配打分并计算 AUC/EF1%/5%,绝不自动生成活性数据。"
                       "未提供完整数据时本模式不产出 ROC 结果。").pack(anchor="w", padx=6, pady=(2, 2))

        note = self._card(f, "说明(科学严谨性)")
        note.pack(fill="x", pady=(0, 4))
        ttk.Label(note, justify="left", wraplength=940, foreground="#555",
                  text=(
                      "• 自对接(self-docking)用共晶配体重对接评估几何精度;勾选“导出对接准确性"
                      "评估报告”后,将额外生成 reports/03_accuracy_assessment.md 与 "
                      "accuracy_assessment.csv(mode-1 RMSD,PASS ≤ 2.0 Å,高精度 ≤ 1.0 Å)。\n"
                      "• 无共晶配体受体无法做回贴验证,结果一律标注为探索性。\n"
                      "• 质子化为文档化简单规则(非 pKa 预测),目标 pH 如上可调。")).pack(anchor="w")

    def _build_cfg(self) -> models.RunConfig:
        missing = self._collect_boxes()
        recs = [{"path": r["path"], "name": r["name"], "chains": r.get("chains") or []}
                for r in self.receptors]
        ligs = [{"path": l["path"], "name": l["name"]} for l in self.ligands]
        return models.RunConfig(
            receptors=recs, ligands=ligs, out_dir=self.out_dir,
            title="DockStudio GUI run",
            boxes=self.boxes,
            exhaustiveness=_safe_int(self.exh, 16),
            refine_exhaustiveness=_safe_int(self.exh_refine, 32),
            n_poses=_safe_int(self.nposes, 9),
            energy_range=_safe_float(self.erange, 4.0),
            top_refine=_safe_int(self.top_ref, 12),
            top_k=_safe_int(self.topk, 5),
            cpu=_safe_int(self.cpu, 2),
            ph=_safe_float(self.ph, 7.4),
            run_selfdock=self.chk_self.get(),
            run_plip=self.chk_plip.get(),
            run_visualization=self.chk_viz.get(),
            run_refine=self.chk_refine.get(),
            write_pse=self.chk_pse.get(),
            run_accuracy_report=self.chk_acc.get(),
            n_workers=_safe_int(self.workers, 1) if hasattr(self, "workers") else 1,
            run_enrichment=self.chk_enrich.get() if hasattr(self, "chk_enrich") else False,
            actives_path=self.act_var.get().strip() if hasattr(self, "act_var") else "",
            decoys_path=self.dec_var.get().strip() if hasattr(self, "dec_var") else "",
            run_html_report=True,
            use_symmetry_rmsd=True,
            delete_bad_res=True,
            overwrite=False)

    # ------------------------------------------------------------- tab 4
    def _build_run_tab(self):
        f = self.tab_run
        top = ttk.Frame(f)
        top.pack(fill="x", pady=4)
        self.start_btn = self._btn(top, "开始 / 断点续跑", self._start_run, bootstyle="success")
        self.start_btn.pack(side="left")
        self._btn(top, "停止(阶段结束后)", self._request_stop, bootstyle="warning").pack(side="left", padx=6)
        self._btn(top, "打开输出目录", self._open_out, bootstyle="secondary-outline").pack(side="left")
        ttk.Label(top, text="提示:中途关闭后再次点击“开始 / 断点续跑”即可从断点继续。",
                  foreground="#888").pack(side="left", padx=12)

        self.progress = ttk.Progressbar(f, mode="determinate")
        self.progress.pack(fill="x", pady=6)
        self.log_text = ScrolledText(f, height=16, state="disabled",
                                     font=("Consolas", 9), relief="flat")
        self.log_text.pack(fill="both", expand=True, pady=(4, 4))
        self.status_var = tk.StringVar(value="就绪")
        ttk.Label(f, textvariable=self.status_var, anchor="w").pack(fill="x")

    def _open_out(self):
        d = self.out_var.get()
        if d and d != OUT_FOLDER_MARKER and os.path.isdir(d):
            if not _open_path(d):
                messagebox.showinfo("路径", d, parent=self.root)

    def _log(self, msg: str):
        self._q.put(("log", msg))

    def _start_run(self):
        if self._worker and self._worker.is_alive():
            return
        if not self.receptors or not self.ligands:
            messagebox.showwarning("输入", "请先添加受体与配体库。", parent=self.root)
            return
        out = self.out_var.get()
        if not out or out == OUT_FOLDER_MARKER:
            messagebox.showwarning("输出", "请选择输出目录。", parent=self.root)
            return
        if not os.path.isdir(out):
            try:
                os.makedirs(out, exist_ok=True)
            except Exception as e:
                messagebox.showerror("输出", f"无法创建输出目录:{e}", parent=self.root)
                return
        self.out_dir = out
        missing = self._collect_boxes()
        if missing:
            answer = messagebox.askyesno(
                "位点缺失",
                f"以下受体无共晶配体且未手动定义盒子:\n{', '.join(missing)}\n"
                "选择“是”=盲对接(探索性);“否”=回到 ②结合位点 手动定义。",
                parent=self.root)
            if answer:
                from dockstudio.core import box as boxmod
                for rec in self.receptors:
                    if rec["name"] in missing:
                        rp = st.parse_structure(rec["path"],
                                                keep_chains=rec.get("chains") or None)
                        self.boxes[rec["name"]] = boxmod.box_blind(rp).to_dict()
            else:
                self._goto(1)
                return
        cfg = self._build_cfg()
        if cfg.run_accuracy_report and not cfg.run_selfdock:
            messagebox.showinfo(
                "准确性评估",
                "已勾选“导出对接准确性评估报告”但未勾选“回贴验证”。准确性评估依赖自对接结果;"
                "若继续,报告将如实标注“不可评估”。", parent=self.root)
        if cfg.run_enrichment:
            bad = [k for k, v in (("actives", cfg.actives_path), ("decoys", cfg.decoys_path))
                   if not v or not os.path.isfile(v)]
            if bad:
                messagebox.showwarning(
                    "富集度验证",
                    "已启用富集度验证但缺少文件: " + ", ".join(bad) +
                    "\n本次将如实跳过 ROC/AUC(不会生成伪造结果)。可返回 ③参数 补齐后重跑。",
                    parent=self.root)
        self._stop_event.clear()
        phases = pipeline.PHASES if cfg.run_refine else \
            [p for p in pipeline.PHASES if p != "refine"]
        self._log("=== 启动 DockStudio 流水线 ===")
        self._log(f"受体: {cfg.receptor_names}  配体文件: {len(cfg.ligands)}  输出: {cfg.out_dir}")
        for r, b in cfg.boxes.items():
            self._log(f"盒子 {r}: {b.get('method')} center={b.get('center')} size={b.get('size')}")
        self._log(f"exhaustiveness={cfg.exhaustiveness} refine={cfg.refine_exhaustiveness} "
                  f"TopK={cfg.top_k} cpu={cfg.cpu} accuracy_report={cfg.run_accuracy_report}")
        self.progress["maximum"] = 100
        self.progress["value"] = 0
        self.status_var.set("运行中…")
        if self._foot_status:
            self._foot_status.set("运行中…")

        def worker():
            try:
                def prog(d):
                    self._q.put(("prog", d))
                summary = pipeline.run_project(cfg, on_log=self._log, on_progress=prog,
                                               phases=phases, time_slice_s=900,
                                               on_cancel=lambda: self._stop_event.is_set())
                self._q.put(("done", summary))
            except Exception as e:
                import traceback
                self._q.put(("error", traceback.format_exc()))

        self._worker = threading.Thread(target=worker, daemon=True)
        self._worker.start()

    def _request_stop(self):
        self._stop_event.set()
        self._log("停止请求已记录(将在当前阶段结束后停止)。")
        self.status_var.set("停止中…")

    def _poll_queue(self):
        try:
            while True:
                kind, payload = self._q.get_nowait()
                if kind == "log":
                    self._log_text_append(payload)
                elif kind == "prog":
                    d = payload
                    done = d.get("done")
                    total = d.get("total")
                    if total:
                        self.progress["value"] = int(100 * done / max(total, 1))
                    self.status_var.set(f"{d.get('stage', '')}  {d.get('current', '')} "
                                        f"(错误 {d.get('errors', 0)})")
                elif kind == "done":
                    s = payload
                    self.status_var.set("完成" if not s.get("interrupted")
                                        else "阶段中断(可再次点开始续跑)")
                    self._log("=== 流水线结束 ===")
                    self.progress["value"] = 100
                    if self._foot_status:
                        self._foot_status.set("完成" if not s.get("interrupted") else "已中断")
                elif kind == "error":
                    self.status_var.set("错误")
                    self._log_text_append("\n" + payload)
                    if self._foot_status:
                        self._foot_status.set("错误")
        except queue.Empty:
            pass
        self.root.after(120, self._poll_queue)

    def _log_text_append(self, msg):
        if not isinstance(msg, str):
            try:
                msg = str(msg)
            except Exception:
                msg = repr(msg)
        self.log_text.config(state="normal")
        self.log_text.insert("end", msg + "\n")
        self.log_text.see("end")
        self.log_text.config(state="disabled")

    # ------------------------------------------------------------- tab 5
    def _build_results_tab(self):
        f = self.tab_results
        top = ttk.Frame(f)
        top.pack(fill="x", pady=2)
        self._btn(top, "刷新结果", self._refresh_results, bootstyle="secondary-outline").pack(side="left")
        self._btn(top, "打开报告目录", lambda: self._open_subdir("reports")).pack(side="left", padx=6)
        self._btn(top, "打开图件目录", lambda: self._open_subdir("results")).pack(side="left", padx=6)
        self._btn(top, "打开准确性评估报告", self._open_accuracy_report,
                  bootstyle="primary").pack(side="left", padx=6)
        self._btn(top, "HTML 总报告", self._open_html_report,
                  bootstyle="success").pack(side="left", padx=6)
        self._btn(top, "3D 查看器目录", lambda: self._open_subdir("html_viewers"),
                  bootstyle="secondary-outline").pack(side="left", padx=6)
        self._btn(top, "打包导出 .zip", self._export_zip,
                  bootstyle="info").pack(side="left", padx=6)

        self.res_tree = ttk.Treeview(f, columns=("receptor", "aff", "src"),
                                     show="headings", height=11)
        self.res_tree.heading("receptor", text="受体")
        self.res_tree.heading("aff", text="结合能 (kcal/mol)")
        self.res_tree.heading("src", text="来源")
        self.res_tree.column("receptor", width=220)
        self.res_tree.column("aff", width=140)
        self.res_tree.column("src", width=160)
        self.res_tree.pack(fill="both", expand=True, pady=6)

        self.res_note = tk.StringVar(value="")
        ttk.Label(f, textvariable=self.res_note, anchor="w", foreground="#666").pack(fill="x")

    def _open_subdir(self, sub):
        base = self.out_var.get()
        if base and base != OUT_FOLDER_MARKER:
            d = os.path.join(base, sub)
            if os.path.isdir(d):
                if _open_path(d):
                    return
                messagebox.showinfo("路径", d, parent=self.root)
                return
        messagebox.showinfo("提示", "尚无结果。", parent=self.root)

    def _open_accuracy_report(self):
        base = self.out_var.get()
        if base and base != OUT_FOLDER_MARKER:
            p = os.path.join(base, "reports", "03_accuracy_assessment.md")
            if os.path.exists(p):
                if _open_path(p):
                    return
                messagebox.showinfo("路径", p, parent=self.root)
                return
        messagebox.showinfo("提示", "尚无准确性评估报告。\n请先在 ③参数 勾选“导出对接准确性评估报告”"
                                   "并运行流水线(需自对接启用)。", parent=self.root)

    def _open_html_report(self):
        base = self.out_var.get()
        if base and base != OUT_FOLDER_MARKER:
            p = os.path.join(base, "reports", "html", "index.html")
            if os.path.exists(p):
                if _open_path(p):
                    return
                messagebox.showinfo("路径", p, parent=self.root)
                return
        messagebox.showinfo("提示", "尚无 HTML 总报告。\n请运行 v2.0 流水线(报告/HTML 阶段)后刷新。",
                            parent=self.root)

    def _export_zip(self):
        from dockstudio.core import project
        base = self.out_var.get()
        if not base or base == OUT_FOLDER_MARKER or not os.path.isdir(base):
            messagebox.showinfo("提示", "请先选择/确认输出目录。", parent=self.root)
            return
        try:
            dest = project.export_run_zip(base)
        except Exception as e:
            messagebox.showerror("导出失败", str(e), parent=self.root)
            return
        if os.path.exists(dest):
            self._q.put(("log", f"已打包导出: {dest}"))
            messagebox.showinfo("导出完成",
                                f"可复现数据包已生成(含全部源文件):\n{dest}", parent=self.root)
        else:
            messagebox.showerror("导出失败", "未生成 zip。", parent=self.root)

    def _refresh_results(self):
        base = self.out_var.get()
        self.res_tree.delete(*self.res_tree.get_children())
        if not base or base == OUT_FOLDER_MARKER:
            self.res_note.set("未选择输出目录")
            return
        k = self.topk.get() if hasattr(self, "topk") else 5
        fp = os.path.join(base, "results", f"final_top{k}.csv")
        if not os.path.exists(fp):
            # fall back to any final_top*.csv (e.g. an older v1.0 layout)
            cands = sorted(f for f in os.listdir(os.path.join(base, "results"))
                           if f.startswith("final_top") and f.endswith(".csv")) \
                if os.path.isdir(os.path.join(base, "results")) else []
            fp = os.path.join(base, "results", cands[-1]) if cands else ""
        if not fp or not os.path.exists(fp):
            self.res_note.set("尚无 final_top*.csv(可能未完成运行)")
            return
        import csv
        with open(fp, encoding="utf-8-sig", newline="") as fh:
            for row in csv.DictReader(fh):
                self.res_tree.insert("", "end",
                                     values=(row.get("receptor"),
                                             row.get("best_affinity_kcal_mol"),
                                             row.get("source")))
        qp = os.path.join(base, "reports", "02_quality_assessment.md")
        ap = os.path.join(base, "reports", "03_accuracy_assessment.md")
        note = "最终 Top-K 见上。"
        if os.path.exists(qp):
            note += " QC 报告: 已生成."
        else:
            note += " QC 报告: 未生成."
        if os.path.exists(ap):
            note += " 对接准确性评估: 已生成."
        self.res_note.set(note)

    # ------------------------------------------------------------- Help viewer
    def _open_help(self):
        HelpViewer(self.root, version=__version__)


class HelpViewer:
    """A simple in-app documentation viewer (v1.1 feature)."""

    # topics: display title -> (relative path from project root, description)
    TOPICS = [
        ("用户手册(快速开始)", "用户手册.md"),
        ("架构设计", os.path.join("docs", "架构设计.md")),
        ("定版记录 v2.0.0", os.path.join("docs", "定版记录_v2.0.0.md")),
        ("定版记录 v1.1.0", os.path.join("docs", "定版记录_v1.1.0.md")),
        ("打包说明", os.path.join("packaging", "打包说明.md")),
        ("README", "README.md"),
        ("第三方声明", "THIRD_PARTY_NOTICES.md"),
    ]

    # shown when the doc file is not available (e.g. a minimal packaged build)
    BUILTIN_HELP = (
        "DockStudio 帮助中心\n"
        "===================\n"
        "DockStudio 是一站式、可复现的批量分子对接平台(仅图形界面)。\n\n"
        "操作流程\n"
        "1. 输入:添加受体(PDB .pdb/.ent 或 mmCIF .cif/.mmcif)与配体库(SDF/MOL),"
        "并用“浏览文件夹…”选择输出目录。\n"
        "2. 结合位点:选中受体 → 自动(共晶配体质心) / 手动(中心+尺寸) / 残基 / 盲对接;"
        "点“预览盒子”确认。\n"
        "3. 参数:按需调整 exhaustiveness、Top-K、CPU 与目标 pH;可勾选导出"
        "“对接准确性评估报告”。\n"
        "4. 运行:点“开始 / 断点续跑”;可随时停止,断点自动续跑。\n"
        "5. 结果:查看 Top-K、报告目录、图件与准确性评估报告。\n\n"
        "科学严谨性\n"
        "• 所有数值来自真实运行;报告只记录真实执行的步骤。\n"
        "• 无共晶配体受体自动标注探索性;质子化为文档化简单规则(非 pKa 预测)。\n"
        "• 发布前请人工复核图件、质子化态与关键结合模式。\n\n"
        "(详细说明文件未随当前安装包含;完整文档见软件发行包 docs/ 目录。)"
    )

    def __init__(self, master, version=""):
        top = tk.Toplevel(master)
        top.title("DockStudio 帮助中心")
        top.geometry("900x640")
        top.transient(master)
        self.top = top
        self.version = version

        bar = ttk.Frame(top, padding=(10, 8))
        bar.pack(fill="x")
        ttk.Label(bar, text="帮助中心", font=("Microsoft YaHei UI", 14, "bold")).pack(side="left")
        ttk.Label(bar, text=f"DockStudio v{version}", foreground="#888").pack(side="right")

        body = ttk.Frame(top)
        body.pack(fill="both", expand=True, padx=10, pady=6)
        body.columnconfigure(1, weight=1)
        body.rowconfigure(0, weight=1)

        left = ttk.Frame(body)
        left.grid(row=0, column=0, sticky="ns", padx=(0, 8))
        ttk.Label(left, text="主题", foreground="#888").pack(anchor="w")
        self._list = tk.Listbox(left, width=30, height=18, exportselection=False)
        for title, _path in self.TOPICS:
            self._list.insert("end", title)
        self._list.pack(fill="y", expand=True, pady=4)
        self._list.bind("<<ListboxSelect>>", lambda e: self._show())

        self.text = tk.Text(body, wrap="word", font=("Microsoft YaHei UI", 10),
                            relief="flat", state="disabled", bg="#fbfcfd")
        self.text.grid(row=0, column=1, sticky="nsew")
        if HAVE_TTB:
            self._btn = ttk.Button(bar, text="在外部打开…", bootstyle="secondary",
                                   command=self._open_external)
        else:
            self._btn = ttk.Button(bar, text="在外部打开…", command=self._open_external)
        self._btn.pack(side="right")

        self._selected = 0
        self._list.selection_set(0)
        self._show()

    def _topic_path(self, idx: int) -> str:
        here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # dockstudio/
        root = os.path.dirname(here)
        rel = self.TOPICS[idx][1]
        return os.path.join(root, rel)

    def _show(self):
        sel = self._list.curselection()
        if not sel:
            return
        idx = int(sel[0])
        self._selected = idx
        path = self._topic_path(idx)
        content = ""
        if os.path.exists(path):
            try:
                with open(path, encoding="utf-8") as fh:
                    content = fh.read()
            except Exception as e:
                content = f"[无法读取 {path}]\n{e}"
        elif idx == 0:
            content = self.BUILTIN_HELP
        else:
            content = f"[文档未随安装附带: {path}]\n\n请在源码发行包的 docs/ 与 packaging/ 目录中查看。"
        self.text.config(state="normal")
        self.text.delete("1.0", "end")
        self.text.insert("1.0", content)
        self.text.config(state="disabled")

    def _open_external(self):
        path = self._topic_path(self._selected)
        if os.path.exists(path):
            if not _open_path(path):
                messagebox.showinfo("路径", path, parent=self.top)
        else:
            messagebox.showinfo("提示", "该文档未随当前安装附带,无法在外部打开。",
                                parent=self.top)


def _make_root():
    """Create the themed (or fallback) root window."""
    if HAVE_TTB:
        try:
            return ttkb.Window(themename="litera", title=f"{APP_TITLE}  v{__version__}")
        except Exception:
            return ttkb.Window(themename="cosmo", title=f"{APP_TITLE}  v{__version__}")
    return tk.Tk()


def main(argv=None) -> int:
    if not _TK_OK:
        # Very defensive: no tkinter present
        raise SystemExit("DockStudio GUI requires a Tk-enabled Python (tkinter).")
    try:
        import multiprocessing
        multiprocessing.freeze_support()
    except Exception:
        pass
    root = _make_root()
    try:
        if HAVE_TTB:
            # ttkbootstrap theme already chosen by Window(themename=...).
            _ = ttkb.Style()
        else:
            style = ttk.Style(root)
            if "clam" in style.theme_names():
                style.theme_use("clam")
    except Exception:
        pass
    DockStudioApp(root)
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
