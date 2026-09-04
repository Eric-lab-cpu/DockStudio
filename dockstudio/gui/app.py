"""DockStudio 中文图形界面 (tkinter).

入口:
    python -m dockstudio
或者打包后的 DockStudio.exe.

界面步骤:
  1) 输入 : 选择受体 PDB 与配体库 SDF + 输出目录 + 链/共晶配体
  2) 位点 : 每个受体定义搜索盒子(自动共晶 / 手动中心尺寸 / 残基 / 盲对接)
  3) 参数 : 打分与计算预算
  4) 运行 : 启动 / 断点续跑 / 实时日志
  5) 结果 : 查看 Top-K 表、打开报告/图片目录

引擎在后台线程运行;界面通过 queue 刷新进度与日志。
"""

from __future__ import annotations

import os
import queue
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from tkinter.scrolledtext import ScrolledText

from dockstudio._version import APP_NAME, APP_TITLE, BRAND, COPYRIGHT_CN, __version__
from dockstudio.core import envinfo, inventory, models, pipeline, structure as st

OUT_FOLDER_MARKER = "输出目录(结果将写入其中)"


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


def _scan_polymer_chains(pdb_path: str) -> list:
    try:
        inv = st.inventory_pdb(pdb_path)
        return [c for c, ci in inv.chains.items() if ci.polymer_res > 0]
    except Exception:
        return []


def _find_cocrystal(pdb_path: str, chains: list) -> list:
    try:
        rows = inventory.inventory_receptor(pdb_path, _stem(pdb_path), chains)
        found = []
        for r in rows:
            if r.get("ligand") and r["chain"] in chains:
                found.append(f"{r['ligand']} ({r['chain']}{r.get('ligand_resseq')}, "
                             f"{r.get('ligand_n_heavy')} heavy)")
        # dedupe preserving order
        seen = set(); out = []
        for x in found:
            if x not in seen:
                seen.add(x); out.append(x)
        return out
    except Exception:
        return []


class DockStudioApp:
    def __init__(self, root: tk.Tk):
        self.root = root
        root.title(f"{APP_TITLE}  v{__version__}")
        root.geometry("1080x800")
        self._icon = None
        self._apply_icon()
        self._q: queue.Queue = queue.Queue()
        self._worker: threading.Thread | None = None
        self._stop_event = threading.Event()

        # ---- data model (equivalent of RunConfig built lazily) ----
        self.receptors = []      # list[dict]: path,name,chains(list)
        self.ligands = []        # list[dict]: path,name
        self.out_dir = ""
        self.boxes = {}          # rec name -> dict(method,center,size,basis,ligand_ref)

        self._build_menu()
        self._build_ui()
        self._build_footer()
        self.root.after(120, self._poll_queue)
        self.nb.bind("<<NotebookTabChanged>>", self._on_tab_changed)
        # environment check note (only block on required tools)
        env = envinfo.gather_environment()
        missing = env.get("missing_required") or []
        miss_opt = env.get("missing_optional") or []
        if missing:
            messagebox.showwarning("依赖提示",
                                   "以下必需组件未检测到,相关功能将不可用:\n" + ", ".join(missing) +
                                   "\n\n请按 打包说明/用户手册 安装依赖。", parent=self.root)
        elif miss_opt:
            note = "可选组件未检测到(仅降级,不影响主流程): " + ", ".join(miss_opt)
            self._q.put(("log", note))

    # ------------------------------------------------------------- helpers
    def _resource_path(self, name: str) -> str:
        """Return an absolute path to a packaged resource file."""
        import sys
        # PyInstaller one-file runtime
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

    # ------------------------------------------------------------- menu
    def _build_menu(self):
        menubar = tk.Menu(self.root)
        m_file = tk.Menu(menubar, tearoff=0)
        m_file.add_command(label="退出", command=self.root.destroy)
        menubar.add_cascade(label="文件", menu=m_file)
        m_help = tk.Menu(menubar, tearoff=0)
        m_help.add_command(label="关于 " + APP_NAME, command=self._about)
        menubar.add_cascade(label="帮助", menu=m_help)
        self.root.config(menu=menubar)

    def _about(self):
        top = tk.Toplevel(self.root)
        top.title("关于")
        top.geometry("480x300")
        top.transient(self.root)
        top.resizable(False, False)
        box = ttk.Frame(top, padding=18)
        box.pack(fill="both", expand=True)
        ttk.Label(box, text=APP_NAME, font=("Microsoft YaHei UI", 15, "bold")).pack(anchor="center")
        ttk.Label(box, text=f"版本 {__version__}", font=("", 10)).pack(pady=(4, 2))
        ttk.Separator(box).pack(fill="x", pady=8)
        ttk.Label(box, text="一站式自动化批量分子对接平台\n"
                            "AutoDock Vina · Meeko · RDKit · PLIP · PyMOL",
                  justify="center").pack(pady=4)
        ttk.Label(box, text=f"本产品由 {BRAND} 提供", font=("", 11, "bold")).pack(pady=(12, 2))
        ttk.Label(box, text=COPYRIGHT_CN, foreground="#666").pack()
        ttk.Button(top, text="关闭", command=top.destroy).pack(pady=10)

    # ------------------------------------------------------------------ UI
    def _build_ui(self):
        nb = ttk.Notebook(self.root)
        nb.pack(fill="both", expand=True, padx=6, pady=6)
        self.nb = nb
        self.tab_inputs = ttk.Frame(nb)
        self.tab_site = ttk.Frame(nb)
        self.tab_params = ttk.Frame(nb)
        self.tab_run = ttk.Frame(nb)
        self.tab_results = ttk.Frame(nb)
        nb.add(self.tab_inputs, text="① 输入")
        nb.add(self.tab_site, text="② 结合位点")
        nb.add(self.tab_params, text="③ 参数")
        nb.add(self.tab_run, text="④ 运行")
        nb.add(self.tab_results, text="⑤ 结果")
        self._build_inputs_tab()
        self._build_site_tab()
        self._build_params_tab()
        self._build_run_tab()
        self._build_results_tab()

    def _build_footer(self):
        bar = ttk.Frame(self.root)
        bar.pack(fill="x", padx=8, pady=(0, 4))
        ttk.Label(bar, text=f"{BRAND} · 版本 {__version__} · 仅供科研使用",
                  foreground="#777").pack(side="left")
        ttk.Label(bar, text=COPYRIGHT_CN, foreground="#aaa").pack(side="right")

    def _on_tab_changed(self, event=None):
        try:
            idx = self.nb.index(self.nb.select())
        except Exception:
            return
        if idx == 1:  # ② 结合位点
            self._refresh_site_lists()

    # ------------------------------------------------------------- tab 1
    def _build_inputs_tab(self):
        f = self.tab_inputs
        left = ttk.LabelFrame(f, text="受体(PDB) - 可添加多个")
        left.pack(side="left", fill="both", expand=True, padx=5, pady=5)
        right = ttk.LabelFrame(f, text="配体库(SDF) - 可添加多个文件")
        right.pack(side="left", fill="both", expand=True, padx=5, pady=5)

        self.rec_list = ttk.Treeview(left, columns=("chains",), height=9)
        self.rec_list.heading("#0", text="受体文件")
        self.rec_list.heading("chains", text="对接链")
        self.rec_list.column("#0", width=330)
        self.rec_list.column("chains", width=120)
        self.rec_list.pack(fill="both", expand=True, padx=4, pady=4)
        bf = ttk.Frame(left); bf.pack(fill="x", padx=4, pady=2)
        ttk.Button(bf, text="添加 PDB…", command=self._add_receptor).pack(side="left")
        ttk.Button(bf, text="删除选中", command=self._del_receptor).pack(side="left", padx=4)
        ttk.Button(bf, text="上移/下移", command=lambda: self._move_receptor(-1)).pack(side="left")
        ttk.Button(bf, text="↓", command=lambda: self._move_receptor(1)).pack(side="left")
        ttk.Button(bf, text="链/共晶…", command=self._edit_receptor).pack(side="left", padx=4)

        self.lig_list = ttk.Treeview(right, columns=("name",), height=9)
        self.lig_list.heading("#0", text="配体库文件")
        self.lig_list.heading("name", text="库名")
        self.lig_list.column("#0", width=300)
        self.lig_list.column("name", width=120)
        self.lig_list.pack(fill="both", expand=True, padx=4, pady=4)
        lf = ttk.Frame(right); lf.pack(fill="x", padx=4, pady=2)
        ttk.Button(lf, text="添加 SDF…", command=self._add_ligand).pack(side="left")
        ttk.Button(lf, text="删除选中", command=self._del_ligand).pack(side="left", padx=4)

        # output dir
        od = ttk.LabelFrame(f, text="输出目录")
        od.pack(side="top", fill="x", padx=5, pady=3)
        row = ttk.Frame(od); row.pack(fill="x", padx=4, pady=4)
        self.out_var = tk.StringVar(value=OUT_FOLDER_MARKER)
        ttk.Entry(row, textvariable=self.out_var).pack(side="left", fill="x", expand=True)
        ttk.Button(row, text="浏览…", command=self._choose_out).pack(side="left", padx=4)
        ttk.Button(od, text="示例: 使用 4DFR 演示数据", command=self._load_demo).pack(anchor="w", padx=4, pady=2)

    def _add_receptor(self):
        paths = filedialog.askopenfilenames(parent=self.root, title="选择受体 PDB",
                                            filetypes=[("PDB", "*.pdb"), ("所有文件", "*.*")])
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
        i = int(sel[0]); j = i + d
        if 0 <= j < len(self.receptors):
            self.receptors[i], self.receptors[j] = self.receptors[j], self.receptors[i]
            self._refresh_rec_list(); self.rec_list.selection_set(str(j))

    def _edit_receptor(self):
        sel = self.rec_list.selection()
        if not sel:
            return
        i = int(sel[0])
        rec = self.receptors[i]
        dlg = tk.Toplevel(self.root); dlg.title(f"受体设置 - {rec['name']}")
        dlg.geometry("440x340"); dlg.transient(self.root)
        ttk.Label(dlg, text="对接链(可多选, Ctrl/Shift):").pack(anchor="w", padx=8, pady=4)
        lb = tk.Listbox(dlg, selectmode="extended", height=10)
        all_chains = _scan_polymer_chains(rec["path"]) or (rec["chains"] or ["A"])
        for c in all_chains:
            lb.insert("end", c)
        for j, c in enumerate(all_chains):
            if c in (rec.get("chains") or []):
                lb.selection_set(j)
        lb.pack(fill="both", expand=True, padx=8)
        coc = _find_cocrystal(rec["path"], all_chains)
        info = ("检测到共晶配体:\n" + "\n".join("  - " + c for c in coc)) if coc else "未在所选链中检测到共晶配体。"
        ttk.Label(dlg, text=info, justify="left").pack(anchor="w", padx=8, pady=4)

        def apply():
            rec["chains"] = [all_chains[x] for x in lb.curselection()] or all_chains
            # default box auto from cocrystal
            if rec["name"] in self.boxes:
                self.boxes[rec["name"]]["basis"] = "auto (edited chains)"
            self._refresh_rec_list()
            dlg.destroy()

        ttk.Button(dlg, text="确定", command=apply).pack(pady=6)

    def _add_ligand(self):
        paths = filedialog.askopenfilenames(parent=self.root, title="选择配体库 SDF",
                                            filetypes=[("SDF", "*.sdf *.mol")])
        for p in paths:
            if any(l["path"] == p for l in self.ligands):
                continue
            self.ligands.append({"path": p, "name": _stem(p)})
        self._refresh_lig_list()

    def _del_ligand(self):
        sel = self.lig_list.selection()
        if sel:
            self.ligands.pop(int(sel[0])); self._refresh_lig_list()

    def _choose_out(self):
        p = filedialog.askdirectory(parent=self.root, title="选择输出目录")
        if p:
            self.out_dir = p; self.out_var.set(p)

    def _load_demo(self):
        demo = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                            "examples", "demo")
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
        self._refresh_rec_list(); self._refresh_lig_list()
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
        ttk.Label(f, text="每个受体的搜索盒子。选择受体后,根据共晶配体/残基自动或手动定义。",
                  wraplength=1000).pack(anchor="w", padx=8, pady=4)
        top = ttk.Frame(f); top.pack(fill="x", padx=8)
        self.site_rec_var = tk.StringVar()
        self.site_combo = ttk.Combobox(top, textvariable=self.site_rec_var, width=40, state="readonly")
        self.site_combo.pack(side="left")
        self.site_combo.bind("<<ComboboxSelected>>", lambda e: self._load_site_for_rec())
        ttk.Button(top, text="刷新(重读 PDB)", command=self._refresh_site_lists).pack(side="left", padx=6)
        ttk.Button(top, text="预览盒子(在结果中显示)", command=self._preview_box).pack(side="left")

        body = ttk.Frame(f); body.pack(fill="both", expand=True, padx=8, pady=6)
        self.site_info = tk.Text(body, height=6, state="disabled", bg="#f5f5f5")
        self.site_info.pack(fill="x")
        stylef = ttk.LabelFrame(body, text="盒子定义")
        stylef.pack(fill="x", pady=6)
        self.box_mode = tk.StringVar(value="auto")
        ttk.Radiobutton(stylef, text="自动: 共晶配体质心 (推荐,仅当存在共晶配体)", value="auto",
                        variable=self.box_mode).pack(anchor="w")
        ttk.Radiobutton(stylef, text="手动: 中心 + 尺寸", value="explicit",
                        variable=self.box_mode).pack(anchor="w")
        ttk.Radiobutton(stylef, text="残基中心(列出如 TRP168 或以逗号分隔编号)", value="residues",
                        variable=self.box_mode).pack(anchor="w")
        ttk.Radiobutton(stylef, text="盲对接(全蛋白,探索性,高假阳性)", value="blind",
                        variable=self.box_mode).pack(anchor="w")
        grid = ttk.Frame(body); grid.pack(fill="x", padx=8)
        self.cx = tk.DoubleVar(value=0.0); self.cy = tk.DoubleVar(value=0.0); self.cz = tk.DoubleVar(value=0.0)
        self.sx = tk.DoubleVar(value=24.0); self.sy = tk.DoubleVar(value=24.0); self.sz = tk.DoubleVar(value=24.0)
        for i, (lab, var) in enumerate([("中心 X", self.cx), ("中心 Y", self.cy), ("中心 Z", self.cz)]):
            ttk.Label(grid, text=lab).grid(row=0, column=i * 2, padx=2, pady=2)
            ttk.Entry(grid, textvariable=var, width=9).grid(row=0, column=i * 2 + 1, padx=2)
        for i, (lab, var) in enumerate([("尺寸 X", self.sx), ("尺寸 Y", self.sy), ("尺寸 Z", self.sz)]):
            ttk.Label(grid, text=lab).grid(row=1, column=i * 2, padx=2, pady=2)
            ttk.Entry(grid, textvariable=var, width=9).grid(row=1, column=i * 2 + 1, padx=2)
        self.site_res_var = tk.StringVar(value="")
        rr = ttk.Frame(body); rr.pack(fill="x", padx=8, pady=4)
        ttk.Label(rr, text="残基列表:").pack(side="left")
        ttk.Entry(rr, textvariable=self.site_res_var, width=50).pack(side="left", padx=4)

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
        info = []
        chains = rec.get("chains") or []
        info.append(f"受体: {name}  对接链: {','.join(chains)}")
        coc = _find_cocrystal(rec["path"], chains)
        info.append("共晶配体: " + ("; ".join(coc) if coc else "无(请手动定义,或选盲对接)"))
        if box:
            info.append(f"现有盒子: method={box.get('method')} center={box.get('center')} "
                        f"size={box.get('size')}")
        self.site_info.config(state="normal")
        self.site_info.delete("1.0", "end")
        self.site_info.insert("1.0", "\n".join(info))
        self.site_info.config(state="disabled")

    def _preview_box(self):
        name = self.site_rec_var.get()
        if not name:
            messagebox.showinfo("提示", "请先在 ①输入 添加受体。", parent=self.root); return
        rec = next((r for r in self.receptors if r["name"] == name), None)
        if not rec:
            return
        mode = self.box_mode.get()
        box = self._compute_box_for(rec, mode, silent=False)
        if box is None:
            return
        self.boxes[name] = box
        self._load_site_for_rec()
        messagebox.showinfo("盒子", f"{name}\n方法: {box['method']}\n"
                           f"中心: {[round(x,2) for x in box['center']]}\n"
                           f"尺寸: {[round(x,2) for x in box['size']]}\n{box.get('basis','')}",
                           parent=self.root)

    def _compute_box_for(self, rec, mode, silent=True):
        from dockstudio.core import box as boxmod
        pdb = rec["path"]; chains = rec.get("chains") or []
        if mode == "auto":
            rows = inventory.inventory_receptor(pdb, rec["name"], chains)
            cand = None
            for r in rows:
                if r.get("ligand") and r["chain"] in chains:
                    if cand is None or (r.get("ligand_n_heavy") or 0) > (cand.get("ligand_n_heavy") or 0):
                        cand = r
            if cand:
                ligr = st.coords_of_ligand(pdb, cand["ligand"], cand["chain"], cand.get("ligand_resseq"))
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
                messagebox.showwarning("位点", "请填写残基列表", parent=self.root); return None
            tokens = [t.strip() for t in txt.replace(";", ",").split(",")]
            recs = st.parse_pdb(pdb, keep_chains=chains or None)
            try:
                b = boxmod.box_from_residues(recs, tokens)
                return b.to_dict()
            except Exception as e:
                messagebox.showerror("残基", str(e), parent=self.root)
                return None
        if mode == "blind":
            recs = st.parse_pdb(pdb, keep_chains=chains or None)
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
                    if cand is None or (r.get("ligand_n_heavy") or 0) > (cand.get("ligand_n_heavy") or 0):
                        cand = r
            if cand:
                ligr = st.coords_of_ligand(rec["path"], cand["ligand"], cand["chain"], cand.get("ligand_resseq"))
                self.boxes[name] = boxmod.box_from_ligand(ligr, ligand_tag=cand["ligand"]).to_dict()
            else:
                missing.append(name)
        return missing

    # ------------------------------------------------------------- tab 3
    def _build_params_tab(self):
        f = self.tab_params
        p = ttk.LabelFrame(f, text="计算参数(与规格书默认一致,可下调以省时)")
        p.pack(fill="x", padx=8, pady=6)
        grid = ttk.Frame(p); grid.pack(padx=8, pady=4)
        self.exh = tk.IntVar(value=16)
        self.exh_refine = tk.IntVar(value=32)
        self.nposes = tk.IntVar(value=9)
        self.erange = tk.DoubleVar(value=4.0)
        self.top_ref = tk.IntVar(value=12)
        self.topk = tk.IntVar(value=5)
        self.cpu = tk.IntVar(value=2)
        self.ph = tk.DoubleVar(value=7.4)
        fields = [("初筛 exhaustiveness", self.exh, 16), ("精修 exhaustiveness", self.exh_refine, 32),
                  ("每配体 pose 数", self.nposes, 9), ("energy_range (kcal/mol)", self.erange, 4.0),
                  ("精修短名单", self.top_ref, 12), ("最终 Top-K", self.topk, 5),
                  ("并行 CPU", self.cpu, 2), ("pH", self.ph, 7.4)]
        for i, (lab, var, _def) in enumerate(fields):
            r, c = divmod(i, 2)
            ttk.Label(grid, text=lab).grid(row=r, column=c * 2, sticky="e", padx=3, pady=3)
            ttk.Entry(grid, textvariable=var, width=8).grid(row=r, column=c * 2 + 1, padx=3)
        toggles = ttk.LabelFrame(f, text="开关")
        toggles.pack(fill="x", padx=8, pady=6)
        self.chk_self = tk.BooleanVar(value=True)
        self.chk_plip = tk.BooleanVar(value=True)
        self.chk_viz = tk.BooleanVar(value=True)
        self.chk_refine = tk.BooleanVar(value=True)
        self.chk_pse = tk.BooleanVar(value=True)
        for i, (var, txt) in enumerate([(self.chk_self, "回贴验证(self-docking)"),
                                        (self.chk_plip, "相互作用分析 PLIP"),
                                        (self.chk_viz, "PyMOL 图件 + .pml"),
                                        (self.chk_refine, "精修(Top-K 前)"),
                                        (self.chk_pse, "输出 .pse 会话")]):
            ttk.Checkbutton(toggles, text=txt, variable=var).grid(row=0, column=i, padx=8)

    def _build_cfg(self) -> models.RunConfig:
        missing = self._collect_boxes()
        recs = [{"path": r["path"], "name": r["name"], "chains": r.get("chains") or []}
                for r in self.receptors]
        ligs = [{"path": l["path"], "name": l["name"]} for l in self.ligands]
        return models.RunConfig(
            receptors=recs, ligands=ligs, out_dir=self.out_dir, title="DockStudio GUI run",
            boxes=self.boxes,
            exhaustiveness=self.exh.get(), refine_exhaustiveness=self.exh_refine.get(),
            n_poses=self.nposes.get(), energy_range=self.erange.get(),
            top_refine=self.top_ref.get(), top_k=self.topk.get(),
            cpu=self.cpu.get(), ph=self.ph.get(),
            run_selfdock=self.chk_self.get(), run_plip=self.chk_plip.get(),
            run_visualization=self.chk_viz.get(), run_refine=self.chk_refine.get(),
            write_pse=self.chk_pse.get(), delete_bad_res=True, overwrite=False)

    # ------------------------------------------------------------- tab 4
    def _build_run_tab(self):
        f = self.tab_run
        row = ttk.Frame(f); row.pack(fill="x", padx=8, pady=6)
        self.start_btn = ttk.Button(row, text="开始 / 断点续跑", command=self._start_run)
        self.start_btn.pack(side="left")
        ttk.Button(row, text="停止(阶段结束后)", command=self._request_stop).pack(side="left", padx=6)
        ttk.Button(row, text="打开输出目录", command=self._open_out).pack(side="left", padx=6)
        self.progress = ttk.Progressbar(f, mode="determinate")
        self.progress.pack(fill="x", padx=8)
        self.log_text = ScrolledText(f, height=18, state="disabled", font=("Consolas", 9))
        self.log_text.pack(fill="both", expand=True, padx=8, pady=6)
        self.status_var = tk.StringVar(value="就绪")
        ttk.Label(f, textvariable=self.status_var, anchor="w").pack(fill="x", padx=8)

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
            messagebox.showwarning("输入", "请先添加受体与配体库。", parent=self.root); return
        out = self.out_var.get()
        if not out or out == OUT_FOLDER_MARKER:
            messagebox.showwarning("输出", "请选择输出目录。", parent=self.root); return
        self.out_dir = out
        missing = self._collect_boxes()
        if missing:
            answer = messagebox.askyesno("位点缺失",
                                         f"以下受体无共晶配体且未手动定义盒子:\n{', '.join(missing)}\n"
                                         "选择‘是’=盲对接(探索性);‘否’=回到 ②结合位点 手动定义。",
                                         parent=self.root)
            if answer:
                from dockstudio.core import box as boxmod
                for rec in self.receptors:
                    if rec["name"] in missing:
                        rp = st.parse_pdb(rec["path"], keep_chains=rec.get("chains") or None)
                        self.boxes[rec["name"]] = boxmod.box_blind(rp).to_dict()
            else:
                self.nb.select(self.tab_site)
                return
        cfg = self._build_cfg()
        self._stop_event.clear()
        phases = pipeline.PHASES if cfg.run_refine else [p for p in pipeline.PHASES if p != "refine"]
        self._log("=== 启动 DockStudio 流水线 ===")
        self._log(f"受体: {cfg.receptor_names}  配体文件: {len(cfg.ligands)}  输出: {cfg.out_dir}")
        for r, b in cfg.boxes.items():
            self._log(f"盒子 {r}: {b.get('method')} center={b.get('center')} size={b.get('size')}")
        self._log(f"exhaustiveness={cfg.exhaustiveness} refine={cfg.refine_exhaustiveness} "
                  f"TopK={cfg.top_k} cpu={cfg.cpu}")
        self.progress["maximum"] = 100
        self.progress["value"] = 0
        self.status_var.set("运行中…")

        def worker():
            try:
                def logf(m):
                    self._log(m)
                def prog(d):
                    total = 100
                    # map coarse stage progress
                    self._q.put(("prog", d))
                summary = pipeline.run_project(cfg, on_log=logf, on_progress=prog,
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
                    # estimate progress: if remaining given, use fraction done
                    done = d.get("done"); total = d.get("total")
                    if total:
                        val = int(100 * done / max(total, 1))
                        self.progress["value"] = val
                    self.status_var.set(f"{d.get('stage','')}  {d.get('current','')} "
                                        f"(错误 {d.get('errors',0)})")
                elif kind == "done":
                    s = payload
                    self.status_var.set("完成" if not s.get("interrupted") else "阶段中断(可再次点开始续跑)")
                    self._log("=== 流水线结束 ===")
                    self.progress["value"] = 100
                elif kind == "error":
                    self.status_var.set("错误")
                    self._log_text_append("\n" + payload)
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
        row = ttk.Frame(f); row.pack(fill="x", padx=8, pady=4)
        ttk.Button(row, text="刷新结果", command=self._refresh_results).pack(side="left")
        ttk.Button(row, text="打开报告目录", command=lambda: self._open_subdir("reports")).pack(side="left", padx=6)
        ttk.Button(row, text="打开图件目录", command=lambda: self._open_subdir("results")).pack(side="left", padx=6)
        self.res_tree = ttk.Treeview(f, columns=("receptor", "aff", "src"), show="headings", height=12)
        self.res_tree.heading("receptor", text="受体")
        self.res_tree.heading("aff", text="结合能 (kcal/mol)")
        self.res_tree.heading("src", text="来源")
        self.res_tree.column("receptor", width=220)
        self.res_tree.column("aff", width=140)
        self.res_tree.column("src", width=160)
        self.res_tree.pack(fill="both", expand=True, padx=8, pady=6)
        self.res_note = tk.StringVar(value="")
        ttk.Label(f, textvariable=self.res_note, anchor="w").pack(fill="x", padx=8)

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

    def _refresh_results(self):
        base = self.out_var.get()
        self.res_tree.delete(*self.res_tree.get_children())
        if not base or base == OUT_FOLDER_MARKER:
            self.res_note.set("未选择输出目录")
            return
        fp = os.path.join(base, "results", "final_top5.csv")
        if not os.path.exists(fp):
            self.res_note.set("尚无 final_top5.csv(可能未完成运行)")
            return
        import csv
        with open(fp, encoding="utf-8-sig", newline="") as fh:
            for row in csv.DictReader(fh):
                self.res_tree.insert("", "end",
                                     values=(row.get("receptor"), row.get("best_affinity_kcal_mol"),
                                             row.get("source")))
        qp = os.path.join(base, "reports", "02_quality_assessment.md")
        self.res_note.set("最终 Top-K 见上。QC 报告: " + ("已生成" if os.path.exists(qp) else "未生成"))


def main(argv=None) -> int:
    root = tk.Tk()
    try:
        style = ttk.Style(root)
        if "vista" in style.theme_names():
            style.theme_use("vista")
    except Exception:
        pass
    DockStudioApp(root)
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
