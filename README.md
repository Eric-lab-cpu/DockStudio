# DockStudio — 一站式自动化批量分子对接平台

> **© 2026 Eric Studio. All rights reserved. / 版权所有 © 2026 Eric Studio。保留所有权利。**
> DockStudio 由 Eric Studio 提供 · **当前版本 v2.0.0**

**DockStudio** 是一个 **图形界面为主、引擎可编程** 的分子对接工具,内置一套严谨、
可复现、可发表级的自动化流程:受体/配体盘点 → 结构准备(Meeko)→ 搜索盒子决策 →
AutoDock Vina 批处理(断点续跑、进程池并行)→ 排名精修 Top-K → 共晶配体回贴验证
(**对称性感知 RMSD**)→ PLIP 相互作用分析 → PyMOL 论文级图件与可编辑
`.pml`/`.pse` → 方法与质量报告 →(可选)**对接准确性评估报告** →(可选)
**ROC/AUC 富集度验证** →(**v2.0**)交互式 HTML 总报告与 **3D 查看器**、一键打包导出、
headless CLI 与跨平台构建脚本。

> **v2.0 亮点(桌面级可发表虚拟筛选工作台)**:① **VS 规模吞吐**——支持数千到十万级
> 配体(进程池多核、分批+流式写结果、按配体断点续跑、吞吐统计);② **ROC/AUC 富集度
> 验证工作流**(用户提供已知活性 + 诱饵时输出 ROC 曲线/AUC/EF1%/5%);③ **交互式
> HTML 总报告** + 每个 Top 复合物的 **3D 查看器**;④ **一键数据打包导出 .zip**;
> ⑤ **headless CLI**(`python -m dockstudio run <project.json>`);⑥ **CITATION.cff**
> 与论文引用说明。同时并入三项前置增强:**对称性感知 RMSD**(自对接主口径,保留贪心
> 对比列)、**SMILES/CSV 配体批量输入**、**多核任务池**。
>
> **v1.1 亮点**:ttkbootstrap 现代界面(左导航 + 顶部操作栏 + 卡片化步骤);
> 输出目录用系统文件夹选择器;受体输入新增 **mmCIF(.cif/.mmcif)**;新增可选
> **分子对接准确性评估报告**;**帮助中心**;多项科学严谨性修正(PDBQT 元素/AutoDock
> 原子类型映射、pH 简单规则显式化、盒内离子/辅因子移除的如实披露、Top-K 文件名与
> 实际 K 一致等)。

## 特性与规格书对照

| 规格书章节 | 实现位置 | 状态 |
| --- | --- | --- |
| §3 环境与工具版本表 | `dockstudio/core/envinfo.py`,报告输出 `environment.json` | 已实现 |
| §4 输入盘点 CSV | `core/inventory.py` → `receptor_inventory.csv`, `ligand_list.csv` | 已实现 |
| §5 结构准备 | `core/structure.py`(**支持 PDB/mmCIF** 读取与清洗),`core/receptor.py`(PDBQT),`core/ligand.py`, `core/cofactor.py`(CCD 共晶参考) | 已实现 |
| §5.4 盒子决策树 | `core/box.py`,GUI ②结合位点;无共晶配体时强制人工选择/盲对接(探索性) | 已实现 |
| §6-7 Vina 对接 + 断点续跑 | `core/docking.py`, `core/batch.py`(`done.flag`、时间片、error.json) | 已实现 |
| §8 排名 → 精修 → Top-K | `core/ranking.py`(初筛矩阵、Top-12 精修、`final_top{K}.csv`/`_final_top{K}.json`,K=cfg.top_k) | 已实现 |
| §9 共晶配体回贴 | `core/selfdock.py`(重原子 RMSD,≤2.0 Å 为 PASS) | 已实现 |
| §10 PLIP | `core/interactions.py`(`plip_summary.json`, `plip_interactions_all.csv`) | 已实现 |
| §11 PyMOL 图件 | `core/visualize.py`(3D/composite/2D/TopK 拼接,ASCII `.pml`, `.pse`) | 已实现 |
| §12 报告与 QC | `core/reports.py`, `core/qc.py`;`reports/01_methods…`, `02_quality…`, `directory_tree.txt` | 已实现 |
| 准确性评估(新增) | `core/accuracy.py`;`reports/03_accuracy_assessment.md`, `accuracy_assessment.csv`, `accuracy_poses.csv` | v1.1 |
| §13 科学决策询问 | GUI 弹窗 + 引擎在“无盒子的受体”时明确报错而非猜测 | 已实现 |
| 对称性感知 RMSD(v2.0) | `core/symrmsd.py` + `selfdock.py`/`accuracy.py`(自对接与准确性报告主口径;贪心值保留为对比列) | v2.0 |
| SMILES/CSV 配体输入(v2.0) | `core/ligand.py`(ETKDG 3D 生成/去重/失败逐分子报告)、`inventory.py`、GUI ① | v2.0 |
| VS 规模吞吐 + 进程池(v2.0) | `core/batch.py`(ProcessPool、分批+流式写结果、按配体断点、吞吐统计)+ `pipeline.py` `run_summary.json` | v2.0 |
| ROC/AUC 富集度验证(v2.0) | `core/enrich.py`(AUC/EF1%/5%、ROC 曲线、SVG/CSV/JSON;仅当用户提供活性+诱饵) | v2.0 |
| 交互式 HTML 报告/3D 查看器(v2.0) | `core/vizhtml.py`(`reports/html/index.html` + `results/html_viewers/*.html`) | v2.0 |
| 数据打包 / .dsproj 工程(v2.0) | `core/project.py`(`.zip` 一键导出,含全部源文件+环境) | v2.0 |
| headless CLI(v2.0) | `__main__.py`:`python -m dockstudio run <project.json>`(逃生门,GUI 定位不变) | v2.0 |
| 引用 / 版本归档(v2.0) | `CITATION.cff` + README“如何引用”;DOI 注册由作者在 Zenodo 完成 | v2.0 |
| GUI | `gui/app.py` — ttkbootstrap 现代界面(左导航/顶部操作栏/卡片化);输出目录系统文件夹选择器;帮助中心 | v1.1/v2.0 |

## 快速开始(GUI)

**方式一(终端用户):** 在 Windows 构建机用 `packaging/online_installer/build_online_installer.bat`
重建后得到
[`dist/DockStudio_Setup_2.0.0.exe`](dist/DockStudio_Setup_2.0.0.exe)(在线安装器,
内含 Eric Studio 版权、图标与官方 Vina Windows 二进制)。安装到用户目录后,**首次启动**
会自动装配 PyMOL/RDKit/Meeko/PLIP 等依赖(需联网),之后即可正常使用。也可在 Windows
上运行 `packaging/online_installer/build_online_installer.bat` 自行重建安装包。

**方式二(开发/构建机):**

```bash
# 安装依赖(详见 packaging/打包说明.md)
pip install -r packaging/requirements-pip.txt     # Linux/macOS
# Windows 推荐 conda: conda env create -f packaging/environment.yml

# 启动
python -m dockstudio
```

界面 5 步(现代左侧导航界面):

1. **①输入**:添加受体 `*.pdb/*.ent/*.cif/*.mmcif`、配体库 `*.sdf/*.mol/*.smi/*.csv`;用“浏览文件夹…”
   选择输出目录。点“示例: 4DFR 演示数据”可载入内置演示。
2. **②结合位点**:选择受体 → 自动(共晶配体质心)/手动中心+尺寸/残基中心/盲对接,
   点“预览盒子”。
3. **③参数**:exhaustiveness、Top-K、CPU、并行任务数、pH、回贴/**准确性评估报告**/PLIP/
   可视化开关;v2.0 可启用“富集度验证模式”并提供已知活性 + 诱饵文件。
4. **④运行**:点“开始 / 断点续跑”,看实时日志;中途关闭软件后重开并再次点击即可续跑。
5. **⑤结果**:查看 Top-K 表、打开报告、HTML 总报告/3D 查看器、打包导出 .zip。

需要帮助时使用顶部 **帮助 → 帮助中心** 在线查看文档。

## 论文引用 DockStudio(How to cite)

发表论文或预印本时可按如下格式引用(各期刊风格可微调):

> Eric Studio. DockStudio: 一站式自动化批量分子对接平台(版本 2.0.0).
> Zenodo, 2026. https://doi.org/10.5281/zenodo.22655142

仓库根含机器可读的 **`CITATION.cff`**(Zenodo 会据此生成引文)。
v2.0.0 已注册 DOI **10.5281/zenodo.22655142**(记录族 concept DOI
`10.5281/zenodo.22655141`)。后续每个新版本:推 GitHub Release → Zenodo 自动铸新
版本 DOI → 把新 DOI 回填到 `CITATION.cff` 的 `doi:` 字段并重新提交。
(软件本身不自动执行 DOI 注册;DOI 由作者在 Zenodo 开启 GitHub 集成后自动/手动完成。)

## 无头冒烟测试(不需要 GUI)

```bash
python tests/run_smoke.py 150 --out examples/demo_out     # 一次跑 ~150 秒
python tests/run_smoke.py 150 --out examples/demo_out     # 再跑几次直到 interrupted=false
```

`examples/demo_out` 已包含一份完整运行示例(89 个文件),可直接当作输出格式范本。

## headless CLI(v2.0,集群/自动化逃生门)

软件保持 **GUI 优先**,但提供命令行运行同一工程文件的能力(方便服务器/脚本复用):

```bash
python -m dockstudio --version
python -m dockstudio run examples/demo_out/config.json --out /tmp/so
python -m dockstudio run my.dsproj --phases prepare,screen,refine,selfdock,analyze,enrich,visualize,report,html
```

工程文件可以是 GUI 生成的 `config.json` 或 v2.0 的 `.dsproj`;`--out` 可覆盖输出目录,
`--phases` 可指定阶段子集,`--time` 设置时间片(秒)。退出码:0=正常结束;1=流水线异常;
2=参数/工程文件错误;3=时间片中断(可续跑)。

## 打包分发(把 PyMOL/Vina 全部带上)

见 **`packaging/打包说明.md`**。推荐 conda-pack 生成绿色免安装目录:
Windows 执行 `packaging/condapack/build_win_condapack.bat`,产物拷给终端用户双击
`DockStudio.bat` 即用。

## 目录结构

```
dockstudio/
├── gui/app.py           # ttkbootstrap 现代中文图形界面(回退标准 ttk)
├── core/                # 引擎(无 GUI 依赖,可编程/可测试)
│   ├── models.py        # RunConfig / BoxDef 数据模型
│   ├── structure.py     # PDB/mmCIF 解析、盘点与清洗(行级 PDB 输出与 gemmi 对齐)
│   ├── inventory.py     # §4 输入盘点 CSV
│   ├── receptor.py      # §5.1 受体准备(meeko,坏残基远处删除/近处询问)
│   ├── ligand.py        # §5.2 配体准备(逐分子、pH 简单规则、PDBQT)
│   ├── cofactor.py      # §5.3 CCD 共晶配体 SDF(坐标取本地 PDB/mmCIF)
│   ├── box.py           # §5.4 盒子决策树
│   ├── docking.py       # §6 Vina 封装(result.json/poses/best/done.flag)
│   ├── batch.py         # §7 批量执行/续跑/时间片
│   ├── ranking.py       # §8 矩阵/精修/Top-K(final_top{K})
│   ├── selfdock.py      # §9 回贴 RMSD(对称性感知主口径 + 贪心对比列)
│   ├── symrmsd.py       # v2.0 对称等价类感知 RMSD(RDKit 自同构)
│   ├── accuracy.py      # 对接准确性评估报告(v1.1)
│   ├── enrich.py        # v2.0 ROC/AUC/EF 富集度验证
│   ├── vizhtml.py       # v2.0 交互式 HTML 总报告 + 3D 查看器
│   ├── project.py       # v2.0 .dsproj 工程 + 一键 .zip 数据打包
│   ├── interactions.py  # §10 PLIP
│   ├── visualize.py     # §11 PyMOL
│   ├── reports.py / qc.py # §12 报告与自动质控
│   └── pipeline.py      # 编排各阶段(引擎核心)
├── examples/            # 演示数据 + demo_out 示例输出
├── tests/               # 冒烟驱动与单测(mmCIF/准确性/PDBQT)
└── packaging/           # conda-pack / PyInstaller / NSIS 打包
```

## 开发与测试

```bash
python -m compileall -q dockstudio          # 语法检查
python tests/run_smoke.py 150 --out /tmp/so  # 端到端
# 单测(轻量,不依赖 vina/PyMOL):
python -m pytest tests -q -k "not smoke"
```

## 重要声明(诚实性)

本工具遵循“**绝不捏造数据**”的原则:

- 所有数值均来自真实运行写入的文件;报告只列出真实执行的步骤。
- 共晶配体回贴 RMSD、PLIP 交互、结合能等若失败或异常,会 **如实报告并给出建议**,
  不会静默省略。
- 无共晶配体的受体全部结果标注为 **探索性**;其回贴验证如实标注“不可评估”。
- **富集度验证(ROC/AUC/EF)只在用户提供已知活性 + 诱饵数据时启用**,所有数值来自真实
  对接得分;缺失数据时如实跳过,绝不自动生成活性/诱饵或伪造曲线。
- **自对接 RMSD 主口径为对称性感知**(等价类内原子互换不计误差);对称性不可用时如实
  退回贪心最近邻并在 `rmsd_method`/`symmetry_used` 列标注;贪心值始终保留作对比。
- 受体清洗会移除水/离子/配体;位于搜索盒内或其附近的被移除离子/辅因子会在报告
  **如实披露清洗提示**。
- 质子化采用文档化简单规则(目标 pH 可调;羧酸 pH≥5 去质子、碱性胺 pH≤8.5 质子化),
  不伪装成严格 pKa 处理。

## 许可

项目代码采用 MIT License(见仓库根目录 LICENSE)。依赖软件各自许可:
AutoDock Vina、Meeko、RDKit、gemmi、ProDy、PLIP、PyMOL(open-source)。
分发前请确认符合各依赖许可证要求(学术使用一般无碍;商业分发请自查)。
