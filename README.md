# DockStudio — 一站式自动化批量分子对接平台

> **© 2026 Eric Studio. All rights reserved. / 版权所有 © 2026 Eric Studio。保留所有权利。**
> DockStudio 由 Eric Studio 提供。

**DockStudio** 是一个 **图形界面为主、引擎可编程** 的分子对接工具,内置一套严谨、
可复现、可发表级的自动化流程:受体/配体盘点 → 结构准备(Meeko)→ 搜索盒子决策 →
AutoDock Vina 批处理(断点续跑)→ 排名精修 Top-K → 共晶配体回贴验证 → PLIP 相互
作用分析 → PyMOL 论文级图件与可编辑 `.pml`/`.pse` → 方法与质量报告。

## 特性与规格书对照

| 规格书章节 | 实现位置 | 状态 |
| --- | --- | --- |
| §3 环境与工具版本表 | `dockstudio/core/envinfo.py`,报告输出 `environment.json` | 已实现 |
| §4 输入盘点 CSV | `core/inventory.py` → `receptor_inventory.csv`, `ligand_list.csv` | 已实现 |
| §5 结构准备 | `core/structure.py`(清洗/altloc/修饰残基),`core/receptor.py`(PDBQT),`core/ligand.py`, `core/cofactor.py`(CCD 共晶参考) | 已实现 |
| §5.4 盒子决策树 | `core/box.py`,GUI ②结合位点;无共晶配体时强制人工选择/盲对接(探索性) | 已实现 |
| §6-7 Vina 对接 + 断点续跑 | `core/docking.py`, `core/batch.py`(`done.flag`、时间片、error.json) | 已实现 |
| §8 排名 → 精修 → Top-K | `core/ranking.py`(初筛矩阵、Top-12 精修、`final_top5.csv`/`_final_top5.json`) | 已实现 |
| §9 共晶配体回贴 | `core/selfdock.py`(重原子 RMSD,<2.0 Å 为 PASS) | 已实现 |
| §10 PLIP | `core/interactions.py`(`plip_summary.json`, `plip_interactions_all.csv`) | 已实现 |
| §11 PyMOL 图件 | `core/visualize.py`(3D/composite/2D/TopK 拼接,ASCII `.pml`, `.pse`) | 已实现 |
| §12 报告与 QC | `core/reports.py`, `core/qc.py`;`reports/01_methods…`, `02_quality…`, `directory_tree.txt` | 已实现 |
| §13 科学决策询问 | GUI 弹窗 + 引擎在“无盒子的受体”时明确报错而非猜测 | 已实现 |

## 快速开始(GUI)

**方式一(终端用户):** 运行仓库内已产出的
[`dist/DockStudio_Setup_1.0.0.exe`](dist/DockStudio_Setup_1.0.0.exe)(在线安装器,
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

界面 5 步:

1. **①输入**:添加受体 `*.pdb`、配体库 `*.sdf`,选输出目录。点“示例:4DFR 演示数据”
   可载入内置演示。
2. **②结合位点**:选择受体 → 自动(共晶配体质心)/手动中心+尺寸/残基中心/盲对接,
   点“预览盒子”。
3. **③参数**:exhaustiveness、Top-K、CPU、pH、回贴/PLIP/可视化开关。
4. **④运行**:点“开始 / 断点续跑”,看实时日志;中途关闭软件后重开并再次点击即可续跑。
5. **⑤结果**:查看 Top-K 表、打开报告与图件目录。

## 无头冒烟测试(不需要 GUI)

```bash
python tests/run_smoke.py 150 --out examples/demo_out     # 一次跑 ~150 秒
python tests/run_smoke.py 150 --out examples/demo_out     # 再跑几次直到 interrupted=false
```

`examples/demo_out` 已包含一份完整运行示例(89 个文件),可直接当作输出格式范本。

## 打包分发(把 PyMOL/Vina 全部带上)

见 **`packaging/打包说明.md`**。推荐 conda-pack 生成绿色免安装目录:
Windows 执行 `packaging/condapack/build_win_condapack.bat`,产物拷给终端用户双击
`DockStudio.bat` 即用。

## 目录结构

```
dockstudio/
├── gui/app.py           # 中文 tkinter 图形界面入口
├── core/                # 引擎(无 GUI 依赖,可编程/可测试)
│   ├── models.py        # RunConfig / BoxDef 数据模型
│   ├── structure.py     # PDB 解析/盘点/清洗(行级格式化,已与 gemmi 输出逐字节核对)
│   ├── inventory.py     # §4 输入盘点 CSV
│   ├── receptor.py      # §5.1 受体准备(meeko,坏残基远处删除/近处询问)
│   ├── ligand.py        # §5.2 配体准备(逐分子、pH 规则、PDBQT)
│   ├── cofactor.py      # §5.3 CCD 共晶配体 SDF(坐标取本地 PDB)
│   ├── box.py           # §5.4 盒子决策树
│   ├── docking.py       # §6 Vina 封装(result.json/poses/best/done.flag)
│   ├── batch.py         # §7 批量执行/续跑/时间片
│   ├── ranking.py       # §8 矩阵/精修/Top-K
│   ├── selfdock.py      # §9 回贴 RMSD
│   ├── interactions.py  # §10 PLIP
│   ├── visualize.py     # §11 PyMOL
│   ├── reports.py / qc.py # §12 报告与自动质控
│   └── pipeline.py      # 编排各阶段(引擎核心)
├── examples/            # 演示数据 + demo_out 示例输出
├── tests/               # 冒烟驱动
└── packaging/           # conda-pack / PyInstaller 打包
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
- 无共晶配体的受体全部结果标注为 **探索性**。
- 质子化采用文档化简单规则,不伪装成严格 pKa 处理。

## 许可

项目代码采用 MIT License(见仓库根目录 LICENSE)。依赖软件各自许可:
AutoDock Vina、Meeko、RDKit、gemmi、ProDy、PLIP、PyMOL(open-source)。
分发前请确认符合各依赖许可证要求(学术使用一般无碍;商业分发请自查)。
