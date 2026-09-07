# DockStudio 跨平台构建(Linux / macOS)— v2.0

> © 2026 Eric Studio · 本目录文件面向“把 DockStudio 从 Windows-only 扩展到
> Linux/macOS”。Windows NSIS 安装包仍由 `packaging/online_installer/` 构建。

## 产物形态

| 平台 | 产物 | 说明 |
| --- | --- | --- |
| Windows | `DockStudio_Setup_2.0.0.exe`(NSIS 在线安装器) | 见 `packaging/online_installer/` |
| Linux | `DockStudioPortable.tar.gz` + `DockStudio.sh` | conda-pack 绿色目录,`bin/python` 自带全引擎工具链 |
| macOS | 同上(需 macOS 构建机实产;脚本相同) | conda-pack 绿色目录 |

引擎是跨平台的(Python + RDKit/Meeko/PLIP + Vina CLI/python-vina),
因此 GUI 与 headless CLI 在三个平台行为一致。打包差异仅来自环境与安装器。

## 构建

```bash
# 1) 建环境(任一平台,需先装 Miniconda/conda-forge)
conda env create -f packaging/environment.yml     # 环境名 dockstudio
conda activate dockstudio

# 2) 产出 conda-pack 绿色目录 + 启动脚本 + MD5
bash packaging/cross/build_condapack_unix.sh dockstudio dist
```

产出:

- `dist/DockStudioPortable.tar.gz`(可拷走/分发)
- `dist/DockStudioPortable/`(绿色目录;启动:`./DockStudio.sh` 或解压后同结构)
- `dist/DockStudioPortable.tar.gz.md5`

### macOS 实产说明(诚实标注)

本仓库开发沙箱为 Linux,**无 macOS 构建机**;上面脚本在 macOS(Apple Silicon 或 Intel)
上同样适用,但需在 macOS 上执行一次以产生并验证 `.tar.gz` 实产。若你的分发对象是
macOS 用户,请在 mac 构建机上跑一次该脚本并做一次冒烟;结果以 mac 本机实测为准。

## 依赖三文件需同步

新增任何 Python 依赖时,请同步更新(否则打包会缺库):

- `packaging/requirements-pip.txt`
- `packaging/environment.yml`
- `packaging/online_installer/setup_env.bat`

## headless CLI(服务器/集群)

跨平台产物里同样可用无头模式:

```bash
./bin/python -m dockstudio --version
./bin/python -m dockstudio run <project.dsproj> --out /tmp/out
```

详见仓库根 `README.md` 的 “headless CLI” 一节。
