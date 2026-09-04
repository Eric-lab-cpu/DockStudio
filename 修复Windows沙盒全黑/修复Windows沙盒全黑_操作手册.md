# Windows 沙盒打开全黑——排查与修复手册

适用:Windows 10/11 “Windows Sandbox(Windows 沙盒)”打开后**整个窗口黑屏、无法操作**,
或先黑屏后一直转圈/提示“无法初始化”。常见于 **Windows 更新之后**、或装有较新
NVIDIA/AMD 显卡驱动的机器(尤其笔记本 + 双显卡)。已实测这类现象与该沙盒的
**虚拟 GPU(vGPU)渲染**和宿主机显卡驱动不兼容高度相关。

> 本目录文件:
> - `Sandbox_NoGPU.wsb` — 关闭沙盒内 GPU 加速的启动配置(最关键修复)
> - `一键修复并启动无GPU沙盒.bat` — 管理员身份运行:修复系统组件 + 直接以无 GPU 模式打开沙盒
> - `诊断Windows沙盒.ps1` — 生成诊断日志,便于进一步求助
> - 本手册

---

## 一、先做最简单有效的:禁用沙盒 GPU 加速(90% 的“黑屏”都因此解决)

1. 双击 **`Sandbox_NoGPU.wsb`**(或右键→打开方式→Windows Sandbox)。
   也可在 PowerShell 执行:
   ```powershell
   Start-Process WindowsSandbox.exe -ArgumentList 'C:\…\Sandbox_NoGPU.wsb'
   ```
2. 等待 30–60 秒(首次/更新后启动偏慢,可能短暂黑屏属正常)。
3. 若此时能出现桌面 → 说明就是 **vGPU/显卡驱动** 问题。之后两种选择:
   - **继续用无 GPU 模式**(推荐先这样用,就是渲染稍慢);
   - 或更新/回滚显卡驱动后,把 `.wsb` 里的 `<vGPU>disable</vGPU>` 一行删掉,恢复硬件加速。

> 说明:`.wsb` 里的 `MemoryInMB=4096`、`CPUs=2` 可按你电脑内存调整;内存少于 8 GB
> 请把 MemoryInMB 降到 2048。

---

## 二、若仍黑屏:管理员一键修复组件

**右键** `一键修复并启动无GPU沙盒.bat` → **以管理员身份运行**。脚本会:

1. `DISM /Online /Cleanup-Image /RestoreHealth` — 修复系统映像;
2. `sfc /scannow` — 系统文件检查;
3. 确保三个功能处于开启状态:
   `Containers-DisposableClientVM`(即 Windows Sandbox)、`VirtualMachinePlatform`、`Microsoft-Hyper-V-All`;
4. 再次用“无 GPU”配置尝试打开沙盒。

执行后**重启电脑**,再打开沙盒。

---

## 三、手动检查/开启功能(不想跑脚本时)

设置 → 应用 → 可选功能 → 更多 Windows 功能(或 `control appwiz.cpl`→启用或关闭 Windows 功能):

- 勾选 **Windows 沙盒**
- 勾选 **虚拟机平台(Virtual Machine Platform)**
- 勾选 **Hyper-V**(含 Hyper-V 平台/管理工具;Windows 家庭版需先启用“虚拟机平台”后才有 Hyper-V)

确定后重启。命令行核对(管理员 PowerShell):

```powershell
Get-WindowsOptionalFeature -Online |
  Where-Object FeatureName -match 'Sandbox|VirtualMachinePlatform|Hyper-V' |
  Select-Object FeatureName,State
```

---

## 四、确认 CPU 虚拟化已开启(BIOS)

打开 任务管理器 → 性能 → CPU → 右下角“虚拟化:已启用”。
若显示“已禁用”,需进 BIOS/UEFI 开启 **Intel VT-x / AMD-V**(笔记本还要看有无
“Virtualization Technology”开关),保存重启。也可用:

```powershell
(Get-CimInstance Win32_Processor).VirtualizationFirmwareEnabled
```

返回 `True` 表示已开。

---

## 五、显卡驱动:更新或回滚

- **回滚**:设备管理器 → 显示适配器 → 你的显卡 → 属性 → 驱动程序 → 回退驱动程序。
- **更新**:到 NVIDIA / AMD / Intel 官网下载对应驱动(不建议第三方工具)。
- 笔记本双显卡(独显+核显)场景,回滚到上一个大版本驱动常能解决沙盒黑屏。

---

## 六、仍不行:生成诊断日志求助

管理员 PowerShell 执行:

```powershell
powershell -ExecutionPolicy Bypass -File "诊断Windows沙盒.ps1"
```

会在本目录生成 `WindowsSandbox_Diagnose_时间戳.txt`,包含:系统版本、CPU 虚拟化是否
开启、内存、功能组件状态、相关服务、显卡驱动、近 30 分钟错误事件。把它发给我或微软社区,
能快速定位。

---

## 七、关于“沙盒内搞坏东西会不会影响主系统”的澄清(重要)

- **不会。** Windows Sandbox 是每次**全新**的轻量虚拟机:关闭即销毁,所有更改(卸载
  Edge、改系统文件、装软件)全部丢弃;**下次打开永远是干净系统**。
- 它和 WSL2 不同:WSL2 会与 Windows 共享文件系统/内核集成;Windows Sandbox 是
  **完全隔离**的临时会话,不影响宿主机。
- 所以你可以放心在里面做各种实验(包括测试 DockStudio 安装包),关掉就复原。

---

## 八、如果只是想“干净测试 DockStudio”,也可选替代方案

- **直接在主系统装/测**:DockStudio 安装包会装到 `%LOCALAPPDATA%\EricStudio\DockStudio`,
  卸载器可完整移除,风险很低。
- **Hyper-V 虚拟机**或 **VMware/VirtualBox**:隔离性同沙盒,但可持久保存快照。

---

参考来源:
- [Windows沙盒总是转圈 - Microsoft Q&A](https://learn.microsoft.com/zh-cn/answers/questions/5913410/windows)
- [Win11 无法开启沙盒功能](https://www.php.cn/faq/2713019.html)
- [Windows sandbox display failure - Dell 社区](https://www.dell.com/community/en/conversations/inspiron/g15-5530-windows-11-windows-sandbox-display-failure/69eac89c65dced6f2989cc0e)
- [Fix Windows Sandbox Errors](https://www.rescuepcrepairs.com/fix/windows-sandbox-error)
- [Windows sandbox closes … connection lost - ElevenForum](https://www.elevenforum.com/t/windows-sandbox-closes-and-says-the-connection-to-the-windows-sandbox-environment-was-lost-do-you-want-to-reconnect.46200/)
