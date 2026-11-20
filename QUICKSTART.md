# NonlinearFEM v0.2.0 Windows 便携版

第二版便携预发布包，面向 Windows 10/11 x64，支持 Q4 与 H8。日常建模、求解和查看云图无需另装 Python、CMake 或编译器。

1. 完整解压 ZIP 到可写文件夹，保留 `_internal/`、`build/` 和其他目录。
2. 双击 `NonlinearFEM-Workbench.exe`。Q4 在主窗口，H8 点击右上角“三维 H8 工作台”。
3. 设置尺寸、网格和材料；可从示例开始，或打开自定义边界窗口选择节点/边/面，逐方向设置约束和载荷。
4. 选择静力或动力，点击运行。结果放在软件目录 `results/` 下。
5. 选择位移、第二 PK 应力、GL 应变或 von Mises（第二 PK）云图。蓝到红为数值由小到大；可固定或手动设定色标、查看极值、导出 CSV。默认变形倍率为 1。

完整说明见 `manuals/USER.md`，原理见 `manuals/THEORY.md`，标准算例见 `manuals/EXAMPLES.md`。也附有 `CONTOURS.md`、`CUSTOM_BOUNDARIES.md` 和 `H8_WORKBENCH.md`。

## 使用范围

- GUI 仍生成矩形 Q4 和长方体 H8 网格，最多 600 节点；没有任意几何造型、接触、塑性或损伤。
- 单位必须一致。动力约束仅支持固定零位移；自动时间步仅依据初始状态。
- 应力为第二 PK，应变为 Green–Lagrange，剪应变采用工程分量。云图为单元平均，尚未平滑。
- 自定义总力在所选节点间均分，不等于按面积积分的均布牵引。

## 运行标准算例

包内附标准输入、六个现有验证程序和内置执行入口。PowerShell 在软件目录执行（输出目录必须尚不存在）：

```powershell
Start-Process .\NonlinearFEM-Workbench.exe -ArgumentList '--standard-cases --output results/standard_check' -Wait
```

完成后查看 `results/standard_check/REPORT.md` 与 `report.json`。这条入口使用内置 Python，不要求另装 Python。源代码开发环境仍可执行 `python tools/run_standard_cases.py --build-dir build`。

## 包内目录

`_internal` 是图形运行环境；`build` 含求解器、动力示例及标准算例所需验证程序；`manuals` 含三本手册与实测附件；`examples` 含模型。`src`、`include`、`gui`、`tests`、`tools` 和 CMake 配置提供对应源码，便于复现和学习。源码编译仍需要开发工具。

遇到问题先看 `RELEASE_NOTES.md` 和用户手册。包未签名，发布检查不等于所有 Windows 机器上的兼容性认证。AI 辅助开发声明与许可证保留在 README、LICENSE 和第三方许可中。
