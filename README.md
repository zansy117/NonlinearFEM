# NonlinearFEM

使用 C++17 实现的二维 Q4 几何非线性有限元研究项目。提供滑块建模界面、文本模型输入、非线性求解和位移结果查看，当前重点是小规模算例和数值验证。

## 当前功能

- 四节点四边形（Q4）单元，基于参考构形计算。
- 变形梯度、Green–Lagrange 应变和第二 Piola–Kirchhoff 应力。
- 当前材料采用平面应力弹性矩阵，将 Green–Lagrange 应变映射到第二 PK 应力；非线性主要来自几何关系。
- 单元内力、切线刚度及多单元全局组装。
- 节点力和指定位移边界条件。
- Newton–Raphson 求解、分步加载及自适应缩步（cutback）。

## 编译

需要支持 C++17 的编译器、CMake 3.16 或更高版本，以及与所选生成器对应的构建工具。Git 仅用于获取和管理代码。

在项目根目录执行：

```powershell
cmake -S . -B build
cmake --build build
```

如果使用 Visual Studio 等多配置生成器，可以执行：

```powershell
cmake --build build --config Debug
```

## 运行测试

默认启用 CTest 测试注册。编译完成后执行：

```powershell
ctest --test-dir build --output-on-failure
```

多配置生成器需指定与构建一致的配置：

```powershell
ctest --test-dir build -C Debug --output-on-failure
```

只运行 2×2 网格求解测试：

```powershell
ctest --test-dir build -R "^test_solver_2x2$" --output-on-failure
```

查看已注册测试：

```powershell
ctest --test-dir build -N
```

当前环境已注册 31 项测试（包含检测到 Python 3 后注册的界面模型与动力输出测试）。以下列出基础测试，动力学测试见后文：

| 范围 | 测试程序 |
| --- | --- |
| 模型输入与导出 | `test_model_io` |
| 运动学、应变、材料 | `test_kinematics`、`test_Gstrain`、`test_material` |
| Q4、组装、边界条件 | `test_q4`、`test_assembly`、`test_boundary` |
| Newton 求解 | `test_solver`、`test_solver_2x2` |
| 悬臂板网格加密验证 | `test_cantilever_convergence` |
| 大变形拉伸解析验证 | `test_uniaxial_tension` |
| 切线刚度验证 | `test_tangent_stiffness`、`test_tangent_stiffness_distorted` |
| 载荷步进 | `test_load_stepping`、`test_load_stepping_cutback`、`test_load_stepping_displacement` |
| Patch test | `test_q4_patch_multi`、`test_q4_patch_stress`、`test_q4_patch_distorted` |
| 全局组装验证 | `test_global_assembly_verification` |

测试涉及解析仿射变形、规则和畸变网格、切线有限差分、全局组装对照、自由度平衡，以及位移约束。2×2 求解测试还将支座反力与外载的整体平衡纳入通过条件。测试通过只表示已有检查满足要求。

## 目录

```text
include/nonlinearfem/  接口和数据结构
src/                  实现代码
tests/                测试与算例
CMakeLists.txt        构建和 CTest 配置
```

## 添加测试

在已有测试文件中增加检查后，重新编译并运行 CTest 即可。

新建测试程序时，在 `CMakeLists.txt` 中用 `add_executable(...)` 列出测试文件及依赖源文件，再将目标名称加入 `NONLINEARFEM_TARGETS`。名称以 `test_` 开头的目标会自动注册到 CTest。修改配置后重新运行 `cmake -S . -B build`。

测试通过时返回 `0`，失败时返回非零值；仅打印失败信息不会使 CTest 判定失败。

## 当前限制

- 主程序支持 `.nfem` 文件；矩形板可以使用 `LaunchGUI.cmd` 滑块界面建模。
- GUI 当前支持五种矩形板边界模板；任意 Q4 网格需要编辑 `.nfem` 文件。尚无任意几何网格生成或应力恢复/应力云图。
- 当前使用稠密全局刚度矩阵和高斯消元，主要面向小规模验证算例。
- 当前材料实现不包含塑性或损伤等历史相关行为。
- 求解器失败路径的专门回归验证、更多独立参考算例及网格收敛研究仍待完善。
## 大变形均匀拉伸验证

`test_uniaxial_tension` 使用单位正方形的 2×2 Q4 网格，厚度为 0.2，杨氏模量为 2.1e5，泊松比为 0.3。左边界限制水平位移，左下角限制竖向位移；右边界施加按参考边长分配的水平节点力，上下边界自由。

对于当前材料模型，令轴向伸长比为 λ，横向无应力条件给出：

- 轴向 Green–Lagrange 应变：Exx = (λ² − 1)/2。
- 横向应变：Eyy = −ν Exx。
- 横向伸长比：λy = sqrt(1 − ν(λ² − 1))。
- 第二 PK 轴向应力：Sxx = Young × Exx。
- 名义轴向应力：Pxx = λ Sxx。
- 总载荷：Force = Pxx × 参考高度 × 参考厚度。
- 解析位移：ux = (λ − 1)X，uy = (λy − 1)Y。

测试由解析关系生成 10 个目标载荷（λ 从 1.05 到 1.50），逐点进行力控制 Newton 求解，并以前一点的数值解作为初值。解析位移只用于比较，不作为求解器的位移边界条件。每点检查所有节点的位移、自由度平衡、整体反力平衡和端部位移单调性，并输出载荷—位移数据。

```powershell
ctest --test-dir build -R "^test_uniaxial_tension$" -V
```

这是当前材料模型的解析一致性验证，不是实际材料在 50% 伸长下的物理准确性证明，也不是网格收敛研究。
## 悬臂板网格加密验证

`test_cantilever_convergence` 使用长 4、高 1、厚 0.2 的平面应力矩形板，杨氏模量 1e5，泊松比 0.3。左端两个位移分量固定；右端沿参考边界施加均匀、方向固定的向下分布力，端点节点分配半权重。所有数值需采用一致单位。

测试采用 4×1、8×2、16×4、24×6 四组网格，总向下力取 0.1、12.5、25、37.5、50。每个目标载荷从零状态开始进行分步加载，默认以目标载荷的 25% 为初始步长，允许自适应缩步。输出参考右边界中点的水平及竖向位移；若该点不是网格节点，使用相邻边界节点插值。

每个算例检查完成全部加载、有限位移、约束误差、自由自由度残差和 X/Y 整体反力平衡。每组网格检查位移随载荷单调增加；每个载荷点检查最后一次网格加密的相对位移变化小于前一次且低于 5%。不同加密比下的变化量用于趋势检查，不用于宣称收敛阶。

小载荷采用 Timoshenko 梁端部位移近似作为参考：

```text
I = t H^3 / 12
G = E / (2 (1 + nu))
delta = F L^3 / (3 E I) + F L / ((5/6) G t H)
```

该比较允许 10% 的差异，因为二维连续体的端部约束和载荷分布与一维梁近似并不完全相同。大变形阶段不使用线性梁理论作为准确解。

首次验证的最大载荷结果（向下位移绝对值）：

| 网格 | 端部中点位移 |
| --- | ---: |
| 4×1 | 0.44567924 |
| 8×2 | 0.57922941 |
| 16×4 | 0.63102510 |
| 24×6 | 0.64214823 |

最后两组网格在最大载荷下相差约 1.73%；小载荷下最细网格位移为 0.00131429，梁近似为 0.00134240，相差约 2.09%。这些结果支持网格稳定趋势，不构成大变形情况下的独立准确性证明。仍需更细网格或独立求解器对照来量化剩余误差。

```powershell
ctest --test-dir build -R "^test_cantilever_convergence$" -V
```
## 滑块建模界面（Windows）

先按上面的步骤编译 C++ 项目。界面还需要 Python 3（包含 tkinter；不需要额外的 pip 包）。在项目根目录双击 `LaunchGUI.cmd` 即可启动，也可以运行：

```powershell
python gui/app.py
```

1. 选择所需的矩形板边界模板。切换模板会恢复对应的示例参数。
2. 拖动滑块设置长度、高度、两个方向的单元数、参考厚度、杨氏模量、泊松比和右端总力。数值框可以精确输入，回车后同步滑块；超出滑块范围的有效数值仍以数值框为准。杨氏模量滑块采用对数刻度。
3. 查看实时原始网格预览：蓝点代表约束节点，橙箭头代表载荷方向。悬臂板右边界为竖向力，正值向上；拉伸板为水平力，正值向右。载荷按参考边界长度分配，总力不是每个节点的力。
4. 点击“运行计算”。计算在独立进程中运行，不会阻塞界面；可以点击“停止”。
5. 成功后直接在窗口内查看变形和位移，切换位移分量或拖动变形放大滑块。原始网格使用灰色虚线，颜色代表各单元节点位移的平均值，不代表应力。
6. “导出模型”保存 `.nfem` 输入；“打开完整结果页面与下载数据”打开本次结果页。

每次计算自动建立 `results/时间戳/`，其中 `model.nfem` 是当次输入，`solution/` 包含 `results.html`、`nodes.csv` 和 `history.csv`。输入与结果一起保存，便于复现。结果目录已加入 Git 忽略规则。

GUI 当前限制为 600 节点，以控制稠密求解器的计算成本。极端参数仍可能不收敛；失败时不显示未完成结果。计算期间修改输入，已完成结果仍属于旧参数，界面会提醒重新计算。窗口不会绘制应力或加载过程动画。

## 文本模型输入与结果导出

也可以不启动界面，直接执行示例：

```powershell
.\build\NonlinearFEM.exe examples/tension.nfem build/tension_output
.\build\NonlinearFEM.exe examples/cantilever.nfem build/cantilever_output
```

多配置生成器使用对应的 `build/Debug/NonlinearFEM.exe` 或 `build/Release/NonlinearFEM.exe`。输出目录必须不存在，程序不会覆盖已有计算。

`.nfem` 是严格的逐行文本格式，支持 `#` 注释。文件以 `NFEM 1` 开头，以 `end` 结束；节点和单元不额外填写 ID，按记录出现顺序从 0 编号。所有单位必须一致。

| 段落 | 格式与说明 |
| --- | --- |
| 材料 | `material E nu`，当前平面应力弹性模型 |
| 厚度 | `thickness t`，正数 |
| 节点 | `nodes N`，后跟 N 行 `x y` |
| 单元 | `elements M`，后跟 M 行 `n0 n1 n2 n3`，节点按逆时针顺序排列 |
| 外载 | `loads K`，后跟 K 行 `node direction value`，重复节点力累加 |
| 位移约束 | `constraints C`，后跟 C 行 `node direction value`，同一自由度不可重复约束 |
| Newton（可选） | `newton maxIterations residualTolerance displacementTolerance` |
| 加载（可选） | `stepping initialIncrement minimumIncrement maximumIncrement cutbackFactor growthFactor maxRetries maxLoadSteps` |

`direction` 为 0（x）或 1（y）。载荷和非零位移约束均按加载系数从 0 比例增加到 1；初始位移为零。此版本支持单一材料、统一厚度、节点死载荷和 Q4 网格，不支持多材料、接触或随动力。

程序检查模型尺寸、节点编号、重复约束、有限数值和初始网格合法性。文本输入最多 2000 节点，但稠密求解器在该规模可能很慢，建议先使用小网格。

成功计算输出：

- `nodes.csv`：参考坐标、位移、位移大小、变形坐标和约束自由度反力。自由自由度反力列填 0，不作为残差输出。
- `history.csv`：已接受加载步的系数、增量、Newton 迭代次数、缩步重试次数和相对残差，不包含逐步位移。
- `results.html`：不依赖网络的最终位移结果页面，可切换分量、调整显示放大倍数并下载 CSV。

`test_model_io` 验证文件输入后解析拉伸解、错误输入拒绝、CSV 位移和反力、输出目录保护。`test_gui_model` 验证模板网格/约束/总力及生成文件经真实 C++ 求解器计算后的拉伸解析解和悬臂板参考位移。
## 显式动力时程与变形播放

动力核心包括 Q4 集中质量、独立内力组装、无阻尼中心差分积分、全局应变能与固定步长分析入口。动力与静力已合并到同一个建模入口；仍保留固定自由振动示例，供独立验证与演示。

编译后双击项目根目录的 `LaunchDynamics.cmd`，点击“计算并查看示例”。也可以执行：

```powershell
cmake --build build --target dynamic_demo
python gui/dynamics.py
```

窗口显示观察节点的 ux/uy 位移曲线、全模型动能/应变能/总能量曲线。拖动时间滑块或点击“播放”查看保存帧，变形放大倍数只改变显示。播放速度按保存帧设置，不等于物理时间。原始网格为灰色虚线，变形网格为蓝色；该窗口不显示应力。

内置示例：长 4、高 1、厚 0.2，4×2 Q4 网格，E=1000、nu=0.3、参考密度=1；左端固定，初始弯曲位移幅值 1e-4，初速度为零，无外力及阻尼。时间步 2e-4，推进 10000 步至 t=2，每 25 步保存一帧，含初始和最终状态共 401 帧。观察节点为右边界中点（节点 9）。单位需一致。

也可直接计算导出，并打开结果目录：

```powershell
.\build\dynamic_demo.exe results/my_dynamic_run
python gui/dynamics.py results/my_dynamic_run
```

输出目录必须不存在。多配置生成器使用对应的 `build/Debug/dynamic_demo.exe` 或 `build/Release/dynamic_demo.exe`。

导出内容：

- `nodes.csv` / `elements.csv`：参考节点与单元连接。
- `history.csv`：观察节点位移、速度、加速度，以及全模型动能、应变能和总能量。
- `frames.csv`：每个保存时刻全部节点的位移，供变形播放使用。

可复用的 `writeDynamicResults()` 接口声明在 `DynamicOutput.h`，输入网格、动力分析结果、观察节点编号和新输出目录。求解与导出分离，不改变时间积分算法。保存间隔越密、网格越大，内存和文件占用越大。

仅无外力、无阻尼且固定约束的自由振动可用机械能漂移衡量积分表现；一般受迫问题尚未计算外力做功，不能把机械能变化直接视为误差。初始刚度估计的时间步不保证大变形全过程稳定，稳定步长也不等于满足精度要求的步长。

动力相关自动测试：`test_lumped_mass`、`test_explicit_dynamics`、`test_q4_explicit`、`test_q4_energy`、`test_global_energy`、`test_multi_q4_explicit`、`test_stable_timestep`、`test_dynamic_analysis` 和 `test_dynamic_output`。最后一项运行真实示例，核对 CSV 与播放数据、观察节点、能量、固定边界和输出覆盖保护。
## 统一静力 / 动力建模入口

双击 `LaunchGUI.cmd`，在顶部“分析类型”选择：

- 静力分析：保留原来的分步加载和 Newton 求解。
- 动力：自由振动：忽略模型外力，从初始扰动及零速度开始运动。
- 动力：渐增加载：总力从 0 线性增长，至计算结束达到左侧填写的总力。
- 动力：突加载荷：从初始时刻持续施加左侧填写的总力。

四种模式共用矩形板几何、网格、材料、厚度和边界模板，调用同一个 `NonlinearFEM` 可执行程序。动力只支持固定零位移约束、无阻尼以及参考构形的集中质量。

动力参数面板：

| 参数 | 说明 |
| --- | --- |
| 参考密度 | 与几何、力和时间单位一致，必须为正 |
| 时间步 | 0 表示自动选取；正数表示手动指定 |
| 计算时长 | 从 t=0 开始计算，必须为正 |
| 每几步保存 | 保存间隔，必须是正整数；初始和最终状态始终保存 |
| 初始扰动幅值 | 悬臂模板为弯曲形位移；拉伸模板为轴向线性位移，可为 0；初速度为零 |
| 观察节点 | 从 0 编号；auto 选择右边界中部节点（奇数高度划分时取偏下节点） |

自动步长在初始扰动状态计算切线和行和上界，取安全系数 0.8 的估计值的 20%，再略微缩小以整步到达终止时刻。手动步长超过初始估计会被拒绝；大变形仍可能使后续稳定条件变化，因此自动选取不等于全过程稳定保证。用于估计的矩阵副本仅在整体非对称误差小于相对 1e-12 时消除舍入级非对称，不修改原内力或切线算法。

计算完成后自动打开时程与变形播放窗口，主窗口显示最终位移。修改输入会清除当前结果显示；运算期间修改参数时，不会将旧参数结果画成新模型的解。每次运行保留 `model.nfem`、`dynamic.cfg` 和结果，`analysis.txt` 记录实际使用步长及步数。

动力命令行入口：

```powershell
.\build\NonlinearFEM.exe --dynamic model.nfem dynamic.cfg results/new_run
```

`dynamic.cfg` 格式如下，8 个参数均必需：

```text
NFEM_DYNAMIC 1
density 1000
time_step 0
duration 1
output_every 10
initial_amplitude 0.001
observation_node 17
load free
initial_mode bending
```

`load` 可选 `free`、`ramp`、`step`；`initial_mode` 可选 `bending`、`axial`。节点 17 对应默认 8×2 矩形网格的右边界中点，其他网格需要自行选择。GUI 的“导出模型”在动力模式会同时保存 `.nfem` 和 `.dynamic.cfg`，可通过命令行复现。

为限制当前内存占用，单次动力分析最多 100 万步、约 200 万个保存节点状态；超限需减少时长、减少网格或增大保存间隔。当前使用全时程内存保存，尚不是流式输出。

`test_unified_analysis` 验证同一主程序的静力结果兼容、三种动力载荷、自动步长、默认悬臂板、终止时刻和配置快照；它不替代动力独立基准解和时间步收敛验证。
动力结果窗口默认使用 1 倍实际位移显示，打开新结果会重置为 1 倍。可点击‘恢复 1 倍’快速还原；放大后位移相对模型尺度过大时，界面会提示显示可能失真。该提示不是单元翻转的数值判定。

## 扩展边界模板

统一建模界面的“模型与边界条件”现在提供五种模板，均支持静力和三种动力载荷模式：

| 模板 | 约束 | 加载 | 动力 auto 观察点 |
| --- | --- | --- | --- |
| 悬臂板 | 左边 ux=uy=0 | 右边竖向均布总力 | 右边中部 |
| 拉伸板 | 左边 ux=0，左下角 uy=0 | 右边水平均布总力 | 右边中部 |
| 固定端轴向板 | 左边 ux=uy=0 | 右边水平均布总力 | 右边中部 |
| 剪切板 | 底边 ux=uy=0 | 顶边水平均布总力 | 顶边中部 |
| 简支梁 | 左下角 ux=uy=0，右下角 uy=0 | 顶边竖向均布总力 | 跨中、高度中部 |

总力按参考边长一致分配到节点，两端半权重；正水平力向右，正竖向力向上。水平总力设负可计算压缩路径，但不代表已经具备屈曲路径追踪能力。剪切板为底部固定、顶部受力的板，不是解析均匀纯剪切边界。简支梁的点支承附近存在局部效应，不宜直接将局部应力解释为网格收敛结果。

动力初始扰动按模板选择：悬臂采用二次弯曲形，轴向板采用线性轴向位移，剪切板采用沿高度线性水平位移，简支梁采用沿长度正弦弯曲形；固定自由度保持零。奇数网格划分没有精确中点节点时，auto 选偏左/偏下的相邻节点，不做观察点插值。自由振动仍忽略模型外力。

`test_gui_model` 已覆盖五种模板的载荷总和、约束、真实静力求解与整体反力平衡；`test_unified_analysis` 覆盖新增模板的动力求解与固定边界。这些验证不等同于每种模板都已完成独立基准解对照。
## AI 辅助开发声明

本项目在开发过程中使用了 OpenAI ChatGPT / Codex 辅助，包括算法讨论、部分代码与测试编写、代码检查、界面与可视化开发，以及文档整理。项目代码由作者与 AI 协作完成，并非全部由作者独立编写，也并非全部由 AI 生成。

作者负责理解、审查和维护项目代码。AI 提供的代码和解释可能存在错误；已有自动测试仅覆盖特定算例与检查条件，不构成对全部功能、数值精度或工程适用性的保证。用于研究结论或实际分析前，仍需进行人工审查、独立基准对照及适用性验证。
## 许可证

本项目采用 [MIT License](LICENSE)，版权署名为 zansy117。允许使用、修改、分发和商用，包括闭源复用，但须按许可证要求保留版权及许可声明。软件按现状提供，不作保证。

## 下载 Windows 便携版

首个预发布版为 **v0.1.0**，面向 Windows 10/11 x64。在 [GitHub Releases](https://github.com/zansy117/NonlinearFEM/releases/tag/v0.1.0) 下载 `NonlinearFEM-v0.1.0-windows-x64.zip`，完整解压到可写目录后，双击 `NonlinearFEM-Workbench.exe`。不需要另行安装 Python 或 C++ 编译工具。请保留 `_internal/` 和 `build/` 文件夹。

见 [快速入门](QUICKSTART.md)、[版本说明](RELEASE_NOTES.md) 和 [第三方许可](THIRD_PARTY_NOTICES.md)。这是未代码签名的学习研究预发布版，不是工程适用性认证。发行包附 SHA-256 校验文件；尚未在干净 Windows 虚拟机上验证。

### 从源码打包（Windows）

先构建并测试 Release 版本。GNU/MinGW 工具链可静态链接其运行库：

```powershell
cmake -S . -B build-release -G Ninja -DCMAKE_BUILD_TYPE=Release "-DCMAKE_EXE_LINKER_FLAGS=-static -static-libgcc -static-libstdc++"
cmake --build build-release
ctest --test-dir build-release --output-on-failure
python -m pip install -r packaging/requirements.txt
python tools/package_windows.py --solver-build build-release
```

以上命令要求 CMake、Ninja 和 GNU C++ 编译器已配置。其他工具链需要调整配置；脚本不会自动下载编译器。打包脚本生成 `dist/NonlinearFEM-v0.1.0-windows-x64/`，拒绝覆盖同名输出目录。它复制当前 Python 的运行环境，因此应使用包含 tkinter 的 Windows x64 Python。使用其他 Python 版本打包时，应复核并更新对应运行环境的第三方许可。

发行验证包括 Release 的 31 项测试，以及实际便携 EXE 在独立目录、移除开发工具 PATH 后的界面初始化、静力/动力求解与结果窗口加载检查。该验证不能代替干净系统兼容性测试。