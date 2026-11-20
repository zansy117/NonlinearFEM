# H8 三维工作台

从项目根目录运行 `LaunchGUI.cmd`，点击右上角“三维 H8 工作台”。原窗口继续提供 Q4 静力与动力分析；两者调用同一个 `build/NonlinearFEM.exe`。

## 使用

1. 选择悬臂体、轴向体或剪切体。网格按用户的 H8 局部节点顺序生成。
2. 用滑块或输入框设置三维尺寸、网格数、材料参数和加载面的总力。手动输入后点击“预览”。所有参数采用一致单位。
3. 选择静力分析、自由振动、渐增加载或突加载荷。动力参数在单独标签页中设置；观察节点自动选择在加载面中部附近。
4. 点击“运行计算”。输出保存在 `results/h8_时间戳/solution`，计算日志为上一层的 `solver.log`。
5. 左键拖动旋转，滚轮缩放。可选位移大小或 ux、uy、uz，输入变形倍率并按回车。默认倍率为 1，灰色虚线表示原网格。
6. 动力结果可播放或拖动时间条，下面显示观察节点三向位移和能量。视图尺寸在同一时程中固定，颜色范围可选自动、固定或手动。

位移按单元节点均值着色，也可选择第二 PK 应力和 GL 应变；定义与平均方式见 [CONTOURS.md](CONTOURS.md)。静力只保存最终构形，不能播放加载过程。可从“打开已有 H8 结果”重新载入 solution 文件夹。

## 输出文件

| 文件 | 内容 |
| --- | --- |
| nodes.csv / elements.csv | 原始三维网格与局部连接顺序 |
| frames.csv | 每个保存时刻、每个节点的 ux、uy、uz |
| history.csv | 观察节点的位移、速度、加速度以及动能、应变能、总能量 |
| reactions.csv | 静力约束自由度上的反力；其他分量写 0 |
| load_history.csv | 静力加载步、迭代次数、残差与重试次数 |
| model.nfem / dynamic.cfg / analysis.txt | 本次模型、动力配置及运行参数 |

受迫振动的动能与应变能之和不要求恒定；这里未累计外力功。

## 命令行与文本模型

```powershell
build\NonlinearFEM.exe --h8 examples\h8_cantilever.nfem results\h8_static
build\NonlinearFEM.exe --dynamic-h8 examples\h8_cantilever.nfem examples\h8_free.cfg results\h8_dynamic
```

输出目录必须尚不存在。建议在项目根目录使用相对路径和英文文件名。

模型头为 `NFEM3D 1`；材料行为 `material E nu`。nodes 为三坐标，elements 为八个从 0 开始的节点编号。loads 和 constraints 的记录格式为 `node direction value`，方向 0/1/2 对应 X/Y/Z。没有 Q4 的 thickness 参数。完整格式见示例。

H8 局部自然坐标顺序：0=(-,-,-)，1=(-,-,+)，2=(+,-,+)，3=(+,-,-)，4=(+,+,-)，5=(-,+,-)，6=(-,+,+)，7=(+,+,+)。不能直接套用其他软件的连接顺序。

## 当前范围

H8 使用现有三维 St. Venant–Kirchhoff 材料、2×2×2 积分、静力 Newton 加载与显式动力接口。入口暂限 600 节点，现有稠密刚度矩阵不适合大模型。自动时间步依据初始状态，不能保证后续大变形中的稳定性；动力约束仅支持固定零位移。

本次验证了编译、三个静力模型的输出读取和动力多帧接入，并运行原有 Q4 测试。H8 的独立数值精度、刚体运动、patch test、切线差分和能量收敛验证仍待开展。当前 GitHub v0.1.0 软件包不含此新增功能，需要使用更新后的本地源码。

自定义节点与边界面选择、约束及载荷设置见 [CUSTOM_BOUNDARIES.md](CUSTOM_BOUNDARIES.md)。
