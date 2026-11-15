# NonlinearFEM

使用 C++17 实现的二维 Q4 几何非线性有限元研究项目。当前重点是求解流程和数值验证，算例的网格、载荷及边界条件在测试代码中定义。

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

目前注册了 18 个测试程序：

| 范围 | 测试程序 |
| --- | --- |
| 运动学、应变、材料 | `test_kinematics`、`test_Gstrain`、`test_material` |
| Q4、组装、边界条件 | `test_q4`、`test_assembly`、`test_boundary` |
| Newton 求解 | `test_solver`、`test_solver_2x2` |
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

- 主程序 `NonlinearFEM` 目前只输出初始化信息；实际求解示例位于 `tests/`。
- 网格、材料和载荷需要在代码中设置，尚未提供通用的模型输入和结果可视化流程。
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