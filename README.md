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

目前注册了 17 个测试程序：

| 范围 | 测试程序 |
| --- | --- |
| 运动学、应变、材料 | `test_kinematics`、`test_Gstrain`、`test_material` |
| Q4、组装、边界条件 | `test_q4`、`test_assembly`、`test_boundary` |
| Newton 求解 | `test_solver`、`test_solver_2x2` |
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