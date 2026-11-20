# v0.2.0 发布验证

日期：2026-11-20，Windows 11 x64。

- Release 构建：33/33 项 CTest 通过。
- Release 标准集：10/10 组通过，数值结果与源码/可执行文件指纹见本目录 report.json 和 REPORT.md。
- 实际 ZIP 完整解压到包含中文和空格的独立目录，清除 PYTHONPATH、PYTHONHOME、TCL_LIBRARY、TK_LIBRARY，并将 PATH 精简为 Windows System32 后检查。
- 冻结工作台通过 Q4/H8 自定义边界应用、静力、动力、全部云图分量、播放与手册文件检查。
- 冻结工作台的 --standard-cases 入口通过全部 10 组，不依赖外部 Python 或编译器 PATH。

这是开发电脑上的隔离路径验证，不代表已在干净虚拟机或所有 Windows 机器上测试。发布包未代码签名。最终软件包另附 BUILD_INFO.txt 记录源码提交与构建环境。
