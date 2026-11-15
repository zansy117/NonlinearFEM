# Third-party notices

The project source is distributed under the MIT License in `LICENSE`.
The portable binary also contains third-party runtime components with their own terms.
No third-party runtime is relicensed under the project's MIT License.

| Component | Notice supplied in this package |
| --- | --- |
| CPython 3.13 runtime and standard library | `licenses/Python-LICENSE.txt` (includes additional incorporated-software notices) |
| Tcl 8.6 | `licenses/Tcl-license.terms` |
| Tk 8.6 | `licenses/Tk-license.terms` |
| PyInstaller bootloader | `licenses/PyInstaller-COPYING.txt`, including the bootloader exception |
| GCC/libgcc/libstdc++ runtime used by the statically linked solver | `licenses/GPL-3.0.txt` and `licenses/GCC-Runtime-Exception-3.1.txt` |
| MinGW-w64 runtime and winpthreads, where incorporated by the toolchain | `licenses/MinGW-w64-COPYING.txt` and `licenses/winpthreads-COPYING.txt` |

The GCC Runtime Library Exception and PyInstaller bootloader exception are included in their full text. The project source remains under its own MIT License.

Upstream sources: https://www.python.org/ ; https://www.tcl-lang.org/ ; https://pyinstaller.org/ ; https://gcc.gnu.org/ ; https://www.mingw-w64.org/

This release was built with Python 3.13, PyInstaller 6.22.3 and a GNU/MinGW-w64 C++ toolchain. Build helpers require PyInstaller only during packaging; users of the portable distribution do not need to install it.
