"""Build a Windows x64 portable directory. Run with PyInstaller installed."""
import argparse
import importlib.metadata
import json
from pathlib import Path
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'gui'))
from runtime import VERSION

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--solver-build',type=Path,default=ROOT/'build-release')
    args=parser.parse_args()
    if sys.platform!='win32':parser.error('Windows packaging must run on Windows.')
    build=args.solver_build.resolve()
    standard=json.loads((ROOT/'examples/standard/cases.json').read_text(encoding='utf-8'))
    binaries=['NonlinearFEM.exe','dynamic_demo.exe']+sorted({name+'.exe' for group in standard['existing_case_groups'] for name in group['executables']})
    for name in binaries:
        if not (build/name).is_file():parser.error(f'Compile Release targets first: {build/name}')
    distribution=ROOT/'dist'/f'NonlinearFEM-v{VERSION}-windows-x64'
    if distribution.exists():parser.error(f'Output already exists: {distribution}')
    subprocess.run([sys.executable,'-m','PyInstaller','--noconfirm','--onedir','--windowed',
        '--name','NonlinearFEM-Workbench','--paths',str(ROOT/'gui'),'--paths',str(ROOT/'tools'),
        '--distpath',str(ROOT/'build-package/dist'),'--workpath',str(ROOT/'build-package/work'),
        '--specpath',str(ROOT/'build-package'),'--hidden-import','package_smoke','--hidden-import','run_standard_cases',str(ROOT/'gui/app.py')],check=True)
    shutil.copytree(ROOT/'build-package/dist/NonlinearFEM-Workbench',distribution)
    (distribution/'build').mkdir()
    for name in binaries:shutil.copy2(build/name,distribution/'build'/name)
    for name in ('examples','manuals','include','src','tests','tools','gui'):
        shutil.copytree(ROOT/name,distribution/name,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    for name in ('LICENSE','README.md','QUICKSTART.md','RELEASE_NOTES.md','THIRD_PARTY_NOTICES.md','CMakeLists.txt','LaunchGUI.cmd','CONTOURS.md','CUSTOM_BOUNDARIES.md','H8_WORKBENCH.md'):
        shutil.copy2(ROOT/name,distribution/name)
    shutil.copytree(ROOT/'packaging/licenses',distribution/'licenses')
    pyinstaller=importlib.metadata.distribution('pyinstaller')
    notice=next(f for f in pyinstaller.files if str(f).replace('\\','/').endswith('.dist-info/licenses/COPYING.txt'))
    shutil.copy2(pyinstaller.locate_file(notice),distribution/'licenses/PyInstaller-COPYING.txt')
    print(distribution)

if __name__=='__main__':main()
