"""Build a Windows x64 portable directory. Run with PyInstaller installed."""
import argparse
import importlib.metadata
from pathlib import Path
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
VERSION='0.1.0'

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--solver-build',type=Path,default=ROOT/'build-release')
    args=parser.parse_args()
    if sys.platform!='win32':parser.error('Windows packaging must run on Windows.')
    build=args.solver_build.resolve()
    for name in ('NonlinearFEM.exe','dynamic_demo.exe'):
        if not (build/name).is_file():parser.error(f'Compile Release targets first: {build/name}')
    distribution=ROOT/'dist'/f'NonlinearFEM-v{VERSION}-windows-x64'
    if distribution.exists():parser.error(f'Output already exists: {distribution}')
    subprocess.run([sys.executable,'-m','PyInstaller','--noconfirm','--onedir','--windowed',
        '--name','NonlinearFEM-Workbench','--paths',str(ROOT/'gui'),
        '--distpath',str(ROOT/'build-package/dist'),'--workpath',str(ROOT/'build-package/work'),
        '--specpath',str(ROOT/'build-package'),'--hidden-import','package_smoke',str(ROOT/'gui/app.py')],check=True)
    shutil.copytree(ROOT/'build-package/dist/NonlinearFEM-Workbench',distribution)
    (distribution/'build').mkdir()
    for name in ('NonlinearFEM.exe','dynamic_demo.exe'):shutil.copy2(build/name,distribution/'build'/name)
    shutil.copytree(ROOT/'examples',distribution/'examples',ignore=shutil.ignore_patterns('*.cpp'))
    for name in ('LICENSE','README.md','QUICKSTART.md','RELEASE_NOTES.md','THIRD_PARTY_NOTICES.md'):
        shutil.copy2(ROOT/name,distribution/name)
    shutil.copytree(ROOT/'packaging/licenses',distribution/'licenses')
    pyinstaller=importlib.metadata.distribution('pyinstaller')
    notice=next(f for f in pyinstaller.files if str(f).replace('\\','/').endswith('.dist-info/licenses/COPYING.txt'))
    shutil.copy2(pyinstaller.locate_file(notice),distribution/'licenses/PyInstaller-COPYING.txt')
    print(distribution)

if __name__=='__main__':main()
