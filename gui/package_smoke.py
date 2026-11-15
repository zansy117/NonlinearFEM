"""Packaged-app smoke check, invoked explicitly by release validation."""
import json
from pathlib import Path
import subprocess
import tempfile
import tkinter as tk

from runtime import application_root, VERSION
from model import PRESETS, model_text
from dynamic_config import configuration
from dynamics import load_results, Viewer

def check(report_path):
    report_path=Path(report_path)
    root=None
    try:
        from app import App
        root=tk.Tk();root.withdraw()
        app=App(root);root.update()
        with tempfile.TemporaryDirectory(prefix="nfem_portable_check_") as temporary:
            folder=Path(temporary)
            model=folder/"model.nfem"
            model.write_text(model_text(PRESETS['tension'],'tension'),encoding='utf-8')
            solver=application_root()/"build"/"NonlinearFEM.exe"
            static=folder/"static"
            command=[str(solver),str(model),str(static)]
            completed=subprocess.run(command,capture_output=True,text=True,timeout=60,
                creationflags=subprocess.CREATE_NO_WINDOW)
            if completed.returncode:raise RuntimeError(completed.stdout+completed.stderr)
            if not (static/'results.html').is_file():raise RuntimeError("Missing static viewer")
            cfg=folder/'dynamic.cfg'
            cfg.write_text(configuration(dict(density='1000',time_step='0',duration='.03',
                output_every='1',initial_amplitude='.0001',observation_node='auto'),
                PRESETS['tension'],'tension','动力：自由振动'),encoding='utf-8')
            output=folder/'dynamic'
            completed=subprocess.run([str(solver),'--dynamic',str(model),str(cfg),str(output)],
                capture_output=True,text=True,timeout=60,creationflags=subprocess.CREATE_NO_WINDOW)
            if completed.returncode:raise RuntimeError(completed.stdout+completed.stderr)
            data=load_results(output)
            if abs(data[2][-1]['time']-.03)>1e-10:raise RuntimeError("Incorrect final time")
            child=tk.Toplevel(root);child.withdraw()
            viewer=Viewer(child);viewer.load(output);root.update()
            if viewer.data is None:raise RuntimeError("Failed to initialize dynamic viewer")
            viewer.close()
            report={'ok':True,'version':VERSION,'application_root':str(application_root()),
                    'static_solve':True,'dynamic_solve':True,'tkinter_ui':True,'frames':len(data[2])}
        root.destroy();root=None
        report_path.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
        return 0
    except Exception as error:
        if root is not None:root.destroy()
        report_path.write_text(json.dumps({'ok':False,'error':str(error)},ensure_ascii=False),encoding='utf-8')
        return 1
