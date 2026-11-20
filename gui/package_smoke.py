"""Explicit release self-check of the frozen workbench and its solver handoff."""
import json
from pathlib import Path
import subprocess
import tempfile
import tkinter as tk
from runtime import application_root,VERSION

def check(report_path):
    report_path=Path(report_path);root=None
    try:
        from app import App
        from h8_app import H8App
        from boundary_editor import BoundaryEditor
        from boundary_data import face_nodes
        from contour_fields import choices
        from dynamics import Viewer as Q4Viewer
        from viewer3d import Viewer as H8Viewer
        import model,h8_model
        errors=[];root=tk.Tk();root.withdraw();root.report_callback_exception=lambda *args:errors.append(str(args))
        app=App(root);child=tk.Toplevel(root);child.withdraw();h8=H8App(child);root.update()
        solver=application_root()/'build/NonlinearFEM.exe';checks={}
        with tempfile.TemporaryDirectory(prefix='nfem_portable_') as directory:
            temporary=Path(directory)
            for dim,module in ((2,model),(3,h8_model)):
                values=dict(model.PRESETS['cantilever']) if dim==2 else dict(h8_model.DEFAULTS)
                values.update(nx=2,ny=1,young=1000)
                if dim==3:values['nz']=1
                data=module.geometry(values,'cantilever');nodes,elements,constraints,loads=data[1:5]
                applied=[];editor=BoundaryEditor(root,nodes,elements,constraints,loads,applied.append)
                editor.window.withdraw();editor.clear();editor.select(face_nodes(nodes,0,0));editor.assign('constraints')
                editor.select(face_nodes(nodes,0,1))
                for d,var in enumerate(editor.directions):var.set(d==1)
                editor.force.set('-.01');editor.assign('loads');editor.apply()
                if len(applied)!=1:raise RuntimeError('Boundary editor did not apply')
                folder=temporary/('q4' if dim==2 else 'h8');folder.mkdir()
                (folder/'model.nfem').write_text(module.model_text(values,'cantilever',applied[0]),encoding='utf-8')
                (folder/'dynamic.cfg').write_text('NFEM_DYNAMIC 1\ndensity 1000\ntime_step 0\nduration .2\noutput_every 1\ninitial_amplitude .001\nobservation_node '+str(len(nodes)-1)+'\nload free\ninitial_mode bending\n',encoding='utf-8')
                for dynamic in (False,True):
                    output='dynamic' if dynamic else 'static'
                    args=[str(solver)]+(['--dynamic' if dim==2 else '--dynamic-h8'] if dynamic else ['--h8'] if dim==3 else [])+['model.nfem']+(['dynamic.cfg'] if dynamic else [])+[output]
                    result=subprocess.run(args,cwd=folder,capture_output=True,text=True,timeout=60,creationflags=subprocess.CREATE_NO_WINDOW)
                    if result.returncode:raise RuntimeError(result.stdout+result.stderr)
                    window=tk.Toplevel(root);window.withdraw()
                    viewer=Q4Viewer(window) if dim==2 else H8Viewer(window)
                    viewer.load(folder/output)
                    if viewer.contour is None:raise RuntimeError('Missing frozen contour module or material snapshot')
                    for field in choices(dim):viewer.field.set(field);viewer.draw();root.update()
                    if dynamic:
                        viewer.toggle();root.update();viewer.stop()
                    if dim==2:viewer.close()
                    else:window.destroy()
                    checks[f'{dim}d_{output}']=True
                checks[f'{dim}d_custom_boundary']=True
            if errors:raise RuntimeError('; '.join(errors))
        for name in ('THEORY.md','USER.md','EXAMPLES.md'):
            if not (application_root()/'manuals'/name).is_file():raise RuntimeError('Missing manual '+name)
        h8.close();root.update();root.destroy();root=None
        report={'ok':True,'version':VERSION,'tkinter_ui':True,'manuals':True,'checks':checks}
        report_path.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');return 0
    except Exception as error:
        if root is not None:root.destroy()
        report_path.write_text(json.dumps({'ok':False,'error':str(error)},ensure_ascii=False),encoding='utf-8');return 1
