import csv
from pathlib import Path
import subprocess
import sys
import tempfile

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'gui'))
from model import PRESETS,model_text
from dynamic_config import configuration
from dynamics import load_results

def require(ok,message):
    if not ok:raise AssertionError(message)

if __name__=='__main__':
    solver=sys.argv[1]
    with tempfile.TemporaryDirectory(prefix='nfem_unified_') as directory:
        root=Path(directory)
        values=dict(PRESETS['tension'])
        model=root/'model.nfem';model.write_text(model_text(values,'tension'))
        controls=dict(density='1000',time_step='0',duration='.05',output_every='3',
                      initial_amplitude='.0001',observation_node='auto')
        results={}
        for index,analysis in enumerate(('动力：自由振动','动力：渐增加载','动力：突加载荷')):
            config=root/f'{index}.cfg';config.write_text(configuration(controls,values,'tension',analysis))
            output=root/f'run{index}'
            completed=subprocess.run([solver,'--dynamic',str(model),str(config),str(output)],capture_output=True,text=True,timeout=60)
            require(completed.returncode==0,completed.stdout+completed.stderr)
            nodes,elements,history,frames=load_results(output)
            require(len(nodes)==9 and len(elements)==4,'Model not reused')
            require(abs(history[-1]['time']-.05)<1e-12,'Requested duration not reached')
            require(abs(history[0]['ux']-.0001)<1e-14,'Initial perturbation not applied')
            require((output/'model.nfem').read_text()==model.read_text(),'Model snapshot missing')
            require((output/'dynamic.cfg').read_text()==config.read_text(),'Configuration snapshot missing')
            require((output/'analysis.txt').is_file(),'Actual time step metadata missing')
            require(all(frame[0]==(0,0) for frame in frames),'Fixed boundary moved')
            if index==0:
                e0=history[0]['total_energy']
                require(max(abs(h['total_energy']-e0)/e0 for h in history)<.05,'Automatic-step free-vibration energy')
            results[index]=history[-1]['ux']
        require(abs(results[1]-results[0])>1e-7 and abs(results[2]-results[1])>1e-7,'Loading modes not distinguished')
        # GUI default cantilever: includes cancellation-level stiffness entries.
        beam=root/'beam.nfem';beam.write_text(model_text(PRESETS['cantilever'],'cantilever'))
        beam_config=root/'beam.cfg'
        beam_config.write_text(configuration(dict(controls,duration='1',initial_amplitude='.001'),PRESETS['cantilever'],'cantilever','动力：自由振动'))
        completed=subprocess.run([solver,'--dynamic',str(beam),str(beam_config),str(root/'beam_output')],capture_output=True,text=True,timeout=60)
        require(completed.returncode==0,completed.stdout+completed.stderr)
        _,_,beam_history,_=load_results(root/'beam_output')
        require(abs(beam_history[-1]['time']-1)<1e-10,'Default cantilever final time')
        for mode in ("shear","clamped_axial","simply_supported"):
            specimen=root/(mode+'.nfem');specimen.write_text(model_text(PRESETS[mode],mode))
            cfg=root/(mode+'.cfg')
            cfg.write_text(configuration(dict(controls,duration='.1'),PRESETS[mode],mode,'动力：自由振动'))
            out=root/(mode+'_dynamic')
            completed=subprocess.run([solver,'--dynamic',str(specimen),str(cfg),str(out)],capture_output=True,text=True,timeout=60)
            require(completed.returncode==0,completed.stdout+completed.stderr)
            _,_,h,frames=load_results(out)
            require(abs(h[-1]['time']-.1)<1e-10,'New template duration')
            from model import geometry
            _,_,_,fixed,_=geometry(PRESETS[mode],mode)
            for frame in frames:
                for node,d,value in fixed:require(frame[node][d]==value,'New template fixed boundary moved')
        # Validate static backwards compatibility with the exact same executable/model.
        static=root/'static'
        completed=subprocess.run([solver,str(model),str(static)],capture_output=True,text=True,timeout=60)
        require(completed.returncode==0,completed.stderr)
        with (static/'nodes.csv').open() as stream:rows=list(csv.DictReader(stream))
        require(abs(float(rows[2]['ux'])-.2)<1e-8,'Static solver changed')
        bad=root/'bad.cfg'
        bad.write_text(configuration(dict(controls,time_step='100'),values,'tension','动力：自由振动'))
        completed=subprocess.run([solver,'--dynamic',str(model),str(bad),str(root/'bad_output')],capture_output=True,text=True,timeout=60)
        require(completed.returncode!=0 and not (root/'bad_output').exists(),'Excessive initial time step accepted')
        for key,value in [('density','0'),('duration','nan'),('output_every','1.2'),('observation_node','999')]:
            caught=False
            try:configuration(dict(controls,**{key:value}),values,'tension','动力：自由振动')
            except ValueError:caught=True
            require(caught,'Invalid GUI control accepted')
    print('PASS: unified static/dynamic executable, three load modes, automatic step, duration, snapshots and GUI controls.')
