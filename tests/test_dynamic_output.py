import csv
from pathlib import Path
import subprocess
import sys
import tempfile

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"gui"))
from dynamics import load_results

def require(condition,message):
    if not condition:raise AssertionError(message)

if __name__=="__main__":
    with tempfile.TemporaryDirectory(prefix="nfem_dynamic_output_") as temporary:
        folder=Path(temporary)/"run"
        completed=subprocess.run([sys.argv[1],str(folder)],capture_output=True,text=True,timeout=60)
        require(completed.returncode==0,completed.stdout+completed.stderr)
        nodes,elements,history,frames=load_results(folder)
        require(len(nodes)==15 and len(elements)==8 and len(history)==401,"Output dimensions")
        require(abs(history[0]['time'])<1e-14 and abs(history[-1]['time']-2)<1e-10,"Output time interval")
        require(abs(history[0]['uy']-1e-4)<1e-14,"Initial observation displacement")
        e0=history[0]['total_energy']
        require(e0>0,"Initial energy")
        for i,row in enumerate(history):
            require(int(row['node'])==9 and int(row['frame'])==i,"Observation metadata")
            require(abs(row['total_energy']-row['kinetic_energy']-row['strain_energy'])<1e-15,"Energy columns")
            require(abs(row['total_energy']-e0)/e0<1e-4,"Energy drift")
            require(abs(row['ux']-frames[i][9][0])<1e-14 and abs(row['uy']-frames[i][9][1])<1e-14,"Frame versus history mismatch")
            for node in (0,5,10):require(frames[i][node]==(0,0),"Fixed boundary moved")
        before=(folder/"history.csv").read_bytes()
        again=subprocess.run([sys.argv[1],str(folder)],capture_output=True,timeout=60)
        require(again.returncode!=0 and (folder/"history.csv").read_bytes()==before,"Existing run overwritten")
        (folder/"frames.csv").write_text("frame,time,node,ux,uy\n",encoding="utf-8")
        rejected=False
        try:load_results(folder)
        except ValueError:rejected=True
        require(rejected,"Incomplete frames accepted")
    print("PASS: real dynamic solve, CSV state/energy, frame consistency, input protection and viewer loading.")
