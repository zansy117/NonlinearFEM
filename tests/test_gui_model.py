import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "gui"))
from model import PRESETS, geometry, model_text


class Models(unittest.TestCase):
    def test_boundary_and_force(self):
        for mode, values in PRESETS.items():
            p,nodes,elements,constraints,loads = geometry(values,mode)
            self.assertEqual(len(nodes),(p["nx"]+1)*(p["ny"]+1))
            self.assertEqual(len(elements),p["nx"]*p["ny"])
            self.assertAlmostEqual(sum(load[2] for load in loads),p["force"])
            expected_count = {"cantilever":2*(p["ny"]+1),"tension":p["ny"]+2,
                              "clamped_axial":2*(p["ny"]+1),"shear":2*(p["nx"]+1),"simply_supported":3}[mode]
            self.assertEqual(len(constraints),expected_count)
            self.assertEqual(len(set((n,d) for n,d,v in constraints)),len(constraints))
            if mode=="shear":
                self.assertTrue(all(nodes[n][1]==0 for n,d,v in constraints))
                self.assertTrue(all(nodes[n][1]==p["height"] and d==0 for n,d,v in loads))
            if mode=="simply_supported":
                self.assertEqual(set(constraints),{(0,0,0.),(0,1,0.),(p["nx"],1,0.)})
                self.assertTrue(all(nodes[n][1]==p["height"] and d==1 for n,d,v in loads))
            self.assertEqual(loads[0][2]*2,loads[1][2])
            for a,b,c,d in elements:
                self.assertGreater((nodes[b][0]-nodes[a][0])*(nodes[d][1]-nodes[a][1]),0)

    def test_invalid(self):
        for key,value in (("nx",1.5),("nx",0),("young",-1),("poisson",.5),
                          ("length",float("nan")),("ny",1000)):
            p=dict(PRESETS["tension"],**{key:value})
            with self.assertRaises(ValueError): model_text(p,"tension")

    def test_generated_model_solver(self):
        solver=Path(sys.argv[1]) if len(sys.argv)>1 else None
        if not solver or not solver.is_file():
            self.skipTest("Solver path not provided")
        for mode,p in PRESETS.items():
            with tempfile.TemporaryDirectory(prefix="nfem_gui_") as directory:
                folder=Path(directory)
                model=folder/"input.nfem"
                model.write_text(model_text(p,mode),encoding="utf-8")
                result=subprocess.run([str(solver),str(model),str(folder/"solution")],
                                      capture_output=True,text=True,timeout=60)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                import csv
                with (folder/"solution"/"nodes.csv").open() as stream:
                    rows=list(csv.DictReader(stream))
                _,_,_,fixed,loads=geometry(p,mode)
                for direction in (0,1):
                    component="reaction_x" if direction==0 else "reaction_y"
                    total=sum(value for node,d,value in loads if d==direction)
                    self.assertAlmostEqual(sum(float(row[component]) for row in rows),-total,delta=1e-4)
                for node,d,value in fixed:
                    self.assertAlmostEqual(float(rows[node]["ux" if d==0 else "uy"]),value,delta=1e-12)
                if mode=="tension":
                    for row in rows:
                        self.assertAlmostEqual(float(row["ux"]),.2*float(row["x"]),delta=1e-8)
                elif mode=="cantilever":
                    self.assertAlmostEqual(float(rows[17]["uy"]),-.57922940616,delta=1e-7)
                elif mode=="shear":
                    self.assertGreater(float(rows[-1]["ux"]),0)
                elif mode=="clamped_axial":
                    self.assertGreater(float(rows[-1]["ux"]),0)
                else:
                    self.assertLess(float(rows[(int(p["ny"])//2)*(int(p["nx"])+1)+int(p["nx"])//2]["uy"]),0)


if __name__=="__main__":
    unittest.main(argv=[sys.argv[0]])
