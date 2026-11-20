"""Verify postprocessing definitions, invariance and solver energy consistency."""
import csv
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'gui'))
from contour_fields import FieldResults,choices
import model,h8_model
from boundary_data import mesh_key

SOLVER=sys.argv.pop() if len(sys.argv)>1 and not sys.argv[-1].startswith('-') else None

class Contours(unittest.TestCase):
    def mesh(self,dim):
        module=model if dim==2 else h8_model
        p=dict(model.PRESETS['cantilever']) if dim==2 else dict(h8_model.DEFAULTS)
        p.update(nx=2,ny=1,length=2,height=1,young=1000,poisson=.25)
        if dim==3:p.update(nz=1,depth=1)
        data=module.geometry(p,'cantilever');return module,p,data[1],data[2]

    def test_affine_extension_and_shear(self):
        for dim in (2,3):
            module,p,nodes,elements=self.mesh(dim)
            for shear in (False,True):
                with self.subTest(dim=dim,shear=shear):
                    a=.03
                    u=[tuple((a*x[1] if shear else a*x[0]) if d==0 else 0 for d in range(dim)) for x in nodes]
                    result=FieldResults(nodes,elements,[u],1000,.25);row=result.recover(0)[0]
                    exx=0 if shear else a+.5*a*a;eyy=.5*a*a if shear else 0;gxy=a if shear else 0
                    self.assertAlmostEqual(row[0],exx,places=12);self.assertAlmostEqual(row[1],eyy,places=12);self.assertAlmostEqual(row[3],gxy,places=12)
                    if dim==2:expected=1000/(1-.25**2)*(exx+.25*eyy)
                    else:expected=400*(exx+eyy)+800*exx
                    self.assertAlmostEqual(row[6],expected,places=10);self.assertAlmostEqual(row[9],400*gxy,places=10)
                    self.assertEqual(len(result.values(0,'von Mises（第二PK）')),len(elements))
                    self.assertTrue(all(v>=0 for v in result.values(0,'von Mises（第二PK）')))

    def test_rigid_rotation_and_translation(self):
        for dim in (2,3):
            _,p,nodes,elements=self.mesh(dim);c=math.cos(.6);s=math.sin(.6)
            u=[tuple([c*x[0]-s*x[1]-x[0]+3,s*x[0]+c*x[1]-x[1]-2]+([.7] if dim==3 else [])) for x in nodes]
            results=FieldResults(nodes,elements,[u],1000,.25)
            for row in results.recover(0):self.assertLess(max(abs(v) for v in row[:13]),1e-10)

    def test_distorted_affine_patch_and_export(self):
        for dim in (2,3):
            _,p,nodes,elements=self.mesh(dim)
            nodes=[tuple([x[0]+.2*x[1],x[1]]+([x[2]+.15*x[0]] if dim==3 else [])) for x in nodes]
            u=[tuple(.02*x[d] for d in range(dim)) for x in nodes]
            result=FieldResults(nodes,elements,[u],1000,.25)
            for row in result.recover(0):
                for d in range(dim):self.assertAlmostEqual(row[d],.0202,places=12)
            with tempfile.TemporaryDirectory() as temp:
                path=Path(temp)/'fields.csv';result.export(path,0,.2)
                with path.open(encoding='utf-8-sig') as stream:rows=list(csv.DictReader(stream))
                self.assertEqual(len(rows),len(elements));self.assertEqual(float(rows[0]['time']),.2)
                self.assertIn('GL_engineering_xy',rows[0]);self.assertIn('von_mises_PK2',rows[0])

    @unittest.skipUnless(SOLVER,'Solver path not supplied')
    def test_q4_static_model_snapshot(self):
        module,p,nodes,elements=self.mesh(2);u=[(.01*x,.005*y) for x,y in nodes]
        custom={'mesh':mesh_key(nodes,elements),'constraints':[(n,d,v) for n,a in enumerate(u) for d,v in enumerate(a)],'loads':[]}
        with tempfile.TemporaryDirectory(prefix='nfem_fields_q4_') as temp:
            root=Path(temp);text=module.model_text(p,'cantilever',custom);(root/'input.nfem').write_text(text)
            result=subprocess.run([str(Path(SOLVER).resolve()),'input.nfem','solution'],cwd=root,capture_output=True,text=True,timeout=30)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertEqual((root/'solution/model.nfem').read_text(),text)
            with (root/'solution/nodes.csv').open() as stream:rows=list(csv.DictReader(stream))
            actual=[(float(r['ux']),float(r['uy'])) for r in rows]
            fields=FieldResults.from_folder(root/'solution',nodes,elements,[actual])
            self.assertIsNotNone(fields)
            for row in fields.recover(0):self.assertAlmostEqual(row[0],.01005,places=12)

    @unittest.skipUnless(SOLVER,'Solver path not supplied')
    def test_h8_energy_matches_existing_solver(self):
        module,p,nodes,elements=self.mesh(3);u=[(.01*x,.005*y,-.002*z) for x,y,z in nodes]
        custom={'mesh':mesh_key(nodes,elements),'constraints':[(n,d,v) for n,a in enumerate(u) for d,v in enumerate(a)],'loads':[]}
        with tempfile.TemporaryDirectory(prefix='nfem_fields_') as temp:
            root=Path(temp);(root/'input.nfem').write_text(module.model_text(p,'cantilever',custom))
            result=subprocess.run([str(Path(SOLVER).resolve()),'--h8','input.nfem','solution'],cwd=root,capture_output=True,text=True,timeout=30)
            self.assertEqual(result.returncode,0,result.stderr)
            with (root/'solution/history.csv').open() as stream:energy=float(next(csv.DictReader(stream))['strain_energy'])
            fields=FieldResults(nodes,elements,[u],1000,.25)
            recovered=sum(.5*sum(row[i]*row[6+i] for i in range(6))*row[13] for row in fields.recover(0))
            self.assertAlmostEqual(energy,recovered,places=10)

if __name__=='__main__':unittest.main()
