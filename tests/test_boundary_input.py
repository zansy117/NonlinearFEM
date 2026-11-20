"""Custom boundary IO regression; uses existing solvers, no new element algorithms."""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'gui'))
from boundary_data import mesh_key,face_nodes,parse_nodes
import model,h8_model

SOLVER=sys.argv.pop() if len(sys.argv)>1 and not sys.argv[-1].startswith('-') else None

class BoundaryInput(unittest.TestCase):
    def test_node_ranges(self):
        self.assertEqual(parse_nodes('0,2,4-6',8),{0,2,4,5,6})
        for bad in ('-1','8','5-2','0-9000000'):
            with self.assertRaises(ValueError):parse_nodes(bad,8)

    def test_export_and_actual_solve(self):
        for module,dim in ((model,2),(h8_model,3)):
            with self.subTest(dimension=dim):
                values=dict(model.PRESETS['cantilever']) if dim==2 else dict(h8_model.DEFAULTS)
                values.update(nx=2,ny=1)
                if dim==3:values['nz']=1
                p,nodes,elements,*_=module.geometry(values,'cantilever')
                fixed=face_nodes(nodes,1,1);loaded=face_nodes(nodes,1,0)
                custom={'mesh':mesh_key(nodes,elements),
                    'constraints':[(n,d,0.) for n in fixed for d in range(dim)],
                    'loads':[(n,0,.01/len(loaded)) for n in loaded]}
                records=module.geometry(values,'cantilever',custom)
                self.assertEqual(records[3],sorted(custom['constraints']))
                self.assertAlmostEqual(sum(v for n,d,v in records[4]),.01)
                with self.assertRaises(ValueError):module.model_text(dict(values,nx=3),'cantilever',custom)
                if not SOLVER:continue
                with tempfile.TemporaryDirectory(prefix='nfem_boundary_') as temp:
                    root=Path(temp);(root/'input.nfem').write_text(module.model_text(values,'cantilever',custom))
                    args=[str(Path(SOLVER).resolve())]+(['--h8'] if dim==3 else [])+['input.nfem','solution']
                    result=subprocess.run(args,cwd=root,capture_output=True,text=True,timeout=30)
                    self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                    with (root/'solution'/'nodes.csv').open() as stream:self.assertEqual(len(list(csv.DictReader(stream))),len(nodes))
                    if dim==3:
                        with (root/'solution'/'reactions.csv').open() as stream:rows=list(csv.DictReader(stream))
                        self.assertAlmostEqual(sum(float(r['rx']) for r in rows),-.01,places=7)
                    else:
                        with (root/'solution'/'nodes.csv').open() as stream:rows=list(csv.DictReader(stream))
                        for n in fixed:self.assertAlmostEqual(float(rows[n]['ux']),0,places=12)

if __name__=='__main__':unittest.main()
