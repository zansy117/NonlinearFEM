"""H8 rectangular mesh; local node order follows the author's convention."""
import math

TEMPLATES={"cantilever":"悬臂体：左端固定，右面竖向力",
           "axial":"轴向体：允许泊松收缩，右面轴向力",
           "shear":"剪切体：底面固定，顶面水平力"}
DEFAULTS=dict(length=4,height=1,depth=1,nx=4,ny=2,nz=2,young=1000,poisson=.3,force=-.05)

def geometry(values,mode,custom=None):
    p={k:float(values[k]) for k in DEFAULTS}
    if mode not in TEMPLATES or not all(math.isfinite(v) for v in p.values()):
        raise ValueError("模型类型或数值无效")
    if any(p[k]<=0 for k in ('length','height','depth','young')) or not -1<p['poisson']<.5:
        raise ValueError("尺寸与杨氏模量须为正，泊松比在 (-1,0.5) 内")
    for key in ('nx','ny','nz'):
        if p[key]<1 or p[key]!=int(p[key]):raise ValueError("网格划分必须为正整数")
        p[key]=int(p[key])
    nx,ny,nz=p['nx'],p['ny'],p['nz']
    if (nx+1)*(ny+1)*(nz+1)>600:raise ValueError("当前 H8 入口限制 600 节点，请先使用小网格")
    node=lambda i,j,k:(k*(ny+1)+j)*(nx+1)+i
    nodes=[(p['length']*i/nx,p['height']*j/ny,p['depth']*k/nz)
           for k in range(nz+1) for j in range(ny+1) for i in range(nx+1)]
    elements=[]
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                elements.append((node(i,j,k),node(i,j,k+1),node(i+1,j,k+1),node(i+1,j,k),
                                 node(i+1,j+1,k),node(i,j+1,k),node(i,j+1,k+1),node(i+1,j+1,k+1)))
    constraints=[];loads=[]
    if mode in ('cantilever','axial'):
        for k in range(nz+1):
            for j in range(ny+1):
                n=node(0,j,k)
                for d in (range(3) if mode=='cantilever' else (0,)):constraints.append((n,d,0.))
                weight=(.5 if k in (0,nz) else 1)*(.5 if j in (0,ny) else 1)/(ny*nz)
                loads.append((node(nx,j,k),1 if mode=='cantilever' else 0,p['force']*weight))
        if mode=='axial':constraints.extend(((node(0,0,0),1,0.),(node(0,0,0),2,0.),(node(0,ny,0),2,0.)))
        observer=node(nx,ny//2,nz//2)
    else:
        for k in range(nz+1):
            for i in range(nx+1):
                constraints.extend((node(i,0,k),d,0.) for d in range(3))
                weight=(.5 if k in (0,nz) else 1)*(.5 if i in (0,nx) else 1)/(nx*nz)
                loads.append((node(i,ny,k),0,p['force']*weight))
        observer=node(nx//2,ny,nz//2)
    from boundary_data import boundary_records
    constraints,loads=boundary_records(nodes,elements,custom,constraints,loads)
    return p,nodes,elements,constraints,loads,observer

def model_text(values,mode,custom=None):
    p,nodes,elements,constraints,loads,observer=geometry(values,mode,custom)
    lines=['NFEM3D 1',f"material {p['young']:.17g} {p['poisson']:.17g}"]
    for key,records in (('nodes',nodes),('elements',elements),('loads',loads),('constraints',constraints)):
        lines.append(f'{key} {len(records)}')
        lines.extend(' '.join(format(v,'.17g') if isinstance(v,float) else str(v) for v in row) for row in records)
    lines+=['newton 40 1e-9 1e-10','stepping 0.25 0.0001 0.25 0.5 1.5 8 100','end']
    return '\n'.join(lines)+'\n'
