"""Recover solver-defined GL strain and second PK stress from saved nodal results.

Reference-volume weighted Gauss averages; no nodal smoothing/extrapolation.
Shear strain fields are engineering shear = 2 E_ij, as in the solver.
"""
import csv
import itertools
import math
from collections import OrderedDict
from pathlib import Path

SIGNS={2:((-1,-1),(1,-1),(1,1),(-1,1)),
       3:((-1,-1,-1),(-1,-1,1),(1,-1,1),(1,-1,-1),(1,1,-1),(-1,1,-1),(-1,1,1),(1,1,1))}
PAIRS=((0,0),(1,1),(2,2),(0,1),(0,2),(1,2))
COMPONENTS=('xx','yy','zz','xy','xz','yz')
DISPLACEMENT={'位移大小':'mag','水平位移 ux':'ux','竖向位移 uy':'uy','ux':'ux','uy':'uy','uz':'uz'}

def choices(dim):
    indices=(0,1,3) if dim==2 else range(6)
    return tuple(['位移大小','ux','uy']+(['uz'] if dim==3 else [])+['von Mises（第二PK）']+
                 ['第二PK S'+COMPONENTS[i] for i in indices]+
                 [('GL E' if i<3 else 'GL 工程剪应变 g')+COMPONENTS[i] for i in indices])

def inverse(matrix):
    n=len(matrix);a=[list(row)+[float(i==j) for j in range(n)] for i,row in enumerate(matrix)]
    scale=max(abs(v) for row in matrix for v in row)
    determinant=1.
    for j in range(n):
        p=max(range(j,n),key=lambda i:abs(a[i][j]))
        if abs(a[p][j])<=1e-13*scale:raise ValueError('后处理：退化单元 Jacobian')
        if p!=j:a[p],a[j]=a[j],a[p];determinant=-determinant
        pivot=a[j][j];determinant*=pivot;a[j]=[v/pivot for v in a[j]]
        for i in range(n):
            if i!=j:
                factor=a[i][j];a[i]=[a[i][k]-factor*a[j][k] for k in range(2*n)]
    if not math.isfinite(determinant) or determinant<=0:raise ValueError('后处理：单元参考构形方向错误')
    return [row[n:] for row in a],determinant

def model_info(directory):
    directory=Path(directory)
    candidates=(directory/'model.nfem',directory.parent/'model.nfem',directory.parent/'input.nfem')
    for path in candidates:
        if not path.is_file():continue
        material=None;elements=None
        lines=[row for line in path.read_text(encoding='utf-8-sig').splitlines() if (row:=line.split('#')[0].strip().split())]
        for i,row in enumerate(lines):
            if row and row[0]=='material':material=tuple(map(float,row[1:]))
            if row and row[0]=='elements':elements=[tuple(map(int,r)) for r in lines[i+1:i+1+int(row[1])]]
        if material is not None and elements is not None:return material,elements
    return None,None

class FieldResults:
    def __init__(self,nodes,elements,frames,young,poisson):
        self.nodes=nodes;self.elements=elements;self.frames=frames;self.dim=len(nodes[0]);self.young=young;self.poisson=poisson
        if self.dim not in SIGNS or not math.isfinite(young) or young<=0 or not math.isfinite(poisson) or not -1<poisson<.5:
            raise ValueError('后处理：材料或维数无效')
        self.geometry=[];self.cache=OrderedDict()
        g=1/math.sqrt(3);dim=self.dim
        for element in elements:
            if len(element)!=2**dim or any(n<0 or n>=len(nodes) for n in element):raise ValueError('后处理：单元连接无效')
            points=[]
            for xi in itertools.product((-g,g),repeat=dim):
                derivatives=[[sign[d]*math.prod(1+sign[k]*xi[k] for k in range(dim) if k!=d)/2**dim for d in range(dim)] for sign in SIGNS[dim]]
                jac=[[sum(nodes[element[n]][i]*derivatives[n][j] for n in range(len(element))) for j in range(dim)] for i in range(dim)]
                inv,det=inverse(jac)
                gradients=[[sum(derivatives[n][k]*inv[k][j] for k in range(dim)) for j in range(dim)] for n in range(len(element))]
                points.append((det,gradients))
            self.geometry.append(points)

    @classmethod
    def from_folder(cls,folder,nodes,elements,frames):
        material,_=model_info(folder)
        if material and len(material)!=2:raise ValueError('模型材料信息必须包含 E 和 nu')
        return cls(nodes,elements,frames,*material) if material else None

    def recover(self,index):
        if index in self.cache:self.cache.move_to_end(index);return self.cache[index]
        u=self.frames[index];dim=self.dim;young=self.young;nu=self.poisson;mu=young/(2*(1+nu));rows=[]
        if len(u)!=len(self.nodes):raise ValueError('后处理：位移与网格不匹配')
        for element,points in zip(self.elements,self.geometry):
            total=[0.]*13;volume=0.
            for weight,gradients in points:
                f=[[float(i==j)+sum(u[n][i]*gradients[a][j] for a,n in enumerate(element)) for j in range(dim)] for i in range(dim)]
                e=[(sum(f[k][i]*f[k][j] for k in range(dim))-float(i==j))*(.5 if i==j else 1.) if i<dim and j<dim else 0. for i,j in PAIRS]
                if dim==2:
                    factor=young/(1-nu*nu);s=[factor*(e[0]+nu*e[1]),factor*(e[1]+nu*e[0]),0.,mu*e[3],0.,0.]
                else:
                    lam=young*nu/((1+nu)*(1-2*nu));trace=sum(e[:3])
                    s=[lam*trace+2*mu*v for v in e[:3]]+[mu*v for v in e[3:]]
                vm=math.sqrt(.5*((s[0]-s[1])**2+(s[1]-s[2])**2+(s[2]-s[0])**2)+3*sum(v*v for v in s[3:]))
                for k,v in enumerate(e+s+[vm]):total[k]+=weight*v
                volume+=weight
            row=[v/volume for v in total]+[volume]
            if not all(math.isfinite(v) for v in row):raise ValueError('后处理结果超出有限数值范围')
            rows.append(row)
        self.cache[index]=rows
        if len(self.cache)>8:self.cache.popitem(last=False)
        return rows

    def values(self,index,field):
        if field in DISPLACEMENT:
            key=DISPLACEMENT[field];u=self.frames[index]
            values=[math.hypot(*v) if key=='mag' else v['xyz'.index(key[-1])] for v in u]
            return [sum(values[n] for n in e)/len(e) for e in self.elements]
        columns={'von Mises（第二PK）':12}
        for i,c in enumerate(COMPONENTS):
            columns['第二PK S'+c]=6+i;columns[('GL E' if i<3 else 'GL 工程剪应变 g')+c]=i
        return [row[columns[field]] for row in self.recover(index)]

    def export(self,path,index,time=0,frame_id=None):
        with Path(path).open('w',encoding='utf-8-sig',newline='') as stream:
            indices=(0,1,3) if self.dim==2 else tuple(range(6))
            writer=csv.writer(stream);writer.writerow(['frame','time','element']+['GL_'+COMPONENTS[i] if i<3 else 'GL_engineering_'+COMPONENTS[i] for i in indices]+['PK2_'+COMPONENTS[i] for i in indices]+['von_mises_PK2','reference_measure'])
            for i,row in enumerate(self.recover(index)):writer.writerow([index if frame_id is None else frame_id,time,i]+[row[j] for j in indices]+[row[6+j] for j in indices]+row[12:])
