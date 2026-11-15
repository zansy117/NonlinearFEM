"""Rectangular Q4 model generation; no GUI dependency."""
import math

PRESETS = {
    "tension": dict(length=1.0, height=1.0, nx=2, ny=2, young=210000.0,
                    poisson=0.3, thickness=0.2, force=11088.0),
    "cantilever": dict(length=4.0, height=1.0, nx=8, ny=2, young=100000.0,
                       poisson=0.3, thickness=0.2, force=-50.0),
}

PRESETS.update({
    "clamped_axial": dict(length=2.,height=1.,nx=8,ny=4,young=100000.,poisson=.3,thickness=.2,force=100.),
    "shear": dict(length=1.,height=1.,nx=4,ny=4,young=100000.,poisson=.3,thickness=.2,force=100.),
    "simply_supported": dict(length=4.,height=1.,nx=8,ny=2,young=100000.,poisson=.3,thickness=.2,force=-50.),
})
TEMPLATES = {
    "cantilever": "悬臂板：左端固定，右端竖向力",
    "tension": "拉伸板：左边水平约束，右端水平力",
    "clamped_axial": "固定端轴向板：左端全固定，右端水平力",
    "shear": "剪切板：底边全固定，顶边水平力",
    "simply_supported": "简支梁：底角支承，顶边竖向均布力",
}
DESCRIPTIONS = {
    "cantilever": "左边全固定，右边均布竖向力；正值向上，负值向下。",
    "tension": "左边 ux=0，左下角 uy=0；右边水平力。允许竖向收缩；负力表示压缩。",
    "clamped_axial": "左边 ux=uy=0，右边水平力。固定端限制横向收缩；负力表示压缩。",
    "shear": "底边 ux=uy=0，顶边水平力；正值向右。不是周期或解析纯剪切边界。",
    "simply_supported": "左下角 ux=uy=0，右下角 uy=0；顶部竖向均布力。支点附近有局部效应。",
}

def observation_node(p, mode):
    if mode=="simply_supported":return (p['ny']//2)*(p['nx']+1)+p['nx']//2
    if mode=="shear":return p['ny']*(p['nx']+1)+p['nx']//2
    return (p['ny']//2)*(p['nx']+1)+p['nx']



def geometry(values, mode):
    if mode not in PRESETS:
        raise ValueError("请选择模型类型")
    p = {key: float(values[key]) for key in PRESETS[mode]}
    if not all(math.isfinite(v) for v in p.values()):
        raise ValueError("输入必须是有限数值")
    for key in ("length", "height", "thickness", "young"):
        if p[key] <= 0:
            raise ValueError("尺寸、厚度和杨氏模量必须大于零")
    if not -1 < p["poisson"] < 0.5:
        raise ValueError("泊松比必须在 -1 与 0.5 之间")
    for key in ("nx", "ny"):
        if p[key] != int(p[key]) or p[key] < 1:
            raise ValueError("网格划分必须是正整数")
        p[key] = int(p[key])
    nx, ny = p["nx"], p["ny"]
    if (nx + 1) * (ny + 1) > 600:
        raise ValueError("当前界面限制为 600 节点，请减少网格划分")
    nodes = [(p["length"] * i / nx, p["height"] * j / ny)
             for j in range(ny + 1) for i in range(nx + 1)]
    elements = []
    for j in range(ny):
        for i in range(nx):
            a = j * (nx + 1) + i
            elements.append((a, a + 1, a + nx + 2, a + nx + 1))
    constraints, loads = [], []
    if mode in ("cantilever", "tension", "clamped_axial"):
        for j in range(ny + 1):
            left=j*(nx+1)
            constraints.append((left,0,0.0))
            if mode!="tension" or j==0:constraints.append((left,1,0.0))
            weight=.5 if j in (0,ny) else 1.
            loads.append((left+nx,1 if mode=="cantilever" else 0,p['force']*weight/ny))
    elif mode=="shear":
        for i in range(nx+1):
            constraints.extend(((i,0,0.0),(i,1,0.0)))
            weight=.5 if i in (0,nx) else 1.
            loads.append((ny*(nx+1)+i,0,p['force']*weight/nx))
    else:
        constraints=[(0,0,0.0),(0,1,0.0),(nx,1,0.0)]
        for i in range(nx+1):
            weight=.5 if i in (0,nx) else 1.
            loads.append((ny*(nx+1)+i,1,p['force']*weight/nx))
    return p, nodes, elements, constraints, loads


def model_text(values, mode):
    p, nodes, elements, constraints, loads = geometry(values, mode)
    lines = ["NFEM 1", f'material {p["young"]:.17g} {p["poisson"]:.17g}',
             f'thickness {p["thickness"]:.17g}', f"nodes {len(nodes)}"]
    for records, next_name in ((nodes, "elements"), (elements, "loads"),
                               (loads, "constraints"), (constraints, None)):
        lines.extend(" ".join(format(v, ".17g") if isinstance(v, float) else str(v)
                              for v in row) for row in records)
        if next_name:
            following = {"elements": elements, "loads": loads, "constraints": constraints}
            lines.append(f"{next_name} {len(following[next_name])}")
    lines += ["newton 40 1e-9 1e-10", "stepping 0.25 0.0001 0.25 0.5 1.5 8 100", "end"]
    return "\n".join(lines) + "\n"
