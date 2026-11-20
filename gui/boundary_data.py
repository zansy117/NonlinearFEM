"""Shared node boundary data for Q4 and H8 model input."""
import math

def mesh_key(nodes,elements):
    return tuple(map(tuple,nodes)),tuple(map(tuple,elements))

def boundary_records(nodes,elements,custom,constraints,loads):
    if custom is None:return constraints,loads
    if custom['mesh']!=mesh_key(nodes,elements):
        raise ValueError('网格已改变：请打开自定义边界，重新设置或恢复示例边界。')
    dim=len(nodes[0]);output=[]
    for name in ('constraints','loads'):
        records={}
        for n,d,v in custom[name]:
            if not isinstance(n,int) or not isinstance(d,int) or not 0<=n<len(nodes) or not 0<=d<dim or not math.isfinite(v):
                raise ValueError('自定义边界存在无效节点、方向或数值')
            key=n,d
            if name=='constraints' and key in records:raise ValueError('存在重复位移约束')
            records[key]=records.get(key,0)+v if name=='loads' else v
            if not math.isfinite(records[key]):raise ValueError('节点力之和超出数值范围')
        output.append([(n,d,v) for (n,d),v in sorted(records.items())])
    return tuple(output)

def face_nodes(nodes,axis,side):
    extreme=(min if side==0 else max)(p[axis] for p in nodes)
    span=max(p[axis] for p in nodes)-min(p[axis] for p in nodes)
    return {i for i,p in enumerate(nodes) if abs(p[axis]-extreme)<=max(span*1e-10,1e-14)}

def parse_nodes(text,count):
    selected=set()
    for item in text.replace('，',',').split(','):
        item=item.strip()
        if not item:continue
        parts=item.split('-')
        if len(parts)>2:raise ValueError('节点格式示例：0,2,5-8')
        a=int(parts[0]);b=int(parts[-1])
        if not 0<=a<=b<count:raise ValueError('节点编号超出范围或区间顺序错误')
        selected.update(range(a,b+1))
    return selected
