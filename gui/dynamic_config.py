"""Validated dynamic controls shared by the GUI and tests."""
import math
from model import geometry, observation_node

def configuration(values, model, mode, analysis):
    p,nodes,*_=geometry(model,mode)
    try:
        numbers={key:float(values[key]) for key in ("density","time_step","duration","output_every","initial_amplitude")}
        if not all(math.isfinite(v) for v in numbers.values()):raise ValueError()
        observer=values["observation_node"]
        observer=observation_node(p,mode) if observer.strip()=="auto" else float(observer)
        if not math.isfinite(observer) or observer!=int(observer) or not 0<=observer<len(nodes):raise ValueError()
        if numbers['density']<=0 or numbers['time_step']<0 or numbers['duration']<=0:raise ValueError()
        every=numbers['output_every']
        if every<1 or every!=int(every) or every>1e6:raise ValueError()
    except (ValueError,KeyError,TypeError) as error:
        raise ValueError("请检查动力参数：密度、时长为正；步长非负；保存间隔为正整数；观察节点有效。") from error
    loading={"动力：自由振动":"free","动力：渐增加载":"ramp","动力：突加载荷":"step"}
    if analysis not in loading:raise ValueError("无效的动力分析模式")
    lines=["NFEM_DYNAMIC 1"]
    lines += [f"{key} {value:.17g}" for key,value in numbers.items()]
    lines += [f"observation_node {int(observer)}",f"load {loading[analysis]}",f"initial_mode { {'cantilever':'bending','tension':'axial','clamped_axial':'axial','shear':'shear','simply_supported':'sine_bending'}[mode] }"]
    return "\n".join(lines)+"\n"
