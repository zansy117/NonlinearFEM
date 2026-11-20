"""Local dynamic result viewer, using only Python's standard library."""
import csv
from datetime import datetime
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from runtime import application_root, VERSION
from contour_fields import FieldResults,choices,model_info,DISPLACEMENT
from contour_controls import ContourControls,color
ROOT = application_root()


def load_results(directory):
    directory = Path(directory)
    if not (directory/'frames.csv').exists():
        with (directory/'nodes.csv').open(encoding='utf-8',newline='') as stream:rows=list(csv.DictReader(stream))
        material,elements=model_info(directory)
        if not rows or not elements:raise ValueError('静力云图需要原模型文件（model.nfem）')
        if 'z' in rows[0] or any(len(e)!=4 or any(n<0 or n>=len(rows) for n in e) for e in elements):raise ValueError('请选择 Q4 结果；三维结果请使用 H8 查看器')
        nodes=[(float(r['x']),float(r['y'])) for r in rows];frame=[(float(r['ux']),float(r['uy'])) for r in rows]
        if not all(math.isfinite(v) for row in nodes+frame for v in row):raise ValueError('结果包含非有限值')
        history=[dict(time=0,node=0,ux=frame[0][0],uy=frame[0][1],kinetic_energy=0,strain_energy=0,total_energy=0)]
        return nodes,elements,history,[frame]
    def read(name):
        with (directory / name).open(encoding="utf-8", newline="") as stream:
            return list(csv.DictReader(stream))
    node_rows=read('nodes.csv')
    if node_rows and 'z' in node_rows[0]:raise ValueError('三维结果请使用 H8 查看器')
    nodes = [(float(r["x"]), float(r["y"])) for r in node_rows]
    elements = [tuple(int(r[f"n{i}"]) for i in range(4)) for r in read("elements.csv")]
    history = [{k: float(v) for k, v in r.items()} for r in read("history.csv")]
    if not nodes or not history or not elements:
        raise ValueError("结果为空")
    for e in elements:
        if any(n < 0 or n >= len(nodes) for n in e):
            raise ValueError("单元编号无效")
    frames = [[None] * len(nodes) for _ in history]
    for r in read("frames.csv"):
        frame, node = int(r["frame"]), int(r["node"])
        if not 0 <= frame < len(frames) or not 0 <= node < len(nodes):
            raise ValueError("帧或节点编号无效")
        if frames[frame][node] is not None:
            raise ValueError("重复的节点记录")
        frames[frame][node] = (float(r["ux"]), float(r["uy"]))
        if abs(float(r["time"]) - history[frame]["time"]) > 1e-12:
            raise ValueError("帧时间不匹配")
    if any(u is None for frame in frames for u in frame):
        raise ValueError("帧记录不完整")
    numbers = [v for n in nodes for v in n]
    numbers += [v for frame in frames for u in frame for v in u]
    numbers += [v for row in history for v in row.values()]
    if not all(math.isfinite(v) for v in numbers):
        raise ValueError("结果包含非有限值")
    if any(history[i]["time"] <= history[i-1]["time"] for i in range(1, len(history))):
        raise ValueError("记录时间未递增")
    return nodes, elements, history, frames


class Viewer:
    def __init__(self, window):
        self.window = window
        window.title("NonlinearFEM · Q4 云图与时程")
        window.geometry("1160x800+30+20")
        window.minsize(920, 700)
        self.data = None
        self.contour = None
        self.process = None
        self.playing = False
        self.timer = None
        panel = ttk.Frame(window, padding=18)
        panel.pack(fill="both", expand=True)
        ttk.Label(panel, text="Q4 结果云图 / 时程", font=("Microsoft YaHei UI", 16, "bold")).pack(anchor="w")
        self.description = tk.StringVar(value="Q4 动力分析 · 单位沿用模型 · 可打开统一建模窗口生成的结果")
        ttk.Label(panel, textvariable=self.description).pack(anchor="w", pady=7)
        buttons = ttk.Frame(panel)
        buttons.pack(fill="x")
        self.run_button = ttk.Button(buttons, text="计算并查看示例", command=self.run)
        self.run_button.pack(side="left")
        ttk.Button(buttons, text="打开已有动力结果", command=self.choose).pack(side="left", padx=8)
        self.status = tk.StringVar(value="点击计算，生成位移、能量和变形时程；也可以打开已导出的结果目录。")
        ttk.Label(panel, textvariable=self.status, wraplength=1080).pack(anchor="w", pady=8)
        self.mesh = tk.Canvas(panel, background="#f6f9fc", highlightthickness=0, height=230)
        self.mesh.pack(fill="both", expand=True)
        row=ttk.Frame(panel);row.pack(fill='x');self.field=tk.StringVar(value='位移大小')
        self.selector=ttk.Combobox(row,textvariable=self.field,state='readonly',width=28,values=('位移大小','ux','uy'));self.selector.pack(side='left');self.selector.bind('<<ComboboxSelected>>',self.draw)
        self.contour_controls=ContourControls(panel,self.draw,self.export_fields)
        controls = ttk.Frame(panel)
        controls.pack(fill="x", pady=8)
        self.play_button = ttk.Button(controls, text="播放", command=self.toggle, state="disabled")
        self.play_button.pack(side="left")
        self.frame = tk.DoubleVar(value=0)
        self.timeline = ttk.Scale(controls, from_=0, to=1, variable=self.frame, command=self.draw)
        self.timeline.pack(side="left", fill="x", expand=True, padx=10)
        self.time_label = tk.StringVar(value="t = —")
        ttk.Label(controls, textvariable=self.time_label, width=20).pack(side="left")
        ttk.Label(controls, text="变形放大").pack(side="left")
        self.factor = tk.DoubleVar(value=1)
        entry = ttk.Entry(controls, textvariable=self.factor, width=8)
        entry.pack(side="left", padx=5)
        ttk.Button(controls, text="恢复 1 倍", command=self.reset_scale).pack(side="left", padx=5)
        entry.bind("<Return>", self.draw)
        entry.bind("<FocusOut>", self.draw)
        ttk.Label(panel, text="灰色虚线：原始网格；彩色：所选结果的单元平均值。播放按保存帧顺序显示。").pack(anchor="w")
        self.display_note = tk.StringVar(value="当前按实际位移显示（1 倍）。")
        ttk.Label(panel, textvariable=self.display_note, foreground="#965513").pack(anchor="w")
        self.curves = tk.Canvas(panel, background="white", highlightthickness=0, height=240)
        self.curves.pack(fill="both", expand=True, pady=(10,0))
        self.mesh.bind("<Configure>", self.draw)
        self.curves.bind("<Configure>", self.draw)
        self.summary = tk.StringVar(value="曲线：观察节点位移与全模型能量；仅无外力、无阻尼示例可用总能量漂移衡量误差。")
        ttk.Label(panel, textvariable=self.summary, wraplength=1080).pack(anchor="w", pady=6)
        window.protocol("WM_DELETE_WINDOW", self.close)

    def choose(self):
        path = filedialog.askdirectory(title="选择含 history.csv、frames.csv 的目录")
        if path:
            self.load(path)

    def load(self, path):
        try:
            data = load_results(path)
            self.stop()
            self.data = data
            self.contour=FieldResults.from_folder(path,data[0],data[1],data[3])
            self.selector.configure(values=choices(2) if self.contour else ('位移大小','ux','uy'))
            self.field.set('位移大小');self.contour_controls.fixed.clear()
            self.factor.set(1)
            self.frame.set(0)
            self.timeline.configure(to=len(data[2])-1)
            self.play_button.configure(state="normal" if len(data[2])>1 else "disabled")
            history = data[2]
            metadata = Path(path)/"analysis.txt"
            if metadata.is_file():
                entries = dict(line.split(maxsplit=1) for line in metadata.read_text().splitlines() if line.strip())
                label = {"free":"自由振动", "ramp":"渐增加载", "step":"突加载荷"}.get(entries.get("load"), "动力分析")
                self.description.set(f"{label} · {len(data[0])} 节点 · {len(data[1])} 单元 · dt={float(entries['time_step']):.5g} · 计算时长 {float(entries['duration']):.5g}")
            else:
                self.description.set(f"动力结果 · {len(data[0])} 节点 · {len(data[1])} 单元 · 单位沿用模型")
            e0 = history[0]["total_energy"]
            drift = max(abs(r["total_energy"]-e0) for r in history)
            drift_text = f"{drift/abs(e0):.3%}" if e0 else "初始能量为零，不计算相对变化"
            self.summary.set(f"观察节点 {int(history[0]['node'])} · {len(history)} 帧 · 总能量相对初值的最大变化：{drift_text}（受迫运动中不等于积分误差）")
            self.status.set(f"已载入：{path}")
            if len(history)==1:
                self.description.set(f'Q4 静力云图 · {len(data[0])} 节点 · {len(data[1])} 单元')
                self.summary.set('静力最终状态；第二PK应力 / GL应变；工程剪应变 = 2 E_xy。')
            self.draw()
        except (OSError, ValueError, KeyError) as error:
            messagebox.showerror("无法读取结果", str(error))

    def run(self):
        if self.process:
            return
        self.stop()
        candidates = [ROOT/"build"/"dynamic_demo.exe", ROOT/"build"/"Debug"/"dynamic_demo.exe",
                      ROOT/"build"/"Release"/"dynamic_demo.exe", ROOT/"build"/"dynamic_demo"]
        solver = next((p for p in candidates if p.is_file()), None)
        if not solver:
            messagebox.showerror("未找到示例程序", "请先运行 cmake --build build --target dynamic_demo")
            return
        try:
            self.output = ROOT/"results"/("dynamic_"+datetime.now().strftime("%Y%m%d_%H%M%S_%f"))
            self.log = tempfile.TemporaryFile(mode="w+b")
            self.output.parent.mkdir(parents=True,exist_ok=True)
            self.process = subprocess.Popen([str(solver), self.output.name], cwd=self.output.parent, stdout=self.log, stderr=subprocess.STDOUT,
                creationflags=subprocess.CREATE_NO_WINDOW if sys.platform=="win32" else 0)
            self.run_button.configure(state="disabled")
            self.status.set("正在运行显式自由振动示例…")
            self.window.after(100, self.poll)
        except OSError as error:
            self.log.close()
            messagebox.showerror("无法计算", str(error))

    def poll(self):
        code = self.process.poll()
        if code is None:
            self.window.after(100, self.poll)
            return
        self.process = None
        self.run_button.configure(state="normal")
        self.log.seek(0)
        message = self.log.read().decode("utf-8", errors="replace")
        self.log.close()
        if code:
            self.status.set("计算失败，未载入未完成结果。")
            messagebox.showerror("计算失败", message[-4000:])
        else:
            self.load(self.output)

    def reset_scale(self):
        self.factor.set(1)
        self.draw()

    def export_fields(self):
        if not self.contour:messagebox.showinfo('导出','结果缺少原模型材料信息。',parent=self.window);return
        path=filedialog.asksaveasfilename(parent=self.window,defaultextension='.csv',initialfile='q4_element_fields.csv')
        if path:
            try:
                index=min(len(self.data[2])-1,max(0,round(self.frame.get())))
                self.contour.export(path,index,self.data[2][index]['time'])
            except (ValueError,OSError,OverflowError) as error:messagebox.showerror('导出失败',str(error),parent=self.window)

    def draw(self, *_):
        if not self.data:
            return
        try:
            factor = self.factor.get()
            if not math.isfinite(factor) or not 0 <= factor <= 1e6:
                return
        except tk.TclError:
            return
        nodes, elements, history, frames = self.data
        span = max(max(x for x,y in nodes)-min(x for x,y in nodes), max(y for x,y in nodes)-min(y for x,y in nodes), 1e-12)
        maximum = max(math.hypot(ux,uy) for frame in frames for ux,uy in frame)
        if factor == 1:
            self.display_note.set("当前按实际位移显示（1 倍）。")
        elif factor == 0:
            self.display_note.set("当前倍率为 0，仅显示参考构形；计算位移未改变。")
        elif factor > 1 and maximum*factor > 0.5*span:
            self.display_note.set(f"注意：{factor:g} 倍显示可能造成网格折叠或形状失真，不代表实际单元翻转。请用 1 倍核对。")
        else:
            self.display_note.set(f"当前位移显示倍率 {factor:g}；形状经过缩放，计算结果未改变。")
        index = min(len(history)-1, max(0, round(self.frame.get())))
        self.time_label.set(f"t = {history[index]['time']:.6g}")
        self.mesh.delete("all")
        w, h = max(200,self.mesh.winfo_width()), max(150,self.mesh.winfo_height())
        # A fixed view across all saved frames prevents animation auto-fit jitter.
        xmin=xmax=nodes[0][0]; ymin=ymax=nodes[0][1]
        for frame in frames:
            for (x,y),(ux,uy) in zip(nodes,frame):
                xmin=min(xmin,x,x+factor*ux);xmax=max(xmax,x,x+factor*ux)
                ymin=min(ymin,y,y+factor*uy);ymax=max(ymax,y,y+factor*uy)
        zoom=min((w-90)/max(xmax-xmin,1e-12),(h-60)/max(ymax-ymin,1e-12))
        project=lambda x,y: (45+(x-xmin)*zoom, h-30-(y-ymin)*zoom)
        try:
            if self.contour:values=self.contour.values(index,self.field.get())
            else:
                field=DISPLACEMENT[self.field.get()]
                nodal=[math.hypot(*u) if field=='mag' else u['xy'.index(field[-1])] for u in frames[index]]
                values=[sum(nodal[n] for n in e)/4 for e in elements]
            low,high=self.contour_controls.bounds(self.field.get(),values)
        except (ValueError,OverflowError) as error:self.display_note.set(str(error));return
        centers=[]
        for i,e in enumerate(elements):
            points=[v for n in e for v in project(nodes[n][0]+factor*frames[index][n][0],nodes[n][1]+factor*frames[index][n][1])]
            self.mesh.create_polygon(points,fill=color(values[i],low,high),outline="#236795",width=1.5)
            centers.append((sum(points[::2])/4,sum(points[1::2])/4))
            original=[v for n in (*e,e[0]) for v in project(*nodes[n])]
            self.mesh.create_line(original,fill="#8796a3",dash=(3,3))
        self.contour_controls.draw(self.mesh,values,centers,self.field.get())
        self.mesh.create_text(10,14,text=self.field.get()+' · 单元平均',anchor='w')
        self.curves.delete("all")
        if len(history)==1:
            self.curves.create_text(20,30,text='静力最终状态，无时间曲线。应力为第二PK；GL剪应变采用工程分量。',anchor='w');return
        w,h=max(400,self.curves.winfo_width()),max(220,self.curves.winfo_height())
        self.plot((55,35,w/2-25,h-40),history,[('ux','#26659d'),('uy','#d66b25')],"观察节点位移",history[index]['time'])
        self.plot((w/2+60,35,w-20,h-40),history,[('kinetic_energy','#d66b25'),('strain_energy','#26659d'),('total_energy','#23825e')],"能量 T / U / T+U",history[index]['time'])

    def plot(self, box, history, series, title, current):
        c=self.curves;x0,y0,x1,y1=box
        t0,t1=history[0]['time'],history[-1]['time']
        vals=[r[key] for key,_ in series for r in history]
        lo,hi=min(0,min(vals)),max(vals)
        padding=max((hi-lo)*.08,abs(hi)*1e-6,1e-16);lo-=padding;hi+=padding
        px=lambda t:x0+(t-t0)/max(t1-t0,1e-16)*(x1-x0)
        py=lambda v:y1-(v-lo)/(hi-lo)*(y1-y0)
        c.create_text(x0,y0-20,text=title,anchor="w",fill="#20384e")
        c.create_rectangle(x0,y0,x1,y1,outline="#c5d0d8")
        for key,color in series:
            points=[v for r in history for v in (px(r['time']),py(r[key]))]
            if len(points)>=4:c.create_line(points,fill=color,width=1.6)
        c.create_line(px(current),y0,px(current),y1,fill="#667684",dash=(3,3))
        for value,y in ((hi,y0),(lo,y1)):
            c.create_text(x0-4,y,text=f"{value:.2g}",anchor="e",fill="#526576")
        c.create_text(x0,y1+14,text=f"{t0:.3g}",anchor="w")
        c.create_text(x1,y1+14,text=f"{t1:.3g}  时间",anchor="e")
        labels={"ux":"ux","uy":"uy","kinetic_energy":"动能 T","strain_energy":"应变能 U","total_energy":"总能量"}
        for i,(key,color) in enumerate(series):
            c.create_text(x0+i*95,y1+31,text=labels[key],anchor="w",fill=color)

    def stop(self):
        self.playing=False
        if self.timer is not None:
            self.window.after_cancel(self.timer)
            self.timer=None
        self.play_button.configure(text="播放")

    def toggle(self):
        if self.playing:
            self.stop()
        elif self.data:
            self.playing=True
            self.play_button.configure(text="暂停")
            if round(self.frame.get())>=len(self.data[2])-1:self.frame.set(0)
            self.tick()

    def tick(self):
        self.timer=None
        if not self.playing or not self.data:return
        index=round(self.frame.get())
        if index>=len(self.data[2])-1:
            self.stop();return
        self.frame.set(index+1);self.draw()
        self.timer=self.window.after(35,self.tick)

    def close(self):
        if self.process and self.process.poll() is None:
            if not messagebox.askyesno("正在计算", "停止计算并关闭窗口？"):return
            self.process.terminate();self.process.wait();self.log.close()
        self.stop();self.window.destroy()


if __name__=="__main__":
    root=tk.Tk(); viewer=Viewer(root)
    if len(sys.argv)==2:viewer.load(sys.argv[1])
    root.mainloop()
