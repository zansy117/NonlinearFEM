"""Local slider-based model builder and result viewer. Python standard library only."""
import csv
from datetime import datetime
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import webbrowser

from model import PRESETS, TEMPLATES, DESCRIPTIONS, geometry, model_text

from runtime import application_root, VERSION
from contour_fields import FieldResults,choices
from contour_controls import ContourControls,color
ROOT = application_root()


class App:
    def __init__(self, window):
        self.window = window
        self.window.title(f"NonlinearFEM v{VERSION} · 矩形板建模与计算")
        self.window.geometry("1250x820")
        self.window.minsize(1060, 760)
        self.process = None
        self.custom = None
        self.contour = None
        self.result = None
        self.last_output = None
        self.vars, self.scales = {}, {}
        self.updating = False
        style = ttk.Style()
        style.theme_use("clam")
        style.configure("TFrame", background="#f2f5f8")
        style.configure("TLabel", background="#f2f5f8", font=("Microsoft YaHei UI", 10))
        style.configure("Title.TLabel", font=("Microsoft YaHei UI", 19, "bold"))
        style.configure("TButton", padding=8)
        container = ttk.Frame(window, padding=22)
        container.pack(fill="both", expand=True)
        ttk.Label(container, text="NonlinearFEM / 模型工作台", style="Title.TLabel").pack(anchor="w")
        ttk.Label(container, text="拖动滑块建模 → 运行计算 → 查看位移与变形   ·   请采用一致单位").pack(anchor="w", pady=(5, 15))
        analysis_row = ttk.Frame(container)
        analysis_row.pack(fill="x", pady=(0, 10))
        ttk.Label(analysis_row, text="分析类型  ").pack(side="left")
        self.analysis = tk.StringVar(value="静力分析")
        analysis = ttk.Combobox(analysis_row, textvariable=self.analysis, state="readonly", width=28,
            values=("静力分析", "动力：自由振动", "动力：渐增加载", "动力：突加载荷"))
        analysis.pack(side="left")
        analysis.bind("<<ComboboxSelected>>", self.analysis_changed)
        ttk.Button(analysis_row,text="自定义 Q4 边界",command=self.edit_boundary).pack(side="left",padx=8)
        ttk.Button(analysis_row,text="恢复示例边界",command=self.reset_boundary).pack(side="left")
        def open_h8():
            from h8_app import H8App
            H8App(tk.Toplevel(self.window))
        ttk.Button(analysis_row, text="三维 H8 工作台", command=open_h8).pack(side="right")
        body = ttk.Frame(container)
        body.pack(fill="both", expand=True)
        left = ttk.Frame(body, width=400)
        left.pack(side="left", fill="y", padx=(0, 20))
        right = ttk.Frame(body)
        right.pack(side="left", fill="both", expand=True)
        ttk.Label(left, text="模型与边界条件").pack(anchor="w")
        self.mode = tk.StringVar(value="悬臂板：左端固定，右端竖向力")
        mode = ttk.Combobox(left, textvariable=self.mode, state="readonly", width=38,
            values=tuple(TEMPLATES.values()))
        mode.pack(fill="x", pady=(5, 9))
        mode.bind("<<ComboboxSelected>>", lambda _: self.preset())
        self.add_slider(left, "length", "长度 L", .1, 10)
        self.add_slider(left, "height", "高度 H", .1, 5)
        self.add_slider(left, "nx", "长度方向单元数", 1, 32, integer=True)
        self.add_slider(left, "ny", "高度方向单元数", 1, 16, integer=True)
        self.add_slider(left, "thickness", "参考厚度 t", .01, 2)
        self.add_slider(left, "young", "杨氏模量 E（滑块按对数变化）", 3, 11.5, logarithmic=True)
        self.add_slider(left, "poisson", "泊松比 ν", -.9, .49)
        self.add_slider(left, "force", "加载边界总力（有正负方向）", -200, 200)
        self.boundary_note = tk.StringVar()
        ttk.Label(left, textvariable=self.boundary_note,
                  wraplength=390).pack(anchor="w", pady=(8, 12))
        buttons = ttk.Frame(left)
        buttons.pack(fill="x")
        self.run_button = ttk.Button(buttons, text="运行计算", command=self.run)
        self.run_button.pack(side="left", fill="x", expand=True)
        self.cancel_button = ttk.Button(buttons, text="停止", command=self.cancel, state="disabled")
        self.cancel_button.pack(side="left", padx=6)
        ttk.Button(buttons, text="导出模型", command=self.save_model).pack(side="left")
        ttk.Button(left, text="恢复当前示例参数", command=self.preset).pack(fill="x", pady=8)
        self.progress = ttk.Progressbar(left, mode="indeterminate")
        self.progress.pack(fill="x")
        self.status = tk.StringVar()
        ttk.Label(left, textvariable=self.status, wraplength=390).pack(anchor="w", pady=10)
        self.dynamic_box = ttk.LabelFrame(right, text="动力参数（固定时间步，无阻尼）", padding=8)
        self.dynamic_vars = {}
        for i, (key, label, value) in enumerate((
            ("density", "参考密度", "1000"), ("time_step", "时间步（0=自动）", "0"),
            ("duration", "计算时长", "1"), ("output_every", "每几步保存", "10"),
            ("initial_amplitude", "初始扰动幅值", "0.001"),
            ("observation_node", "观察节点（auto）", "auto"))):
            var = tk.StringVar(value=value)
            self.dynamic_vars[key] = var
            ttk.Label(self.dynamic_box, text=label).grid(row=(i//3)*2, column=i%3, sticky="w", padx=8)
            ttk.Entry(self.dynamic_box, textvariable=var, width=19).grid(row=(i//3)*2+1, column=i%3, sticky="ew", padx=8, pady=3)
            var.trace_add("write", lambda *_: self.changed())
        ttk.Label(self.dynamic_box, text="自由振动忽略外力；加载模式使用左侧总力。初始扰动可设为 0。自动步长只依据初始刚度。",
                  wraplength=650).grid(row=4,column=0,columnspan=3,sticky="w",padx=8,pady=5)
        self.summary = tk.StringVar()
        self.summary_label = ttk.Label(right, textvariable=self.summary)
        self.summary_label.pack(anchor="w", pady=(0, 7))
        self.canvas = tk.Canvas(right, background="white", highlightthickness=0)
        self.canvas.pack(fill="both", expand=True)
        self.canvas.bind("<Configure>", lambda _: self.draw())
        self.canvas.bind("<Motion>", self.hover)
        self.detail = tk.StringVar(value="模型预览：灰色网格 / 蓝色约束 / 橙色载荷箭头")
        ttk.Label(right, textvariable=self.detail, wraplength=690).pack(anchor="w", pady=7)
        controls = ttk.Frame(right)
        controls.pack(fill="x")
        self.field = tk.StringVar(value="位移大小")
        field = ttk.Combobox(controls, textvariable=self.field, state="readonly", width=24,
                             values=choices(2))
        field.pack(side="left")
        field.bind("<<ComboboxSelected>>", lambda _: self.draw())
        ttk.Label(controls, text="  变形放大").pack(side="left")
        self.factor = tk.DoubleVar(value=1)
        ttk.Scale(controls, from_=0, to=10, variable=self.factor, command=lambda _: self.draw()).pack(side="left", fill="x", expand=True)
        self.factor_text = tk.StringVar(value="1.0×")
        ttk.Label(controls, textvariable=self.factor_text, width=6).pack(side="left")
        self.contour_controls=ContourControls(right,self.draw,self.export_fields)
        self.open_button = ttk.Button(right, text="查看完整云图与导出数据", command=self.open_result, state="disabled")
        self.open_button.pack(fill="x", pady=8)
        self.window.protocol("WM_DELETE_WINDOW", self.close)
        self.preset()

    def is_dynamic(self):
        return self.analysis.get() != "静力分析"

    def reset_boundary(self):
        if self.process:return
        self.custom=None;self.boundary_note.set(DESCRIPTIONS[self.key()]);self.changed()

    def edit_boundary(self):
        if self.process:self.status.set('请先停止计算，再编辑边界。');return
        try:
            p,nodes,elements,constraints,loads=geometry(self.values(),self.key())
            from boundary_data import boundary_records,mesh_key
            if self.custom is not None and self.custom['mesh']!=mesh_key(nodes,elements):
                if not messagebox.askyesno('网格已改变','旧边界不适用于当前网格。以当前示例边界重新编辑？',parent=self.window):return
            else:constraints,loads=boundary_records(nodes,elements,self.custom,constraints,loads)
            from boundary_editor import BoundaryEditor
            def apply(custom):
                self.custom=custom;self.changed()
            BoundaryEditor(self.window,nodes,elements,constraints,loads,apply)
        except ValueError as error:messagebox.showerror('无法编辑',str(error),parent=self.window)

    def analysis_changed(self, *_):
        if self.is_dynamic():
            self.dynamic_box.pack(fill="x", pady=(0,8), before=self.summary_label)
        else:
            self.dynamic_box.pack_forget()
        self.open_button.configure(text="查看云图、时程与变形播放" if self.is_dynamic() else "查看完整云图与导出数据")
        self.changed()

    def dynamic_config(self):
        if any(v!=0 for n,d,v in geometry(self.values(),self.key(),self.custom)[3]):
            raise ValueError('动力分析仅支持零位移约束，请修改自定义边界。')
        from dynamic_config import configuration
        return configuration({k:v.get() for k,v in self.dynamic_vars.items()},self.values(),self.key(),self.analysis.get())

    def key(self):
        return next(key for key,label in TEMPLATES.items() if label==self.mode.get())

    def add_slider(self, parent, key, label, lo, hi, integer=False, logarithmic=False):
        row = ttk.Frame(parent)
        row.pack(fill="x", pady=2)
        ttk.Label(row, text=label).grid(row=0, column=0, sticky="w")
        var = tk.StringVar()
        self.vars[key] = var
        entry = ttk.Entry(row, textvariable=var, width=12)
        entry.grid(row=0, column=1, sticky="e")
        entry.bind("<Return>", lambda _: self.sync_entries())
        entry.bind("<FocusOut>", lambda _: self.sync_entries())
        scale = ttk.Scale(row, from_=lo, to=hi,
            command=lambda value, k=key, i=integer, log=logarithmic: self.slide(k, value, i, log))
        scale.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(2, 0))
        row.columnconfigure(0, weight=1)
        self.scales[key] = (scale, integer, logarithmic)
        var.trace_add("write", lambda *_: self.changed())

    def slide(self, key, value, integer, logarithmic):
        if self.updating:
            return
        number = float(value)
        number = 10 ** number if logarithmic else round(number) if integer else number
        self.vars[key].set(str(int(number)) if integer else f"{number:.6g}")

    def values(self):
        try:
            return {key: float(var.get()) for key, var in self.vars.items()}
        except ValueError as error:
            raise ValueError("请在数值框中填写有效数字") from error

    def sync_entries(self):
        self.updating = True
        try:
            for key, (scale, _, logarithmic) in self.scales.items():
                value = float(self.vars[key].get())
                if math.isfinite(value) and (not logarithmic or value > 0):
                    scale.set(math.log10(value) if logarithmic else value)
        except ValueError:
            pass
        finally:
            self.updating = False
        self.changed()

    def preset(self):
        self.custom=None
        self.updating = True
        limit=max(200,abs(PRESETS[self.key()]["force"])*4)
        self.scales["force"][0].configure(from_=-limit,to=limit)
        self.boundary_note.set(DESCRIPTIONS[self.key()]+" 节点力按边长分配。")
        for key, value in PRESETS[self.key()].items():
            self.vars[key].set(str(value))
        self.updating = False
        self.sync_entries()

    def changed(self):
        if self.updating:
            return
        self.result = None
        self.contour = None
        self.contour_controls.fixed.clear()
        self.open_button.configure(state="disabled")
        try:
            p, nodes, elements, constraints, loads = geometry(self.values(), self.key(), self.custom)
            if self.custom is not None:self.boundary_note.set(f'自定义边界：{len(constraints)} 项位移约束、{len(loads)} 项节点力。总力滑块仅作用于示例边界。')
            if self.is_dynamic():
                self.dynamic_config()
            self.summary.set(f"输入模型  /  {p['nx']} × {p['ny']} 单元  ·  {len(nodes)} 节点  ·  {len(elements)} 单元")
            if self.process is None:
                self.status.set("参数已更新。点击运行计算，生成对应的结果。")
                self.run_button.configure(state="normal")
        except ValueError as error:
            self.status.set(str(error))
            self.run_button.configure(state="disabled")
        self.draw()

    def draw(self):
        if not hasattr(self, "canvas"):
            return
        self.canvas.delete("all")
        self.screen = []
        try:
            p, nodes, elements, constraints, loads = geometry(self.values(), self.key(), self.custom)
        except ValueError:
            return
        w, h = max(self.canvas.winfo_width(), 200), max(self.canvas.winfo_height(), 200)
        factor = self.factor.get()
        self.factor_text.set(f"{factor:.1f}×")
        positions = [(x + factor*self.result[i][0], y + factor*self.result[i][1])
                     for i, (x, y) in enumerate(nodes)] if self.result else nodes
        all_points = nodes + positions
        xmin, xmax = min(x for x,y in all_points), max(x for x,y in all_points)
        ymin, ymax = min(y for x,y in all_points), max(y for x,y in all_points)
        zoom = min((w-100)/max(xmax-xmin,1e-12), (h-100)/max(ymax-ymin,1e-12))
        project = lambda x,y: (50+(x-xmin)*zoom, h-50-(y-ymin)*zoom)
        self.screen = [project(x,y) for x,y in positions]
        try:
            if self.result and self.contour is None:self.contour=FieldResults(nodes,elements,[self.result],p['young'],p['poisson'])
            means=self.contour.values(0,self.field.get()) if self.result else []
            low,high=self.contour_controls.bounds(self.field.get(),means) if means else (0,0)
        except (ValueError,OverflowError) as error:self.detail.set(str(error));return
        import colorsys
        for i,e in enumerate(elements):
            points = [v for n in e for v in self.screen[n]]
            fraction = (means[i]-low)/(high-low) if high!=low else .5
            fill=color(means[i],low,high) if self.result else '#e7eef5'
            self.canvas.create_polygon(points, fill=fill, outline="#607990", width=1)
        if self.result:
            for e in elements:
                points = [v for n in (*e,e[0]) for v in project(*nodes[n])]
                self.canvas.create_line(points, fill="#9aa5af", dash=(3,3))
            self.canvas.create_text(18,20,anchor="w",text=f"{self.field.get()} · 单元平均",fill="#24445d")
            centers=[(sum(self.screen[n][0] for n in e)/4,sum(self.screen[n][1] for n in e)/4) for e in elements]
            self.contour_controls.draw(self.canvas,means,centers,self.field.get())
        else:
            fixed = set(n for n,d,v in constraints)
            for n in fixed:
                x,y = self.screen[n]
                self.canvas.create_oval(x-4,y-4,x+4,y+4,fill="#2463a4",outline="")
            if self.is_dynamic() and self.analysis.get()=="动力：自由振动":
                loads = []
            for n,d,value in loads:
                x,y = self.screen[n]
                if value:
                    delta = 25 if value>0 else -25
                    self.canvas.create_line(x,y,x+delta*(d==0),y-delta*(d==1),arrow="last",fill="#d16b21",width=2)
            self.canvas.create_text(18,20,anchor="w",text="x →   y ↑   /   原始构形与边界条件",fill="#24445d")

    def hover(self, event):
        if not self.screen:
            return
        n = min(range(len(self.screen)), key=lambda i: math.hypot(event.x-self.screen[i][0],event.y-self.screen[i][1]))
        if self.result:
            ux,uy = self.result[n]
            self.detail.set(f"节点 {n}  ·  ux = {ux:.6g}  ·  uy = {uy:.6g}  ·  |u| = {math.hypot(ux,uy):.6g}")
        else:
            self.detail.set(f"最近节点 {n}  ·  蓝色点：约束节点  ·  橙色箭头：载荷方向")

    def save_model(self):
        try:
            text = model_text(self.values(), self.key(), self.custom)
            path = filedialog.asksaveasfilename(defaultextension=".nfem", filetypes=[("有限元模型", "*.nfem")])
            if path:
                config = self.dynamic_config() if self.is_dynamic() else None
                Path(path).write_text(text, encoding="utf-8")
                if config is not None:
                    Path(path).with_suffix(".dynamic.cfg").write_text(config,encoding="utf-8")
                self.status.set("模型文件已保存。")
        except (ValueError,OSError) as error:
            messagebox.showerror("无法保存", str(error))

    def export_fields(self):
        if not self.result or self.contour is None:
            messagebox.showinfo('应力应变导出','请先完成计算。',parent=self.window);return
        path=filedialog.asksaveasfilename(parent=self.window,defaultextension='.csv',initialfile='q4_element_fields.csv')
        if path:
            try:
                time=0;frame_id=0
                if getattr(self,'last_dynamic',False) and self.last_output:
                    with (self.last_output/'history.csv').open(encoding='utf-8',newline='') as stream:history=list(csv.DictReader(stream))
                    time=float(history[-1]['time']);frame_id=int(history[-1]['frame'])
                self.contour.export(path,0,time,frame_id)
            except (OSError,ValueError,OverflowError) as error:messagebox.showerror('导出失败',str(error),parent=self.window)

    def run(self):
        if self.process:
            return
        try:
            text = model_text(self.values(), self.key(), self.custom)
            config = self.dynamic_config() if self.is_dynamic() else None
            candidates = [ROOT/"build"/"NonlinearFEM.exe",ROOT/"build"/"Debug"/"NonlinearFEM.exe",
                          ROOT/"build"/"Release"/"NonlinearFEM.exe",ROOT/"build"/"NonlinearFEM"]
            solver = next((p for p in candidates if p.is_file()), None)
            if not solver:
                raise ValueError("尚未找到求解程序，请先按照 README 编译项目。")
            directory = ROOT/"results"/datetime.now().strftime("%Y%m%d_%H%M%S_%f")
            directory.mkdir(parents=True)
            model = directory/"model.nfem"
            model.write_text(text, encoding="utf-8")
            self.active_output = directory/"solution"
            command = [str(solver),model.name,self.active_output.name]
            if config is not None:
                config_path = directory/"dynamic.cfg"
                config_path.write_text(config,encoding="utf-8")
                command = [str(solver),"--dynamic",model.name,config_path.name,self.active_output.name]
            self.running_analysis = self.analysis.get()
            self.running_config = config
            self.log = tempfile.TemporaryFile(mode="w+b")
            self.process = subprocess.Popen(command,cwd=directory,
                stdout=self.log,stderr=subprocess.STDOUT,
                creationflags=subprocess.CREATE_NO_WINDOW if sys.platform=="win32" else 0)
            self.cancelled = False
            self.result = None
            self.run_button.configure(state="disabled")
            self.open_button.configure(state="disabled")
            self.cancel_button.configure(state="normal")
            self.progress.start(15)
            self.running_values, self.running_mode = self.values(), self.key()
            self.running_model_text = text
            self.status.set("正在计算… 可以停止；调整参数后需要重新计算。")
            self.window.after(100,self.poll)
        except (OSError,ValueError) as error:
            messagebox.showerror("无法运行",str(error))

    def poll(self):
        code = self.process.poll()
        if code is None:
            self.window.after(100,self.poll)
            return
        self.process = None
        self.progress.stop()
        self.cancel_button.configure(state="disabled")
        self.log.seek(0)
        message = self.log.read().decode("utf-8",errors="replace")
        self.log.close()
        self.changed()
        if self.cancelled:
            self.status.set("计算已停止，没有显示未完成的结果。")
        elif code:
            self.status.set("计算未成功；可减小载荷或调整模型。")
            messagebox.showerror("计算未成功",message[-4000:])
        else:
            self.last_output = self.active_output
            try:
                same = self.values()==self.running_values and self.key()==self.running_mode and self.analysis.get()==self.running_analysis and model_text(self.values(),self.key(),self.custom)==self.running_model_text
                if self.running_config is not None:
                    same = same and self.dynamic_config()==self.running_config
                    from dynamics import load_results
                    nodes, elements, history, frames = load_results(self.last_output)
                    self.last_dynamic = True
                    if same:
                        self.result = frames[-1]
                        self.summary.set(f"动力计算完成 / {len(history)} 帧 · t={history[-1]['time']:.6g} · 此处显示最终状态")
                        self.draw()
                    self.status.set("动力计算完成，点击下方按钮查看时程与播放。" if same else "旧参数动力计算完成；当前输入已改变，请重新运行。")
                    self.open_button.configure(state="normal", text="查看上次动力时程与变形播放")
                    if same:
                        self.open_result()
                    return
                self.last_dynamic = False
                if same:
                    with (self.last_output/"nodes.csv").open(newline="",encoding="utf-8") as stream:
                        rows = list(csv.DictReader(stream))
                    self.result = [(float(row["ux"]),float(row["uy"])) for row in rows]
                    maximum = max(math.hypot(*u) for u in self.result)
                    self.summary.set(f"计算完成  /  最大位移 {maximum:.6g}  ·  原始网格为灰色虚线")
                    self.draw()
                    self.status.set(f"计算成功。结果保存在：\n{self.last_output.relative_to(ROOT)}")
                else:
                    self.status.set("上一组参数已计算完成；当前输入已改变，请重新运行。")
                self.open_button.configure(state="normal")
            except (OSError,ValueError,KeyError) as error:
                self.status.set(f"读取结果失败：{error}")

    def cancel(self):
        if self.process and self.process.poll() is None:
            self.cancelled = True
            self.process.terminate()

    def open_result(self):
        if self.last_output and getattr(self,"last_dynamic",False):
            from dynamics import Viewer
            child = tk.Toplevel(self.window)
            viewer = Viewer(child)
            viewer.load(self.last_output)
        elif self.last_output:
            from dynamics import Viewer
            child=tk.Toplevel(self.window);viewer=Viewer(child);viewer.load(self.last_output)

    def close(self):
        if self.process and self.process.poll() is None:
            if not messagebox.askyesno("计算进行中", "关闭窗口并停止当前计算？"):
                return
            self.process.terminate()
            self.process.wait()
            self.log.close()
        self.window.destroy()


if __name__ == "__main__":
    if len(sys.argv)>1 and sys.argv[1]=='--standard-cases':
        from run_standard_cases import main
        sys.argv=[sys.argv[0]]+sys.argv[2:]
        sys.exit(main(default_project=ROOT))
    if len(sys.argv)==3 and sys.argv[1]=="--self-test":
        from package_smoke import check
        sys.exit(check(sys.argv[2]))
    root = tk.Tk()
    App(root)
    root.mainloop()
