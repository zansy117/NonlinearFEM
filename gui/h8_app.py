"""H8 input and results panel in the shared NonlinearFEM workbench."""
from datetime import datetime
from pathlib import Path
import os
import subprocess
import tkinter as tk
from tkinter import ttk,filedialog,messagebox
from h8_model import DEFAULTS,TEMPLATES,geometry,model_text
from viewer3d import Viewer,surface_faces
from runtime import application_root

class H8App:
    def __init__(self,window):
        self.window=window;window.title('NonlinearFEM · 三维 H8 工作台');window.geometry('1200x800+20+20')
        self.process=None;self.output=None;self.log=None
        self.custom=None
        body=ttk.Frame(window,padding=12);body.pack(fill='both',expand=True)
        left=ttk.Frame(body,width=350);left.pack(side='left',fill='y',padx=(0,12))
        right=ttk.Frame(body);right.pack(side='left',fill='both',expand=True)
        ttk.Label(left,text='三维 H8 / 与 Q4 共用求解程序',font=('Microsoft YaHei UI',13,'bold')).pack(anchor='w')
        self.mode=tk.StringVar(value=TEMPLATES['cantilever'])
        box=ttk.Combobox(left,textvariable=self.mode,values=tuple(TEMPLATES.values()),state='readonly',width=40)
        box.pack(fill='x',pady=8);box.bind('<<ComboboxSelected>>',self.preview)
        self.analysis=tk.StringVar(value='静力分析')
        ttk.Combobox(left,textvariable=self.analysis,state='readonly',values=('静力分析','自由振动','渐增加载','突加载荷')).pack(fill='x')
        tabs=ttk.Notebook(left);tabs.pack(fill='x',pady=8)
        model=ttk.Frame(tabs,padding=8);dynamic=ttk.Frame(tabs,padding=8)
        tabs.add(model,text='模型参数');tabs.add(dynamic,text='动力参数')
        self.values={}
        specs=(('length','长度 X',.1,10),('height','高度 Y',.1,5),('depth','厚度 Z',.1,5),
               ('nx','X 方向单元',1,12),('ny','Y 方向单元',1,6),('nz','Z 方向单元',1,6),
               ('young','杨氏模量',100,10000),('poisson','泊松比',-.5,.49),('force','总力（有符号）',-1,1))
        for key,label,low,high in specs:
            row=ttk.Frame(model);row.pack(fill='x',pady=2)
            ttk.Label(row,text=label,width=14).pack(side='left')
            var=tk.StringVar(value=str(DEFAULTS[key]));self.values[key]=var
            ttk.Entry(row,textvariable=var,width=10).pack(side='right')
            def slide(value,k=key):
                v=round(float(value)) if k in ('nx','ny','nz') else float(value)
                self.values[k].set(str(v) if isinstance(v,int) else f'{v:.5g}');self.preview()
            scale=ttk.Scale(model,from_=low,to=high,command=slide);scale.set(DEFAULTS[key]);scale.pack(fill='x')
        self.dynamic={}
        for key,label,value in (('density','参考密度','1000'),('time_step','时间步（0 自动）','0'),
            ('duration','时长','1'),('output_every','每几步保存','10'),('initial_amplitude','初始扰动幅值','0.001')):
            ttk.Label(dynamic,text=label).pack(anchor='w');var=tk.StringVar(value=value);self.dynamic[key]=var
            ttk.Entry(dynamic,textvariable=var).pack(fill='x',pady=(0,8))
        ttk.Label(dynamic,text='无阻尼。自动步长依据初始刚度；大变形时仍需核查稳定性。自由振动忽略外力。',wraplength=300).pack()
        ttk.Label(left,text='一致单位；网格最多 600 节点。编辑数值后点“预览”。',wraplength=330).pack(anchor='w')
        buttons=ttk.Frame(left);buttons.pack(fill='x',pady=8)
        ttk.Button(buttons,text='预览',command=self.preview).pack(side='left')
        self.run_button=ttk.Button(buttons,text='运行计算',command=self.run);self.run_button.pack(side='left',padx=4)
        self.stop_button=ttk.Button(buttons,text='停止',command=self.cancel,state='disabled');self.stop_button.pack(side='left')
        ttk.Button(left,text='导出 H8 模型',command=self.save_model).pack(fill='x',pady=3)
        row=ttk.Frame(left);row.pack(fill='x')
        ttk.Button(row,text='自定义约束与载荷',command=self.edit_boundary).pack(side='left',fill='x',expand=True)
        ttk.Button(row,text='恢复示例边界',command=self.reset_boundary).pack(side='left')
        ttk.Button(left,text='打开已有 H8 结果',command=self.open_results).pack(fill='x',pady=3)
        ttk.Button(left,text='打开输出文件夹',command=self.open_folder).pack(fill='x',pady=3)
        self.status=tk.StringVar(value='准备就绪');ttk.Label(left,textvariable=self.status,wraplength=330).pack(anchor='w',pady=10)
        self.viewer=Viewer(right);self.preview();window.protocol('WM_DELETE_WINDOW',self.close)

    def inputs(self):
        mode=next(k for k,v in TEMPLATES.items() if v==self.mode.get())
        return {k:v.get() for k,v in self.values.items()},mode

    def reset_boundary(self):
        if self.process:return
        self.custom=None;self.preview()

    def edit_boundary(self):
        if self.process:self.status.set('请先停止计算，再编辑边界。');return
        try:
            values,mode=self.inputs();p,nodes,elements,constraints,loads,*_=geometry(values,mode)
            from boundary_data import boundary_records,mesh_key
            if self.custom is not None and self.custom['mesh']!=mesh_key(nodes,elements):
                if not messagebox.askyesno('网格已改变','旧边界不适用于当前网格。以当前示例边界重新编辑？',parent=self.window):return
            else:constraints,loads=boundary_records(nodes,elements,self.custom,constraints,loads)
            from boundary_editor import BoundaryEditor
            def apply(custom):
                self.custom=custom;self.preview()
            BoundaryEditor(self.window,nodes,elements,constraints,loads,apply)
        except ValueError as error:messagebox.showerror('无法编辑',str(error),parent=self.window)

    def preview(self,*_):
        if not hasattr(self,'viewer'):return
        try:
            values,mode=self.inputs();p,nodes,elements,constraints,loads,*_=geometry(values,mode,self.custom)
            self.viewer.stop();self.viewer.data=None;self.viewer.nodes=nodes;self.viewer.elements=elements
            self.viewer.faces=surface_faces(elements);self.viewer.play.configure(state='disabled')
            self.viewer.boundary_constraints=[(n,d) for n,d,v in constraints];self.viewer.boundary_loads=loads
            self.viewer.status.set(f'参考网格 · {len(nodes)} 节点 · {len(elements)} 单元')
            self.viewer.time.set('参考构形');self.viewer.draw();self.status.set(('自定义边界（总力滑块不生效）' if self.custom is not None else '示例边界')+f'：{len(constraints)} 项约束，{len(loads)} 项节点力。蓝点：约束；红箭头：力。')
        except (ValueError,KeyError) as error:self.status.set(str(error))

    def save_model(self):
        try:
            text=model_text(*self.inputs(),self.custom)
            path=filedialog.asksaveasfilename(parent=self.window,defaultextension='.nfem',filetypes=[('H8 model','*.nfem')])
            if path:Path(path).write_text(text,encoding='utf-8')
        except ValueError as error:messagebox.showerror('模型无效',str(error),parent=self.window)

    def run(self):
        if self.process:return
        try:
            values,mode=self.inputs();text=model_text(values,mode,self.custom);observer=geometry(values,mode,self.custom)[-1]
            if self.analysis.get()!='静力分析' and any(v!=0 for n,d,v in geometry(values,mode,self.custom)[3]):
                raise ValueError('动力分析仅支持零位移约束，请修改自定义边界。')
            root=application_root();exe=root/'NonlinearFEM.exe'
            if not exe.is_file():exe=root/'build'/'NonlinearFEM.exe'
            if not exe.is_file():raise ValueError('没有找到求解程序，请先编译 NonlinearFEM。')
            run=root/'results'/('h8_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'));run.mkdir(parents=True)
            model=run/'input.nfem';model.write_text(text,encoding='utf-8');output=run/'solution'
            args=[str(exe),'--h8',model.name,output.name]
            if self.analysis.get()!='静力分析':
                config=run/'dynamic.cfg';load={'自由振动':'free','渐增加载':'ramp','突加载荷':'step'}[self.analysis.get()]
                lines=['NFEM_DYNAMIC 1']+[f'{k} {v.get()}' for k,v in self.dynamic.items()]
                lines += [f'observation_node {observer}',f'load {load}',f'initial_mode {"axial" if mode=="axial" else "shear" if mode=="shear" else "bending"}']
                config.write_text('\n'.join(lines)+'\n',encoding='utf-8')
                args=[str(exe),'--dynamic-h8',model.name,config.name,output.name]
            self.output=output;self.log=(run/'solver.log').open('w',encoding='utf-8')
            self.running_model_text=text
            self.process=subprocess.Popen(args,cwd=run,stdout=self.log,stderr=subprocess.STDOUT,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
            self.run_button.configure(state='disabled');self.stop_button.configure(state='normal');self.status.set('正在计算……可随时停止。')
            self.window.after(150,self.poll)
        except (OSError,ValueError) as error:
            if self.log:self.log.close();self.log=None
            messagebox.showerror('无法运行',str(error),parent=self.window)

    def poll(self):
        if not self.process:return
        code=self.process.poll()
        if code is None:self.window.after(150,self.poll);return
        self.log.close();self.log=None;self.process=None
        self.run_button.configure(state='normal');self.stop_button.configure(state='disabled')
        if code:
            log=(self.output.parent/'solver.log').read_text(encoding='utf-8',errors='replace')
            self.status.set('计算未完成；请查看 solver.log。');messagebox.showerror('求解未完成',log[-3000:],parent=self.window);return
        try:
            if model_text(*self.inputs(),self.custom)!=self.running_model_text:
                self.status.set('旧参数计算完成；当前模型已改变，请重新运行。旧结果可通过“打开已有 H8 结果”查看。');return
            self.viewer.load(self.output);self.status.set(f'完成：{self.output}')
        except (OSError,ValueError,KeyError) as error:messagebox.showerror('结果读取失败',str(error),parent=self.window)

    def open_results(self):
        if self.process:self.status.set('请先等待当前计算完成，或点击停止。');return
        folder=filedialog.askdirectory(parent=self.window,title='选择 H8 solution 文件夹')
        if not folder:return
        try:self.viewer.load(folder);self.output=Path(folder);self.status.set(f'已载入：{folder}')
        except (OSError,ValueError,KeyError) as error:messagebox.showerror('H8 结果无效',str(error),parent=self.window)

    def open_folder(self):
        if self.output:os.startfile(str(self.output if self.output.exists() else self.output.parent))

    def cancel(self):
        if self.process:
            self.process.terminate();self.process.wait();self.log.close();self.log=None;self.process=None
            self.run_button.configure(state='normal');self.stop_button.configure(state='disabled');self.status.set('已停止；未完成结果不作为求解结果。')

    def close(self):
        self.cancel();self.viewer.stop();self.window.destroy()
