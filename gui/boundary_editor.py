"""Node/face selection and DOF editing for the common Q4/H8 workbench."""
import math
import tkinter as tk
from tkinter import ttk,messagebox
from boundary_data import mesh_key,face_nodes,parse_nodes,boundary_records

class BoundaryEditor:
    def __init__(self,parent,nodes,elements,constraints,loads,on_apply):
        self.nodes=nodes;self.elements=elements;self.dim=len(nodes[0]);self.on_apply=on_apply
        self.original=(list(constraints),list(loads));self.constraints={(n,d):v for n,d,v in constraints}
        self.loads={(n,d):v for n,d,v in loads};self.selected=set();self.screen=[]
        self.window=tk.Toplevel(parent);self.window.title('自定义边界 · '+('Q4' if self.dim==2 else 'H8'))
        self.window.geometry('1180x780+25+25');self.window.transient(parent)
        self.window.grab_set()
        body=ttk.Frame(self.window,padding=10);body.pack(fill='both',expand=True)
        left=ttk.Frame(body,width=350);left.pack(side='left',fill='y',padx=(0,10))
        right=ttk.Frame(body);right.pack(side='left',fill='both',expand=True)
        tabs=ttk.Notebook(left);tabs.pack(fill='both',expand=True)
        selection=ttk.Frame(tabs,padding=6);settings=ttk.Frame(tabs,padding=6)
        tabs.add(selection,text='① 选择节点');tabs.add(settings,text='② 边界设置')
        left=selection
        ttk.Label(left,text='1. 选择节点或边界',font=('Microsoft YaHei UI',12,'bold')).pack(anchor='w')
        ttk.Label(left,text='单击切换节点；三维中拖动旋转。重叠节点可旋转后选择，或输入编号。',wraplength=330).pack(anchor='w',pady=5)
        row=ttk.Frame(left);row.pack(fill='x')
        self.selection_mode=tk.StringVar(value='替换选择')
        ttk.Combobox(row,textvariable=self.selection_mode,values=('替换选择','追加选择'),state='readonly',width=12).pack(side='left')
        ttk.Button(row,text='清空选择',command=lambda:self.select(set(),True)).pack(side='left',padx=5)
        faces=ttk.Frame(left);faces.pack(fill='x',pady=5)
        for axis in range(self.dim):
            for side in range(2):
                ttk.Button(faces,text='XYZ'[axis]+(' 最小边/面' if side==0 else ' 最大边/面'),
                    command=lambda a=axis,s=side:self.select(face_nodes(self.nodes,a,s))).grid(row=axis,column=side,sticky='ew')
        self.ids=tk.StringVar();ttk.Entry(left,textvariable=self.ids).pack(fill='x',pady=3)
        ttk.Button(left,text='按编号选择（如 0,2,5-8）',command=self.select_ids).pack(fill='x')
        self.selection_note=tk.StringVar();ttk.Label(left,textvariable=self.selection_note,wraplength=330).pack(anchor='w',pady=7)
        ttk.Button(left,text='下一步：设置边界',command=lambda:tabs.select(settings)).pack(fill='x',pady=8)
        left=settings
        ttk.Label(left,textvariable=self.selection_note,wraplength=330).pack(anchor='w',pady=7)
        ttk.Label(left,text='2. 设置选中节点',font=('Microsoft YaHei UI',12,'bold')).pack(anchor='w')
        row=ttk.Frame(left);row.pack(fill='x',pady=4);self.directions=[]
        for d in range(self.dim):
            var=tk.BooleanVar(value=True);self.directions.append(var)
            ttk.Checkbutton(row,text='XYZ'[d],variable=var).pack(side='left')
        self.displacement=tk.StringVar(value='0');self.force=tk.StringVar(value='0')
        row=ttk.Frame(left);row.pack(fill='x');ttk.Label(row,text='指定位移').pack(side='left');ttk.Entry(row,textvariable=self.displacement,width=13).pack(side='right')
        ttk.Button(left,text='设置勾选方向的位移约束',command=lambda:self.assign('constraints')).pack(fill='x',pady=3)
        row=ttk.Frame(left);row.pack(fill='x');ttk.Label(row,text='每个勾选方向的总节点力').pack(side='left');ttk.Entry(row,textvariable=self.force,width=13).pack(side='right')
        ttk.Button(left,text='设置总力（在选中节点间均分）',command=lambda:self.assign('loads')).pack(fill='x',pady=3)
        ttk.Label(left,text='设置会覆盖选中自由度的旧值。总力按节点均分，并非按边长/面积积分的均布面力。动力约束仅支持零位移。',wraplength=330).pack(anchor='w',pady=3)
        row=ttk.Frame(left);row.pack(fill='x')
        for name,label in (('constraints','移除所选约束'),('loads','移除所选力')):
            ttk.Button(row,text=label,command=lambda n=name:self.remove(n)).pack(side='left')
        self.table=ttk.Treeview(left,columns=('kind','node','direction','value'),show='headings',height=6)
        for key,label,width in (('kind','类型',60),('node','节点',55),('direction','方向',45),('value','值',100)):
            self.table.heading(key,text=label);self.table.column(key,width=width,stretch=True)
        self.table.pack(fill='x',pady=5)
        ttk.Button(left,text='清除全部边界（从空模型开始）',command=self.clear).pack(fill='x',pady=2)
        self.note=tk.StringVar();ttk.Label(left,textvariable=self.note,wraplength=330).pack(anchor='w')
        ttk.Button(left,text='应用到模型并关闭',command=self.apply).pack(fill='x',pady=8)
        if self.dim==3:
            from viewer3d import Viewer
            self.viewer=Viewer(right,nodes,elements);self.viewer.on_pick=self.toggle_node
            self.viewer.show_nodes=True
        else:
            self.viewer=None;self.canvas=tk.Canvas(right,bg='#f5f8fb',highlightthickness=0);self.canvas.pack(fill='both',expand=True)
            self.canvas.bind('<Configure>',lambda _:self.draw());self.canvas.bind('<Button-1>',self.pick)
        self.refresh()

    def select(self,ids,replace=False):
        self.selected=set(ids) if replace or self.selection_mode.get()=='替换选择' else self.selected|set(ids)
        self.refresh()

    def select_ids(self):
        try:self.select(parse_nodes(self.ids.get(),len(self.nodes)))
        except ValueError as error:messagebox.showerror('节点选择',str(error),parent=self.window)

    def toggle_node(self,n):
        if n in self.selected:self.selected.remove(n)
        else:self.selected.add(n)
        self.refresh()

    def pick(self,event):
        if not self.screen:return
        n=min(range(len(self.screen)),key=lambda i:math.hypot(event.x-self.screen[i][0],event.y-self.screen[i][1]))
        if math.hypot(event.x-self.screen[n][0],event.y-self.screen[n][1])<=14:self.toggle_node(n)

    def assign(self,kind):
        try:
            if not self.selected:raise ValueError('请先选择节点')
            directions=[i for i,v in enumerate(self.directions) if v.get()]
            if not directions:raise ValueError('请勾选至少一个方向')
            value=float((self.displacement if kind=='constraints' else self.force).get())
            if not math.isfinite(value):raise ValueError('数值必须有限')
            if kind=='loads':value/=len(self.selected)
            records=getattr(self,kind)
            for n in self.selected:
                for d in directions:records[n,d]=value
            self.refresh()
        except ValueError as error:messagebox.showerror('边界设置',str(error),parent=self.window)

    def remove(self,kind):
        records=getattr(self,kind)
        for key in list(records):
            if key[0] in self.selected and self.directions[key[1]].get():del records[key]
        self.refresh()

    def clear(self):self.constraints.clear();self.loads.clear();self.refresh()

    def refresh(self):
        ids=sorted(self.selected);self.selection_note.set(f'已选 {len(ids)} 节点：'+','.join(map(str,ids[:18]))+('…' if len(ids)>18 else ''))
        self.table.delete(*self.table.get_children())
        for name,label in (('constraints','位移'),('loads','节点力')):
            for (n,d),v in sorted(getattr(self,name).items()):
                if not ids or n in self.selected:self.table.insert('', 'end',values=(label,n,'XYZ'[d],f'{v:.6g}'))
        self.note.set(f'{len(self.constraints)} 项位移约束，{len(self.loads)} 项节点力；未勾选的方向不受本次操作影响。')
        self.draw()

    def draw(self):
        if self.viewer:
            self.viewer.selected_nodes=set(self.selected);self.viewer.boundary_constraints=list(self.constraints)
            self.viewer.boundary_loads=[(n,d,v) for (n,d),v in self.loads.items()];self.viewer.draw();return
        c=self.canvas;c.delete('all');w,h=max(300,c.winfo_width()),max(200,c.winfo_height())
        xmin,ymin=[min(p[d] for p in self.nodes) for d in range(2)];xmax,ymax=[max(p[d] for p in self.nodes) for d in range(2)]
        scale=min((w-80)/max(xmax-xmin,1e-14),(h-80)/max(ymax-ymin,1e-14))
        self.screen=[(40+(x-xmin)*scale,h-40-(y-ymin)*scale) for x,y in self.nodes]
        for e in self.elements:c.create_polygon([v for n in e for v in self.screen[n]],fill='#e7eef5',outline='#607990')
        fixed={n for n,d in self.constraints}
        for n,(x,y) in enumerate(self.screen):
            color='#ec9b20' if n in self.selected else '#2463a4' if n in fixed else '#667785'
            c.create_oval(x-4,y-4,x+4,y+4,fill=color,outline='white')
            if n in self.selected or len(self.nodes)<=80:c.create_text(x+6,y-8,text=str(n),anchor='w')
        for (n,d),v in self.loads.items():
            if v:
                x,y=self.screen[n];delta=25 if v>0 else -25
                c.create_line(x,y,x+delta*(d==0),y-delta*(d==1),arrow='last',fill='#ba4e39',width=2)
        c.create_text(12,15,text='橙色：选中 / 蓝色：有位移约束 / 红箭头：节点力',anchor='w')

    def apply(self):
        custom={'mesh':mesh_key(self.nodes,self.elements),
                'constraints':[(n,d,v) for (n,d),v in sorted(self.constraints.items())],
                'loads':[(n,d,v) for (n,d),v in sorted(self.loads.items())]}
        boundary_records(self.nodes,self.elements,custom,[],[])
        self.on_apply(custom);self.window.destroy()
