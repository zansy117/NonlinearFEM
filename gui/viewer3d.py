"""Orthographic, depth-sorted H8 surface viewer. No external Python packages."""
import csv
import colorsys
import math
from pathlib import Path
import tkinter as tk
from tkinter import ttk,messagebox,filedialog
from contour_fields import FieldResults,choices,DISPLACEMENT
from contour_controls import ContourControls,color as contour_color

# Local faces for the author's H8 node order; orientation is irrelevant to this viewer.
FACES=((0,1,2,3),(5,4,7,6),(0,3,4,5),(1,6,7,2),(0,5,6,1),(3,2,7,4))

def surface_faces(elements):
    faces={}
    for e in elements:
        for local in FACES:
            face=tuple(e[i] for i in local);key=tuple(sorted(face))
            if key in faces:faces[key]=None
            else:faces[key]=face
    return [f for f in faces.values() if f is not None]

def load_results(folder):
    folder=Path(folder)
    def read(name):
        with (folder/name).open(encoding='utf-8',newline='') as stream:return list(csv.DictReader(stream))
    rows=read('nodes.csv');nodes=[tuple(float(r[k]) for k in ('x','y','z')) for r in rows]
    if [int(r['node']) for r in rows]!=list(range(len(nodes))):raise ValueError('节点编号不连续')
    elements=[tuple(int(r['n'+str(i)]) for i in range(8)) for r in read('elements.csv')]
    history=[{k:float(v) for k,v in row.items()} for row in read('history.csv')]
    if not nodes or not elements or not history:raise ValueError('三维结果为空')
    if any(n<0 or n>=len(nodes) for e in elements for n in e):raise ValueError('单元节点越界')
    frames=[[None]*len(nodes) for _ in history]
    for row in read('frames.csv'):
        f,n=int(row['frame']),int(row['node'])
        if not 0<=f<len(frames) or not 0<=n<len(nodes) or frames[f][n] is not None:raise ValueError('帧或节点记录无效')
        if abs(float(row['time'])-history[f]['time'])>1e-12:raise ValueError('帧时间不一致')
        frames[f][n]=tuple(float(row[k]) for k in ('ux','uy','uz'))
    if any(u is None for frame in frames for u in frame):raise ValueError('帧记录不完整')
    numbers=[v for n in nodes for v in n]+[v for frame in frames for u in frame for v in u]+[v for row in history for v in row.values()]
    if not all(math.isfinite(v) for v in numbers):raise ValueError('结果存在非有限数值')
    if any(history[i]['time']<=history[i-1]['time'] for i in range(1,len(history))):raise ValueError('时刻未递增')
    metadata={}
    if (folder/'analysis.txt').is_file():
        metadata=dict(line.split(maxsplit=1) for line in (folder/'analysis.txt').read_text().splitlines() if line.strip())
    return nodes,elements,history,frames,metadata

class Viewer:
    def __init__(self,window,nodes=None,elements=None,folder=None):
        self.window=window
        if hasattr(window,'title'):window.title('NonlinearFEM · H8 三维变形与时程')
        self.data=None;self.timer=None;self.playing=False;self.angle=-.6;self.elevation=.45;self.zoom=1.;self.drag_point=None
        self.on_pick=None;self.show_nodes=False;self.selected_nodes=set();self.boundary_constraints=[];self.boundary_loads=[];self.press_point=None;self.screen_nodes=[]
        self.contour=None
        panel=ttk.Frame(window,padding=12);panel.pack(fill='both',expand=True)
        ttk.Label(panel,text='H8 三维结果 / 拖动旋转 · 滚轮缩放',font=('Microsoft YaHei UI',15,'bold')).pack(anchor='w')
        self.status=tk.StringVar(value='参考网格预览');ttk.Label(panel,textvariable=self.status).pack(anchor='w',pady=5)
        controls=ttk.Frame(panel);controls.pack(fill='x')
        self.field=tk.StringVar(value='位移大小')
        self.selector=ttk.Combobox(controls,textvariable=self.field,state='readonly',width=24,values=('位移大小','ux','uy','uz'));self.selector.pack(side='left');self.selector.bind('<<ComboboxSelected>>',self.draw)
        ttk.Label(controls,text='  变形倍率 ').pack(side='left');self.factor=tk.StringVar(value='1')
        entry=ttk.Entry(controls,textvariable=self.factor,width=8);entry.pack(side='left');entry.bind('<Return>',self.draw);entry.bind('<FocusOut>',self.draw)
        ttk.Button(controls,text='恢复视图',command=self.reset).pack(side='left',padx=6)
        self.show_original=tk.BooleanVar(value=True);ttk.Checkbutton(controls,text='原始网格',variable=self.show_original,command=self.draw).pack(side='left')
        self.contour_controls=ContourControls(panel,self.draw,self.export_fields)
        self.canvas=tk.Canvas(panel,bg='#f5f8fb',height=340,highlightthickness=0);self.canvas.pack(fill='both',expand=True,pady=8)
        self.canvas.bind('<Configure>',self.draw);self.canvas.bind('<ButtonPress-1>',self.begin_drag);self.canvas.bind('<B1-Motion>',self.rotate);self.canvas.bind('<MouseWheel>',self.wheel)
        self.canvas.bind('<ButtonRelease-1>',self.release)
        bar=ttk.Frame(panel);bar.pack(fill='x');self.play=ttk.Button(bar,text='播放',command=self.toggle);self.play.pack(side='left')
        self.index=tk.DoubleVar(value=0);self.timeline=ttk.Scale(bar,from_=0,to=1,variable=self.index,command=self.draw);self.timeline.pack(side='left',fill='x',expand=True,padx=10)
        self.time=tk.StringVar(value='参考构形');ttk.Label(bar,textvariable=self.time,width=20).pack(side='left')
        self.note=tk.StringVar(value='表面节点位移平均值着色，不是应力云图。');ttk.Label(panel,textvariable=self.note,wraplength=850).pack(anchor='w',pady=5)
        self.curves=tk.Canvas(panel,bg='white',height=185,highlightthickness=0);self.curves.pack(fill='x');self.curves.bind('<Configure>',self.draw)
        self.nodes=nodes or [];self.elements=elements or [];self.faces=surface_faces(self.elements)
        window.bind('<Destroy>',self.destroy,add='+')
        if folder:self.load(folder)
        else:self.play.configure(state='disabled');self.draw()

    def load(self,folder):
        self.stop()
        self.data=load_results(folder);self.nodes,self.elements,self.history,self.frames,self.metadata=self.data
        self.contour=FieldResults.from_folder(folder,self.nodes,self.elements,self.frames)
        self.selector.configure(values=choices(3) if self.contour else ('位移大小','ux','uy','uz'))
        self.field.set('位移大小');self.contour_controls.fixed.clear()
        self.face_elements={tuple(sorted(face)):i for i,e in enumerate(self.elements) for local in FACES for face in [tuple(e[j] for j in local)]}
        self.max_displacement=max(math.hypot(*u) for frame in self.frames for u in frame)
        self.limits={'位移大小':(0,self.max_displacement)}
        for d,key in enumerate(('ux','uy','uz')):
            self.limits[key]=(min(u[d] for frame in self.frames for u in frame),max(u[d] for frame in self.frames for u in frame))
        self.faces=surface_faces(self.elements);self.index.set(0);self.factor.set('1');self.timeline.configure(to=max(1,len(self.frames)-1))
        self.play.configure(state='normal' if len(self.frames)>1 else 'disabled')
        kind='静力最终状态' if self.metadata.get('analysis')=='static' else '动力时程'
        self.status.set(f'{kind} · {len(self.nodes)} 节点 · {len(self.elements)} 单元 · {len(self.frames)} 帧')
        self.draw()

    def reset(self):
        self.angle=-.6;self.elevation=.45;self.zoom=1.;self.factor.set('1');self.draw()
    def export_fields(self):
        if not self.contour or not self.data:messagebox.showinfo('导出','需要载入包含原模型材料信息的计算结果。',parent=self.window);return
        path=filedialog.asksaveasfilename(parent=self.window,defaultextension='.csv',initialfile='h8_element_fields.csv')
        if path:
            try:
                index=min(len(self.frames)-1,max(0,round(self.index.get())))
                self.contour.export(path,index,self.history[index]['time'])
            except (ValueError,OSError,OverflowError) as error:messagebox.showerror('导出失败',str(error),parent=self.window)
    def begin_drag(self,event):self.drag_point=(event.x,event.y);self.press_point=self.drag_point
    def release(self,event):
        if self.on_pick and self.press_point and math.hypot(event.x-self.press_point[0],event.y-self.press_point[1])<4:
            candidates=[i for i,p in enumerate(self.screen_nodes) if math.hypot(event.x-p[0],event.y-p[1])<14]
            if candidates:
                nearest=min(math.hypot(event.x-self.screen_nodes[i][0],event.y-self.screen_nodes[i][1]) for i in candidates)
                candidates=[i for i in candidates if math.hypot(event.x-self.screen_nodes[i][0],event.y-self.screen_nodes[i][1])<=nearest+2]
                self.on_pick(max(candidates,key=lambda i:self.screen_nodes[i][2]))
        self.drag_point=None;self.press_point=None
    def rotate(self,event):
        if self.drag_point:
            x,y=self.drag_point;self.angle+=(event.x-x)*.008;self.elevation+=(event.y-y)*.008
            self.drag_point=(event.x,event.y);self.draw()
    def wheel(self,event):self.zoom=max(.15,min(8,self.zoom*(1.1 if event.delta>0 else 1/1.1)));self.draw()

    def draw(self,*_):
        if not self.nodes:return
        try:factor=float(self.factor.get())
        except ValueError:self.note.set('请输入有限的非负变形倍率，按回车应用。');return
        if not math.isfinite(factor) or not 0<=factor<=1e6:self.note.set('倍率应在 0～1000000 之间。');return
        self.canvas.delete('all')
        index=min(len(self.frames)-1,max(0,round(self.index.get()))) if self.data else 0
        u=self.frames[index] if self.data else [(0,0,0)]*len(self.nodes)
        # Use a fixed spatial bound across saved states; never auto-fit each frame.
        extent=max(max(n[d] for n in self.nodes)-min(n[d] for n in self.nodes) for d in range(3))
        center=[(max(n[d] for n in self.nodes)+min(n[d] for n in self.nodes))/2 for d in range(3)]
        max_u=self.max_displacement if self.data else 0
        w,h=max(self.canvas.winfo_width(),300),max(self.canvas.winfo_height(),200)
        scale=.6*min(w,h)*self.zoom/max(extent+2*factor*max_u,1e-12)
        ca,sa=math.cos(self.angle),math.sin(self.angle);ce,se=math.cos(self.elevation),math.sin(self.elevation)
        def project(point):
            x,y,z=[point[d]-center[d] for d in range(3)]
            xx=ca*x+sa*z;zz=-sa*x+ca*z;yy=ce*y-se*zz;depth=se*y+ce*zz
            return w/2+scale*xx,h/2-scale*yy,depth
        deformed=[project(tuple(p[d]+factor*a[d] for d in range(3))) for p,a in zip(self.nodes,u)]
        original=[project(p) for p in self.nodes]
        self.screen_nodes=original
        try:
            if self.data and self.contour:
                cell_values=self.contour.values(index,self.field.get())
                face_values=[cell_values[self.face_elements[tuple(sorted(face))]] for face in self.faces]
            else:
                field=self.field.get() if self.field.get() in DISPLACEMENT else '位移大小'
                values=[a['xyz'.index(DISPLACEMENT[field][-1])] if field!='位移大小' else math.hypot(*a) for a in u]
                face_values=[sum(values[n] for n in face)/4 for face in self.faces]
                cell_values=[sum(values[n] for n in e)/8 for e in self.elements]
            lo,hi=self.contour_controls.bounds(self.field.get(),cell_values) if self.data else (0,0)
        except (ValueError,OverflowError) as error:self.note.set(str(error));return
        for i in sorted(range(len(self.faces)),key=lambda i:sum(deformed[n][2] for n in self.faces[i])/4):
            face=self.faces[i];t=(face_values[i]-lo)/(hi-lo) if hi!=lo else .5
            fill=contour_color(face_values[i],lo,hi)
            self.canvas.create_polygon([v for n in face for v in deformed[n][:2]],fill=fill if self.data else '#d5e5f1',outline='#395a73')
        if self.show_original.get():
            edges=set()
            for face in self.faces:
                for a,b in zip(face,face[1:]+face[:1]):edges.add(tuple(sorted((a,b))))
            for a,b in edges:self.canvas.create_line(*original[a][:2],*original[b][:2],fill='#8997a1',dash=(3,3))
        if not self.data:
            fixed={n for n,d in self.boundary_constraints}
            for n in sorted(range(len(original)),key=lambda n:original[n][2]):
                if not self.show_nodes and n not in fixed and n not in self.selected_nodes:continue
                x,y,_=original[n];color='#ec9b20' if n in self.selected_nodes else '#2463a4' if n in fixed else '#667785'
                self.canvas.create_oval(x-4,y-4,x+4,y+4,fill=color,outline='white')
                if n in self.selected_nodes:self.canvas.create_text(x+6,y-8,text=str(n),anchor='w')
            for n,d,value in self.boundary_loads:
                if value:
                    point=list(self.nodes[n]);point[d]+=extent*.15*(1 if value>0 else -1)
                    self.canvas.create_line(*original[n][:2],*project(point)[:2],arrow='last',width=2,fill='#ba4e39')
        for label,point,color in [('X',(extent*.25,0,0),'#bc4b3e'),('Y',(0,extent*.25,0),'#37965a'),('Z',(0,0,extent*.25),'#3b6db2')]:
            origin=project(center);tip=project(tuple(center[d]+point[d] for d in range(3)))
            dx,dy=tip[0]-origin[0],tip[1]-origin[1];self.canvas.create_line(45,h-40,45+dx,h-40+dy,fill=color,arrow='last',width=2);self.canvas.create_text(45+dx,h-40+dy,text=label,fill=color,anchor='s')
        warning=' 放大可能导致显示折叠，请用 1 倍核对。' if factor>1 and factor*max_u>.5*extent else ''
        self.note.set(f'{self.field.get()} · 单元平均 · 色标 {lo:.5g} ～ {hi:.5g} · 倍率 {factor:g}。第二PK应力；GL剪应变为工程分量。'+warning)
        if self.data:
            centers=[(sum(deformed[n][0] for n in e)/8,sum(deformed[n][1] for n in e)/8) for e in self.elements]
            self.contour_controls.draw(self.canvas,cell_values,centers,self.field.get())
        if self.data:
            self.time.set('静力最终状态' if self.metadata.get('analysis')=='static' else f"t={self.history[index]['time']:.6g}")
        self.draw_curves(index)

    def draw_curves(self,index):
        c=self.curves;c.delete('all')
        if not self.data or len(self.history)<2:
            c.create_text(15,25,anchor='w',text='静力结果：只有最终位移；反力见 reactions.csv，加载记录见 load_history.csv。');return
        w,h=max(c.winfo_width(),500),185
        for plot,series in enumerate(((('ux','#ba4e39'),('uy','#238257'),('uz','#306baa')),(('kinetic_energy','#ba4e39'),('strain_energy','#306baa'),('total_energy','#238257')))):
            x0=55+plot*w/2;x1=(plot+1)*w/2-20;y0=30;y1=145
            vals=[r[key] for r in self.history for key,_ in series];lo=min(0,min(vals));hi=max(vals);pad=max((hi-lo)*.08,1e-16);lo-=pad;hi+=pad
            t0=self.history[0]['time'];t1=self.history[-1]['time'];px=lambda t:x0+(t-t0)/(t1-t0)*(x1-x0);py=lambda v:y1-(v-lo)/(hi-lo)*(y1-y0)
            c.create_text(x0,13,anchor='w',text='观察节点位移' if plot==0 else '动能 / 应变能 / 总能量（受迫时不恒定）')
            c.create_rectangle(x0,y0,x1,y1,outline='#ccd6dd')
            for key,color in series:
                c.create_line([v for row in self.history for v in (px(row['time']),py(row[key]))],fill=color,width=1.5)
            c.create_line(px(self.history[index]['time']),y0,px(self.history[index]['time']),y1,dash=(3,3),fill='#777')
            c.create_text(x0-3,y0,text=f'{hi:.2g}',anchor='e');c.create_text(x0-3,y1,text=f'{lo:.2g}',anchor='e')
            for i,(key,color) in enumerate(series):c.create_text(x0+i*95,165,anchor='w',text=key,fill=color)
    def toggle(self):
        if self.playing:self.stop();return
        if not self.data or len(self.frames)<2:return
        self.playing=True;self.play.configure(text='暂停')
        if round(self.index.get())>=len(self.frames)-1:self.index.set(0)
        self.tick()
    def tick(self):
        self.timer=None
        if not self.playing:return
        i=round(self.index.get())
        if i>=len(self.frames)-1:self.stop();return
        self.index.set(i+1);self.draw();self.timer=self.window.after(40,self.tick)
    def stop(self):
        self.playing=False
        if self.timer:self.window.after_cancel(self.timer);self.timer=None
        if self.play.winfo_exists():self.play.configure(text='播放')
    def destroy(self,event):
        if event.widget==self.window:self.stop()
