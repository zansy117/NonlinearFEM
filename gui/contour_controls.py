"""Shared numerical color legend, range controls and extrema markers."""
import colorsys
import math
import tkinter as tk
from tkinter import ttk

def color(value,low,high):
    t=max(0,min(1,(value-low)/(high-low))) if high!=low else .5
    return '#%02x%02x%02x'%tuple(round(v*255) for v in colorsys.hsv_to_rgb((1-t)*2/3,.85,.95))

class ContourControls:
    def __init__(self,parent,redraw,export=None):
        row=ttk.Frame(parent);row.pack(fill='x',pady=3)
        self.mode=tk.StringVar(value='自动当前帧');self.low=tk.StringVar(value='0');self.high=tk.StringVar(value='1');self.extrema=tk.BooleanVar(value=True);self.fixed={}
        box=ttk.Combobox(row,textvariable=self.mode,state='readonly',width=13,values=('自动当前帧','固定首次范围','手动范围'));box.pack(side='left');box.bind('<<ComboboxSelected>>',lambda _:redraw())
        for label,var in (('下限',self.low),('上限',self.high)):
            ttk.Label(row,text=label).pack(side='left');entry=ttk.Entry(row,textvariable=var,width=9);entry.pack(side='left');entry.bind('<Return>',lambda _:redraw());entry.bind('<FocusOut>',lambda _:redraw())
        ttk.Checkbutton(row,text='极值',variable=self.extrema,command=redraw).pack(side='left')
        ttk.Button(row,text='重置范围',command=lambda:self.reset(redraw)).pack(side='left')
        if export:ttk.Button(row,text='导出应力应变',command=export).pack(side='left',padx=4)

    def reset(self,redraw):self.fixed.clear();redraw()
    def bounds(self,field,values):
        if self.mode.get()=='手动范围':
            low,high=float(self.low.get()),float(self.high.get())
            if not math.isfinite(low) or not math.isfinite(high) or low>=high:raise ValueError('手动色标必须有限，且下限小于上限')
            return low,high
        bounds=min(values),max(values)
        if self.mode.get()=='固定首次范围':bounds=self.fixed.setdefault(field,bounds)
        return bounds

    def draw(self,canvas,values,centers,field):
        low,high=self.bounds(field,values);w=max(300,canvas.winfo_width())
        for i in range(40):
            fill=color(i,0,39);canvas.create_rectangle(w-188+i*4,12,w-184+i*4,24,fill=fill,outline=fill)
        canvas.create_text(w-188,35,text=f'{low:.5g}',anchor='w',fill='#263f52')
        canvas.create_text(w-28,35,text=f'{high:.5g}',anchor='e',fill='#263f52')
        if self.extrema.get() and values:
            ids=(min(range(len(values)),key=values.__getitem__),max(range(len(values)),key=values.__getitem__))
            for slot,i in enumerate(ids):
                x,y=centers[i];fill='#1e5cac' if slot==0 else '#af3029'
                canvas.create_oval(x-4,y-4,x+4,y+4,outline='white',fill=fill)
                canvas.create_text(x+7,y+(-12 if slot==0 else 12),text=f'{"MIN" if slot==0 else "MAX"} e{i}: {values[i]:.4g}',anchor='w',fill=fill,font=('Segoe UI',9,'bold'))
        return low,high
