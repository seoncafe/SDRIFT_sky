#!/usr/bin/env python3
"""Python desktop controls for the ideal-orbit airglow viewer."""
from datetime import datetime, timezone, timedelta
import queue
import threading
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from airglow_orbit import AirglowOrbitModel, ORBIT_TYPES, BANDS, DEFAULT_DATA, utc_time


class Viewer:
    def __init__(self, root, data_dir=DEFAULT_DATA, kind='total'):
        self.root=root; self.sky=AirglowOrbitModel(data_dir); self.result=None; self.busy=False
        self.messages=queue.Queue(); self.vars={}
        root.title('S-DRIFT | orbit sky background'); root.geometry('1380x1000')
        controls=ttk.Frame(root,padding=10); controls.pack(fill='x')
        def row():
            f=ttk.Frame(controls); f.pack(fill='x',pady=4); return f
        def entry(frame,key,label,value,values=None,width=12):
            ttk.Label(frame,text=label).pack(side='left',padx=(0,5))
            var=tk.StringVar(value=value); self.vars[key]=var
            widget=(ttk.Entry(frame,textvariable=var,width=width) if values is None else
                    ttk.Combobox(frame,textvariable=var,values=values,state='readonly',width=width))
            widget.pack(side='left',padx=(0,12)); return widget
        self.utc=datetime.now(timezone.utc).replace(microsecond=0)
        f=row(); entry(f,'date','UTC',self.utc.isoformat().replace('+00:00','Z'),width=23)
        ttk.Button(f,text='Now',command=self.now).pack(side='left',padx=4)
        for text,sec in [('-1 h',-3600),('-1 min',-60),('+1 min',60),('+1 h',3600)]:
            ttk.Button(f,text=text,command=lambda n=sec:self.shift(n)).pack(side='left',padx=3)
        f=row(); entry(f,'kind','Map','Zodiacal + DGL + Airglow' if kind=='total' else 'Airglow',('Airglow','Zodiacal + DGL + Airglow'),28)
        f=row(); entry(f,'orbit','Orbit','dawn-dusk',ORBIT_TYPES,16)
        entry(f,'height','H (km)','650',width=7); entry(f,'orbital_phase_deg','북쪽 기준 궤도 위상 φ (°)','0',width=7)
        entry(f,'band','Band','VIS',BANDS,6); entry(f,'model','Model','combined',('combined','glow'),10)
        entry(f,'nside','NSIDE','128',('16','32','64','128'),5)
        f=row()
        entry(f,'min_elv_deg','Min ELV (deg)','20.0',width=6)
        self.evolve=tk.BooleanVar(value=False)
        ttk.Checkbutton(f,text='Evolve phase from epoch',variable=self.evolve).pack(side='left',padx=5)
        entry(f,'epoch','Phase epoch UTC','2026-10-17T00:00:00Z',width=23)
        f=row(); self.sun=tk.BooleanVar(value=True); self.moon=tk.BooleanVar(value=True)
        ttk.Checkbutton(f,text='Sun avoidance',variable=self.sun).pack(side='left')
        entry(f,'sun_separation','Min separation (deg)','90',width=6)
        ttk.Checkbutton(f,text='Moon avoidance',variable=self.moon).pack(side='left')
        entry(f,'moon_separation','Min separation (deg)','60',width=6)
        combo=entry(f,'contour','Contour (mag)','0.5',tuple(f'{i/10:.1f}' for i in range(1,11)),5)
        combo.bind('<<ComboboxSelected>>',lambda _:self.draw())
        f=row(); self.button=ttk.Button(f,text='Calculate map',command=self.update);self.button.pack(side='left')
        for kind in ['pdf','png','fits.gz']:
            ttk.Button(f,text='Save '+kind.upper(),command=lambda k=kind:self.save(k)).pack(side='left',padx=5)
        ttk.Label(controls,text='Ideal solar-following circles. Both: φ=0 projected north, 180 south. NM: φ=90 midnight, 270 noon. Zenith: 0 up, 90 horizontal, 180 down.').pack(anchor='w')
        ttk.Label(controls,text='Mean GLOW/PALACE partial bands; no UV absorption or twilight history. Moon is geocentric (no spacecraft parallax).').pack(anchor='w')
        self.status=tk.StringVar(value='Preparing map'); ttk.Label(root,textvariable=self.status,padding=8).pack(fill='x')
        self.figure=Figure(figsize=(13,7.5)); self.canvas=FigureCanvasTkAgg(self.figure,master=root)
        self.canvas.get_tk_widget().pack(fill='both',expand=True)
        NavigationToolbar2Tk(self.canvas,root).update()
        self.canvas.mpl_connect('motion_notify_event',self.cursor)
        root.after(100,self.update); root.after(100,self.poll)

    def now(self):
        self.vars['date'].set(datetime.now(timezone.utc).isoformat(timespec='seconds').replace('+00:00','Z'))
        self.update()

    def shift(self,seconds):
        try:
            stamp=datetime.fromtimestamp(utc_time(self.vars['date'].get())[0],timezone.utc)+timedelta(seconds=seconds)
            self.vars['date'].set(stamp.isoformat().replace('+00:00','Z')); self.update()
        except ValueError as exc: messagebox.showerror('UTC',str(exc))

    def update(self):
        if self.busy:return
        try:
            v={k:x.get() for k,x in self.vars.items()}; settings={k:float(v[k]) for k in ['height','orbital_phase_deg','min_elv_deg','sun_separation','moon_separation']}
            settings.update(kind='total' if v['kind']=='Zodiacal + DGL + Airglow' else 'airglow',orbit=v['orbit'],band=v['band'],model=v['model'],sun_avoidance=self.sun.get(),moon_avoidance=self.moon.get(),epoch=v['epoch'] if self.evolve.get() else None)
            utc_time(v['date']); nside=int(v['nside'])
        except ValueError as exc: messagebox.showerror('Input',str(exc)); return
        self.busy=True; self.button.configure(state='disabled'); self.status.set('Calculating… (first Numba compilation may take several seconds)')
        def work():
            try:self.messages.put(('ok',self.sky.sky_map(v['date'],nside=nside,**settings)))
            except Exception as exc:self.messages.put(('error',str(exc)))
        threading.Thread(target=work,daemon=True).start()

    def poll(self):
        try:
            state,value=self.messages.get_nowait();self.busy=False;self.button.configure(state='normal')
            if state=='error':self.status.set('Calculation failed');messagebox.showerror('Airglow',value)
            else:
                self.result=value;self.draw();self.status.set(f"{value['engine']} | φ={value['state']['orbital_phase_deg']:.2f} deg | SZA={value['state']['sza']:.2f} deg | {int(value['view'].sum()):,} allowed pixels")
        except queue.Empty:pass
        except Exception as exc:messagebox.showerror('Plot',str(exc))
        self.root.after(100,self.poll)

    def draw(self):
        if self.result is None:return
        self.sky.plot(self.result,contour_step=float(self.vars['contour'].get()),figure=self.figure);self.canvas.draw()

    def cursor(self,event):
        if self.result is None or event.inaxes is not self.figure.axes[0] or event.xdata is None:return
        import numpy as np
        import healpy as hp
        ra=(-np.rad2deg(event.xdata))%360; dec=np.rad2deg(event.ydata)
        i=hp.ang2pix(self.result['nside'],ra,dec,lonlat=True);r=self.result
        self.status.set(f"{r['mag'][i]:.3f} mag at (RA, DEC) = ({ra:.2f}, {dec:.2f}) | zenith {r['theta'][i]:.2f} deg | ELV {r['elv'][i]:.2f} deg"+(' | excluded' if not r['view'][i] else ' | Earth' if r['earth'][i] else ' | partial model' if not r['valid'][i] else ''))

    def save(self,kind):
        if self.result is None:return
        r=self.result;stamp=r['utc'].replace('-','').replace(':','').replace('T','_').rstrip('Z')
        name=f"{r['kind']}_{r['band']}_nside{r['nside']}_{stamp}_{r['orbit']}_H{r['height']:g}_phi{r['state']['orbital_phase_deg']:.3f}_elv{r['min_elv_deg']:g}_{r['model']}"+ (f"_moon{r['moon_separation']:g}" if r['moon_avoidance'] else '')+('_anti-sun' if r['sun_avoidance'] else '')
        path=filedialog.asksaveasfilename(defaultextension='.'+kind,initialfile=name+'.'+kind)
        if path:
            if kind=='fits.gz':self.sky.save_fits(r,path)
            else:self.figure.savefig(path,dpi=160)
            self.status.set('Saved '+path)


def main(data_dir=DEFAULT_DATA, kind='total'):
    root=tk.Tk();Viewer(root,data_dir,kind=kind);root.mainloop()

if __name__=='__main__':main()
