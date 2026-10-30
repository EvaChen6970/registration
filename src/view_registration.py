#!/usr/bin/env python3
import argparse
from pathlib import Path
import tifffile,numpy as np
p=argparse.ArgumentParser();p.add_argument('--registered',required=True);p.add_argument('--codex',required=True);p.add_argument('--level',type=int,default=-1);p.add_argument('--x',type=int);p.add_argument('--y',type=int);p.add_argument('--size',type=int,default=2048);a=p.parse_args()
try: import napari
except ImportError: raise SystemExit('Install napari in the registration environment')
def read(path,level):
 with tifffile.TiffFile(path) as t: return t.series[0].levels[level].asarray()
he=read(a.registered,a.level);cx=read(a.codex,-1);cx=cx[0] if cx.ndim==3 else cx
if a.x is not None and a.y is not None:
 s=a.size; he=he[a.y:a.y+s,a.x:a.x+s];cx=cx[a.y:a.y+s,a.x:a.x+s]
v=napari.Viewer();v.add_image(cx,name='CODEX DAPI',colormap='gray');v.add_image(he,name='Registered H&E',rgb=he.ndim==3,opacity=.55);napari.run()
