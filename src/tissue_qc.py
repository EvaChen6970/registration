#!/usr/bin/env python3
import argparse,csv
from pathlib import Path
import numpy as np,tifffile,cv2
p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--out',required=True);p.add_argument('--min-tissue-fraction',type=float,default=.02);a=p.parse_args()
root=Path(a.root); files=[]
for group in ('Normal_CODEX_Images','AML_NSM_CODEX_Images','Normal_H&E_Images','AML_NSM_H&E_Images'): files+=list((root/group).glob('*.qptiff'))
rows=[]
for f in files:
 try:
  with tifffile.TiffFile(f) as t: x=t.series[0].levels[-1].asarray()
  if x.ndim>=3: x=x[0] if x.shape[0]<10 else x[...,:3].mean(-1)
  x=np.asarray(x,dtype=np.float32); nz=x[x>0]; lo,hi=np.percentile(nz,[1,99]) if nz.size else (0,1); u=np.uint8(np.clip((x-lo)/max(hi-lo,1),0,1)*255)
  tissue=(u<220) if 'H&E' in str(f) else (u>20)
  rows.append({'file':str(f),'width':u.shape[1],'height':u.shape[0],'tissue_fraction':float(tissue.mean()),'laplacian_variance':float(cv2.Laplacian(u,cv2.CV_32F).var()),'keep':int(tissue.mean()>=a.min_tissue_fraction)})
 except Exception as e: rows.append({'file':str(f),'error':repr(e),'keep':0})
out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
with out.open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=sorted({k for r in rows for k in r}));w.writeheader();w.writerows(rows)
print(f'wrote {len(rows)} rows to {out}')
