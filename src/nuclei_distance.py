#!/usr/bin/env python3
import argparse,json
import numpy as np,tifffile
from scipy.spatial import cKDTree
p=argparse.ArgumentParser(description='Nearest HE-to-DAPI nucleus centroid distances after affine registration.')
p.add_argument('--he-labels',required=True);p.add_argument('--dapi-labels',required=True);p.add_argument('--transform',required=True);p.add_argument('--pixel-size-um',type=float);p.add_argument('--out',required=True);a=p.parse_args()
def centroids(path):
 x=tifffile.imread(path); ids,inv,c=np.unique(x,return_inverse=True,return_counts=True); sums=np.column_stack(np.unravel_index(np.bincount(inv,weights=np.indices(x.shape)[0].ravel()),x.shape)) if False else None
 from scipy.ndimage import center_of_mass
 return np.array(center_of_mass(np.ones_like(x),x,ids[ids!=0]))[:,::-1]
he=centroids(a.he_labels);dy=centroids(a.dapi_labels)
j=json.load(open(a.transform));m=np.array(j['full_resolution_affine_2x3']);h=np.c_[he,np.ones(len(he))]@m.T
nearest,idx=cKDTree(dy).query(h);unit='px'
if a.pixel_size_um: nearest*=a.pixel_size_um;unit='um'
np.savetxt(a.out,np.c_[h,dy[idx],nearest],delimiter=',',header=f'he_x_codex,he_y_codex,dapi_x,dapi_y,nearest_distance_{unit}',comments='')
print(f'HE nuclei={len(he)} DAPI nuclei={len(dy)} median={np.median(nearest):.3f} {unit}')
