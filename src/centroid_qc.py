#!/usr/bin/env python3
"""Extract DAPI/HE nuclear centroids from registered images and score distances."""
import argparse, csv, json
from pathlib import Path
import cv2
import numpy as np
import tifffile
from scipy.spatial import cKDTree


def read_level(path, level, sp_channel=0, he=False):
    with tifffile.TiffFile(path) as t:
        s=t.series[0]
        if level < 0 or level >= len(s.levels):
            raise ValueError(f'{path}: requested level {level}, available 0..{len(s.levels)-1}')
        arr=s.levels[level].asarray()
        axes=s.axes
        shape=tuple(s.shape)
        y_axis, x_axis = axes.index('Y'), axes.index('X')
        full=(int(shape[y_axis]),int(shape[x_axis]))
    if he:
        if arr.ndim == 2: rgb=cv2.cvtColor(arr.astype(np.uint8), cv2.COLOR_GRAY2RGB)
        else: rgb=np.asarray(arr[..., :3])
        return rgb, full
    if arr.ndim == 2: dapi=arr
    else:
        channel_axis=axes.index('C') if 'C' in axes else None
        if channel_axis is None: raise ValueError(f'Cannot find channel axis in SP axes={axes}, shape={arr.shape}')
        dapi=np.take(arr,sp_channel,axis=channel_axis)
        remaining=axes.replace('C','')
        while dapi.ndim > 2:
            dapi=np.take(dapi,0,axis=0)
            remaining=remaining[1:]
        if remaining != 'YX':
            y_axis,x_axis=remaining.index('Y'),remaining.index('X')
            dapi=np.moveaxis(dapi,(y_axis,x_axis),(0,1))
    return np.asarray(dapi), full


def clean_components(mask, min_area, max_area):
    mask=np.uint8(mask>0)
    mask=cv2.morphologyEx(mask, cv2.MORPH_OPEN, np.ones((3,3),np.uint8))
    n, labels, stats, cents=cv2.connectedComponentsWithStats(mask, 8)
    keep=(stats[:,cv2.CC_STAT_AREA]>=min_area)&(stats[:,cv2.CC_STAT_AREA]<=max_area)
    keep[0]=False
    ids=np.flatnonzero(keep)
    return cents[ids], labels, stats[ids,cv2.CC_STAT_AREA]


def dapi_centroids(a, min_area, max_area):
    x=a.astype(np.float32); nz=x[x>0]; lo,hi=np.percentile(nz,[1,99.7]) if nz.size else (0,1)
    u=np.uint8(np.clip((x-lo)/max(hi-lo,1),0,1)*255)
    u=cv2.GaussianBlur(u,(0,0),1)
    _,m=cv2.threshold(u,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
    return clean_components(m,min_area,max_area)


def he_centroids(rgb, min_area, max_area):
    rgb=np.asarray(rgb).astype(np.uint8)
    gray=cv2.cvtColor(rgb,cv2.COLOR_RGB2GRAY)
    # Hematoxylin nuclei are usually blue/purple; combine blue excess with darkness.
    r,g,b=[rgb[...,i].astype(np.float32) for i in range(3)]
    blue_excess=np.clip(b-0.5*(r+g)+128,0,255).astype(np.uint8)
    dark=255-gray
    _,m1=cv2.threshold(blue_excess,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
    _,m2=cv2.threshold(dark.astype(np.uint8),0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
    m=cv2.bitwise_and(m1,m2)
    return clean_components(m,min_area,max_area)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--he-registered',required=True,help='Registered HE image in SP coordinates')
    p.add_argument('--sp',required=True,help='Original SP/CODEX image containing DAPI')
    p.add_argument('--out-dir',required=True)
    p.add_argument('--sp-channel',type=int,default=0)
    p.add_argument('--level',type=int,default=3,help='Common pyramid level used for both images; choose a level with matching HE/SP dimensions (default: 3)')
    p.add_argument('--min-area',type=int,default=8,help='Minimum object area at the lowest pyramid level')
    p.add_argument('--max-area',type=int,default=2000)
    p.add_argument('--max-distance',type=float,default=25,help='Pass threshold in lowest-level pixels')
    p.add_argument('--pixel-size-um',type=float,default=None,help='Pixel size at the analyzed pyramid level')
    a=p.parse_args(); out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True)
    he,he_full=read_level(a.he_registered,a.level,he=True); sp,sp_full=read_level(a.sp,a.level,a.sp_channel,he=False)
    if he.shape[:2]!=sp.shape[:2]:
        raise ValueError(f'HE/SP dimensions differ at level {a.level}: HE={he.shape[:2]}, SP={sp.shape[:2]}. Choose a common level with matching dimensions.')
    resize_note='no resize; same analyzed thumbnail grid'
    dapi,_,da=dapi_centroids(sp,a.min_area,a.max_area)
    hen,_,ha=he_centroids(he,a.min_area,a.max_area)
    if len(dapi)==0 or len(hen)==0: raise RuntimeError(f'No nuclei detected: DAPI={len(dapi)}, HE={len(hen)}; adjust thresholds/area')
    dist,idx=cKDTree(dapi).query(hen)
    reverse_dist,reverse_idx=cKDTree(hen).query(dapi)
    unit='px'; scale=1.0
    if a.pixel_size_um: dist=dist*a.pixel_size_um; reverse_dist=reverse_dist*a.pixel_size_um; scale=a.pixel_size_um; unit='um'
    order=np.argsort(dist)
    with (out/'centroid_pairs.csv').open('w',newline='') as f:
        w=csv.writer(f); w.writerow(['he_x','he_y','dapi_x','dapi_y',f'nearest_distance_{unit}','pass'])
        for i in order:
            w.writerow([hen[i,0],hen[i,1],dapi[idx[i],0],dapi[idx[i],1],dist[i],int(dist[i]<=a.max_distance*scale)])
    summary={'he_detected':int(len(hen)),'dapi_detected':int(len(dapi)),
             'he_to_dapi_median_distance_'+unit:float(np.median(dist)),
             'he_to_dapi_mean_distance_'+unit:float(np.mean(dist)),'he_to_dapi_p90_distance_'+unit:float(np.percentile(dist,90)),
             'dapi_to_he_median_distance_'+unit:float(np.median(reverse_dist)),
             'dapi_to_he_mean_distance_'+unit:float(np.mean(reverse_dist)),'dapi_to_he_p90_distance_'+unit:float(np.percentile(reverse_dist,90)),
             'pass_threshold_'+unit:float(a.max_distance*scale),
             'he_to_dapi_pass_fraction':float(np.mean(dist<=a.max_distance*scale)),
             'dapi_to_he_pass_fraction':float(np.mean(reverse_dist<=a.max_distance*scale)),
             'he_full_hw':he_full,'sp_full_hw':sp_full,'pyramid_level':a.level,'he_analyzed_shape':list(he.shape[:2]),'sp_analyzed_shape':list(sp.shape[:2]),'resize_note':resize_note,
             'warning':'This is a screening QC metric; it depends on HE/DAPI segmentation quality and pyramid level.'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)); print(json.dumps(summary,indent=2))
    he_vis=he.copy(); dapi_u=np.uint8(np.clip((sp.astype(np.float32)-np.percentile(sp,[1])[0])/max(float(np.percentile(sp,[99])[0]-np.percentile(sp,[1])[0]),1),0,1)*255); dapi_vis=cv2.cvtColor(dapi_u,cv2.COLOR_GRAY2BGR)
    for x,y in hen: cv2.circle(he_vis,(int(round(x)),int(round(y))),3,(255,0,0),1)
    for x,y in dapi: cv2.circle(dapi_vis,(int(round(x)),int(round(y))),3,(0,0,255),1)
    cv2.imwrite(str(out/'he_detection.jpg'),cv2.cvtColor(he_vis,cv2.COLOR_RGB2BGR)); cv2.imwrite(str(out/'dapi_detection.jpg'),dapi_vis)

if __name__=='__main__': main()
