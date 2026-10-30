#!/usr/bin/env python3
"""Register one HE WSI to one SP/CODEX WSI using pyramidal SIFT and affine RANSAC."""
import argparse
import json
import subprocess
from pathlib import Path

import cv2
import numpy as np
import tifffile


def normalize_u8(a):
    a = np.asarray(a)
    nz = a[a > 0]
    lo, hi = np.percentile(nz, [1, 99.7]) if nz.size else (0, 1)
    return np.uint8(np.clip((a.astype(np.float32) - lo) / max(hi - lo, 1), 0, 1) * 255)


def read_thumbnail(path, channel, is_he):
    with tifffile.TiffFile(path) as tif:
        series = tif.series[0]
        full_shape = tuple(series.shape)
        thumb = series.levels[-1].asarray()
    if is_he:
        if thumb.ndim == 2:
            rgb = cv2.cvtColor(normalize_u8(thumb), cv2.COLOR_GRAY2RGB)
        else:
            rgb = np.asarray(thumb[..., :3])
        gray = normalize_u8(255 - cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY))
        full_hw = full_shape[:2]
    else:
        if thumb.ndim == 2:
            plane = thumb
        elif thumb.ndim >= 3:
            plane = thumb[channel]
            while plane.ndim > 2:
                plane = plane[0]
        else:
            raise ValueError(f'Unsupported SP thumbnail shape: {thumb.shape}')
        gray = normalize_u8(plane)
        rgb = cv2.cvtColor(gray, cv2.COLOR_GRAY2RGB)
        full_hw = full_shape[-2:]
    return rgb, gray, tuple(map(int, full_hw))


def structural(a):
    a = cv2.GaussianBlur(a, (0, 0), 2)
    a = cv2.morphologyEx(a, cv2.MORPH_GRADIENT, np.ones((7, 7), np.uint8))
    return cv2.normalize(a, None, 0, 255, cv2.NORM_MINMAX)


def estimate(fixed, moving, max_features):
    sift = cv2.SIFT_create(nfeatures=max_features, contrastThreshold=0.01, edgeThreshold=20)
    kf, df = sift.detectAndCompute(structural(fixed), None)
    km, dm = sift.detectAndCompute(structural(moving), None)
    if df is None or dm is None:
        raise RuntimeError('No SIFT descriptors found')
    pairs = cv2.BFMatcher().knnMatch(dm, df, k=2)
    good = [m for m, n in pairs if m.distance < 0.75 * n.distance]
    if len(good) < 8:
        raise RuntimeError(f'Only {len(good)} cross-modal SIFT matches')
    src = np.float32([km[m.queryIdx].pt for m in good])
    dst = np.float32([kf[m.trainIdx].pt for m in good])
    mat, inliers = cv2.estimateAffine2D(src, dst, method=cv2.RANSAC, ransacReprojThreshold=8,
                                        maxIters=20000, confidence=0.999, refineIters=30)
    if mat is None or int(inliers.sum()) < 5:
        raise RuntimeError('Affine estimation failed or had fewer than 5 inliers')
    mask = inliers.ravel().astype(bool)
    err = np.linalg.norm(src @ mat[:, :2].T + mat[:, 2] - dst, axis=1)
    return mat, len(good), int(mask.sum()), float(np.median(err[mask]))


def full_matrix(mat, fixed_thumb_hw, moving_thumb_hw, fixed_full_hw, moving_full_hw):
    sf = np.diag([fixed_full_hw[1] / fixed_thumb_hw[1], fixed_full_hw[0] / fixed_thumb_hw[0], 1.0])
    smi = np.diag([moving_thumb_hw[1] / moving_full_hw[1], moving_thumb_hw[0] / moving_full_hw[0], 1.0])
    return (sf @ np.vstack([mat, [0, 0, 1]]) @ smi)[:2]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--he', required=True, help='Moving H&E WSI (OME-TIFF/QPTIFF readable by tifffile and libvips)')
    p.add_argument('--sp', required=True, help='Fixed spatial-proteomics WSI, e.g. CODEX QPTIFF')
    p.add_argument('--out-dir', required=True)
    p.add_argument('--sp-channel', type=int, default=0, help='SP pyramid channel used for feature matching (default: 0/DAPI)')
    p.add_argument('--max-features', type=int, default=12000)
    p.add_argument('--vips', default='vips')
    p.add_argument('--overwrite', action='store_true')
    a = p.parse_args()
    he, sp = Path(a.he).resolve(), Path(a.sp).resolve()
    if not he.is_file() or not sp.is_file(): p.error('Both --he and --sp must exist')
    out = Path(a.out_dir).resolve(); out.mkdir(parents=True, exist_ok=True)
    stem = he.stem + '_to_' + sp.stem
    registered = out / f'{stem}_HE_to_SP.ome.tif'
    transform = out / f'{stem}_HE_to_SP.json'
    overlay_path = out / f'{stem}_overlay.jpg'
    if registered.exists() and not a.overwrite: p.error(f'{registered} exists; use --overwrite')
    sp_rgb, fixed, sp_full = read_thumbnail(sp, a.sp_channel, False)
    he_rgb, moving, he_full = read_thumbnail(he, 0, True)
    thumb_mat, matches, inliers, med = estimate(fixed, moving, a.max_features)
    full_mat = full_matrix(thumb_mat, fixed.shape, moving.shape, sp_full, he_full)
    det = float(np.linalg.det(full_mat[:, :2]))
    if det <= 0: raise RuntimeError(f'Reflected/degenerate affine transform (det={det})')
    warped = cv2.warpAffine(he_rgb, thumb_mat, (fixed.shape[1], fixed.shape[0]), borderValue=(255,255,255))
    overlay = cv2.addWeighted(sp_rgb, .5, warped, .5, 0)
    contact = np.hstack([sp_rgb, warped, overlay])
    cv2.imwrite(str(overlay_path), cv2.cvtColor(contact, cv2.COLOR_RGB2BGR), [cv2.IMWRITE_JPEG_QUALITY, 92])
    a0,b0,tx=full_mat[0]; c0,d0,ty=full_mat[1]
    vips_dst=f'{registered}[tile,pyramid,subifd,bigtiff,compression=deflate,tile-width=512,tile-height=512]'
    cmd=[a.vips,'affine',str(he),vips_dst,f'{a0} {b0} {c0} {d0}','--odx',str(tx),'--ody',str(ty),
         '--oarea',f'0 0 {sp_full[1]} {sp_full[0]}','--interpolate','bicubic','--background','255 255 255']
    subprocess.run(cmd,check=True)
    payload={'moving_he':str(he),'fixed_sp':str(sp),'mapping':'HE full-resolution pixel coordinates to SP full-resolution pixel coordinates',
             'sp_channel':a.sp_channel,'thumbnail_affine_2x3':thumb_mat.tolist(),'full_resolution_affine_2x3':full_mat.tolist(),
             'sp_full_hw':sp_full,'he_full_hw':he_full,'matches':matches,'inliers':inliers,
             'inlier_fraction':inliers/matches,'median_inlier_error_thumbnail_px':med,'qc_overlay':str(overlay_path)}
    transform.write_text(json.dumps(payload,indent=2))
    print(f'registered={registered}\ntransform={transform}\nqc_overlay={overlay_path}\nmatches={matches} inliers={inliers} median_error_thumbnail_px={med:.3f}')

if __name__ == '__main__': main()
