#!/usr/bin/env python3
"""Compare two integer instance-label masks already sampled on the same pixel grid."""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import tifffile


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--he-labels', required=True, help='Instance-label TIFF segmented from registered HE in SP coordinates')
    p.add_argument('--sp-labels', required=True, help='SP/CODEX instance-label TIFF')
    p.add_argument('--out-prefix', required=True)
    p.add_argument('--chunk-rows', type=int, default=512)
    p.add_argument('--match-iou', type=float, default=0.1)
    a = p.parse_args()
    ha, sa = tifffile.memmap(a.he_labels), tifffile.memmap(a.sp_labels)
    if ha.shape != sa.shape or ha.ndim != 2:
        p.error(f'label masks must be 2D with identical shape; got HE={ha.shape}, SP={sa.shape}')
    max_sp = int(np.max(sa))
    if max_sp < 0:
        p.error('SP labels contain no nonnegative background/instance labels')
    base = max_sp + 1
    he_area, sp_area = defaultdict(int), defaultdict(int)
    intersections = defaultdict(int)
    for y in range(0, ha.shape[0], a.chunk_rows):
        h = np.asarray(ha[y:y+a.chunk_rows]).astype(np.int64, copy=False).ravel()
        s = np.asarray(sa[y:y+a.chunk_rows]).astype(np.int64, copy=False).ravel()
        hv, hc = np.unique(h[h > 0], return_counts=True)
        sv, sc = np.unique(s[s > 0], return_counts=True)
        he_area.update(dict(zip(hv.tolist(), hc.tolist())))
        sp_area.update(dict(zip(sv.tolist(), sc.tolist())))
        both = (h > 0) & (s > 0)
        if both.any():
            keys, counts = np.unique(h[both] * base + s[both], return_counts=True)
            for key, count in zip(keys.tolist(), counts.tolist()):
                intersections[key] += count
    candidates = []
    for key, inter in intersections.items():
        hid, sid = divmod(key, base)
        union = he_area[hid] + sp_area[sid] - inter
        iou = inter / union
        candidates.append((iou, hid, sid, inter, union))
    candidates.sort(reverse=True)
    used_h, used_s, matches = set(), set(), []
    for iou, hid, sid, inter, union in candidates:
        if iou < a.match_iou: break
        if hid not in used_h and sid not in used_s:
            used_h.add(hid); used_s.add(sid)
            matches.append((hid, sid, inter, union, iou, 2*inter/(he_area[hid]+sp_area[sid])))
    prefix = Path(a.out_prefix); prefix.parent.mkdir(parents=True, exist_ok=True)
    with Path(str(prefix) + '_pairs.csv').open('w', newline='') as f:
        w = csv.writer(f); w.writerow(['he_label','sp_label','intersection_px','union_px','iou','dice'])
        w.writerows(matches)
    all_he, all_sp = set(he_area), set(sp_area)
    summary = {
        'shape': list(ha.shape), 'he_nuclei': len(all_he), 'sp_nuclei': len(all_sp),
        'matched_pairs': len(matches), 'match_iou_threshold': a.match_iou,
        'he_match_fraction': len(used_h)/max(len(all_he),1),
        'sp_match_fraction': len(used_s)/max(len(all_sp),1),
        'mean_matched_iou': float(np.mean([m[4] for m in matches])) if matches else 0.0,
        'median_matched_iou': float(np.median([m[4] for m in matches])) if matches else 0.0,
        'mean_matched_dice': float(np.mean([m[5] for m in matches])) if matches else 0.0,
        'interpretation': 'One-to-one instance matches selected greedily by descending IoU; pixel-grid agreement is required.'
    }
    Path(str(prefix) + '_summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

if __name__ == '__main__': main()
