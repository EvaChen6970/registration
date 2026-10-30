#!/usr/bin/env python3
import argparse, csv, json
from pathlib import Path

p=argparse.ArgumentParser(); p.add_argument('--root',required=True); a=p.parse_args(); root=Path(a.root); reg=root/'registration'
manifest=reg/'registration_manifest.csv'
rows=list(csv.DictReader(manifest.open())) if manifest.exists() else []
imgs=list((reg/'registered_h&e').rglob('*.ome.tif'))
trans=list((reg/'transforms').rglob('*_HE_to_CODEX.json'))
print(f'manifest_pairs={len(rows)} registered_outputs={len(imgs)} transforms={len(trans)}')
from collections import Counter
print('manifest_statuses=',dict(Counter(r.get('status','') for r in rows)))
print('missing_output_rows:')
for r in rows:
 if r.get('status','').startswith('registered'):
  sid=r.get('sample',''); found=list((reg/'registered_h&e').rglob(f'{sid}_HE_to_CODEX.ome.tif'))
  if not found: print(sid, 'expected', reg/'registered_h&e'/f'{sid}_HE_to_CODEX.ome.tif')
print('rejected_or_failed:')
for r in rows:
 if not r.get('status','').startswith('registered'): print(r.get('sample'),r.get('status'),r.get('error',''))
summary=reg/'FINAL_QC_SUMMARY.json'
if summary.exists(): print('FINAL_QC_SUMMARY=',json.dumps(json.loads(summary.read_text()),indent=2))
