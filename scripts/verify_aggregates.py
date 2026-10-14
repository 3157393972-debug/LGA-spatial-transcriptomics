#!/usr/bin/env python3
"""Standalone sanity checks on archived aggregate CSVs; NOT training."""
import csv
from pathlib import Path
from collections import defaultdict
root=Path(__file__).resolve().parents[1]
def read(filename):
    with (root/'results/verified'/filename).open(newline='') as h:
        return list(csv.DictReader(h))
rows=read('CANONICAL_ALL_54.csv')
assert len(rows)==54, len(rows)
assert {r['model'] for r in rows}=={'graphsage','pca'}
seen={(r['model'],r['source'],r['target'],r['seed']) for r in rows}
assert len(seen)==54
for model,expected_base,expected_auto in [('graphsage',0.2857,0.4528),('pca',0.4325,0.4733)]:
    group=[r for r in rows if r['model']==model]
    base=sum(float(r['baseline']) for r in group)/len(group)
    after=sum(float(r['Auto_conditional']) for r in group)/len(group)
    assert abs(base-expected_base)<0.0002,(model,'baseline',base)
    assert abs(after-expected_auto)<0.0002,(model,'Auto',after)
    print(f'{model}: n={len(group)}, baseline={base:.6f}, fixed-Auto={after:.6f}')
r6=read('FIG6_108_ROTATIONS.csv')
assert len(r6)==108, len(r6)
print('figure6:',len(r6),'rotation records (81 GraphSAGE + 27 PCA)')
print('Descriptive CSV verification PASSED; this does not verify raw-data reproduction.')
