#!/usr/bin/env python3
import json,sys,math
from pathlib import Path
p=Path(sys.argv[1]);d=json.loads(p.read_text())
required=('lt_orig','lt_unrot_lga','lt_rot','lt_lga')
recs=d['results']
assert len(recs)==3 and [r['seed'] for r in recs]==[42,123,456]
for key in required:
    vals=[float(r[key]) for r in recs]
    assert all(math.isfinite(v) and 0<=v<=1 for v in vals),(key,vals)
    mu=sum(vals)/len(vals)
    assert abs(mu-float(d['aggregate'][key]))<1e-10,(key,mu,d['aggregate'][key])
print('MOUSE RESULT VALID:',p,'aggregate',d['aggregate'])
# These are historical expectations, not exact-fail assertions across machines:
old={'lt_orig':0.1777539341917024,'lt_rot':0.1579640450212869,'lt_lga':0.6164520743919885}
for key,expected in old.items():
    delta=float(d['aggregate'][key])-expected
    print('Historical comparison:',key,'difference=',delta)
    if abs(delta)>0.02:print('WARNING: historical deviation >0.02; investigate before scientific interpretation')
print('CRITICAL: this remains pseudo-label validation, NOT independent anatomical truth.')
