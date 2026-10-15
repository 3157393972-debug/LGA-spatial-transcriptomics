#!/usr/bin/env python3
"""Unified, actual-latent independent PCA/GraphSAGE controlled-rotation evaluation for Figure 6.
Uses fixed K=7, pooled depth quantiles, tau=0.1, canonical minimum-change orthogonal completion,
target-to-source LT with held-out anatomical labels used only in evaluation.
"""
import os,re,csv,json,hashlib,importlib.util
from pathlib import Path
from collections import defaultdict
import numpy as np
ROOT=Path(os.environ.get('BIB_FULL9_ROOT','./BIB_PAPER1_INDEPENDENT_FULL_9PAIRS')).expanduser().resolve();OUT=Path(os.environ.get('BIB_FIG6_OUT','./fig6')).expanduser().resolve();OUT.mkdir(parents=True,exist_ok=True)
from pathlib import Path
_REPO_SOURCE = Path(__file__).resolve().parent / 'paper1_independent_jsonfixed.py'
_ARCHIVED_SOURCE = ROOT / 'paper1_independent_jsonfixed.py'
_SOURCE = _REPO_SOURCE if _REPO_SOURCE.is_file() else _ARCHIVED_SOURCE
if not _SOURCE.is_file():
    raise FileNotFoundError(f'No LGA independent source at {_REPO_SOURCE} or {_ARCHIVED_SOURCE}')
spec=importlib.util.spec_from_file_location('original',_SOURCE);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
pairs=[('151507','151508'),('151508','151509'),('151509','151510'),('151669','151670'),('151670','151671'),('151671','151672'),('151673','151674'),('151674','151675'),('151675','151676')]

def donor_of(s):return 1 if s in {'151507','151508','151509'} else 2 if s in {'151669','151670','151671'} else 3
cache={}
def load(sid,seed,model):
 k=(sid,seed,model)
 if k not in cache:
  with np.load(ROOT/f'{sid}_seed{seed}_{model}_independent.npz',allow_pickle=False) as d:cache[k]={x:d[x] for x in ['H','coords','labels']}
 return cache[k]
def make_rotation(n,seed_str):
 seed=int(hashlib.sha256(seed_str.encode()).hexdigest()[:8],16)
 g=np.random.default_rng(seed).normal(size=(n,n));Q,R=np.linalg.qr(g);Q*=np.where(np.diag(R)>=0,1,-1)[None,:]
 return Q,seed
rows=[]
for model in ['graphsage','pca']:
 for s,t in pairs:
  for seed in ([42,123,456] if model=='graphsage' else [42]):
   ds=load(s,seed,model);dt=load(t,seed,model)
   hs=ds['H'];ht=dt['H']
   axis=m.axis_from_source(ds['coords']);gs,gt,edges=m.auto_groups(ds['coords'],dt['coords'],axis,K=7,scheme='pooled_quantile')
   baseline=m.full_target_lt(hs,ht,ds['labels'],dt['labels'])['LT']
   for rotation_i in range(3):
    rot,rot_seed=make_rotation(hs.shape[1],f'{model}:{s}->{t}:seed={seed}:rotation={rotation_i}:v1')
    hsr=hs.astype(float)@rot
    rotlt=m.full_target_lt(hsr,ht,ds['labels'],dt['labels'])['LT']
    A,T,diag=m.harmonize(hsr,ht,gs,gt,whitening='conditional',tau=.1,completion='min_change_sensitivity_only')
    autolt=m.full_target_lt(A,T,ds['labels'],dt['labels'])['LT']
    row={'model':model,'donor':donor_of(s),'pair':f'{s}->{t}','seed':seed,'rotation':rotation_i,'rotation_seed':rot_seed,'baseline_LT':baseline,'rotated_LT':rotlt,'auto_LT':autolt,'delta_LT':autolt-rotlt,'rank':diag['rank'],'whitened':diag['whitened']}
    rows.append(row)
  print(model,s,t,'complete',len(rows),flush=True)
with (OUT/'FIG6_108_ROTATIONS.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
paired=defaultdict(list)
for r in rows:paired[(r['model'],r['donor'],r['pair'])].append(r)
pr=[]
for (model,donor,pair),rs in sorted(paired.items()):
 pr.append({'model':model,'donor':donor,'pair':pair,'n_rotated':len(rs),**{k:float(np.mean([x[k] for x in rs])) for k in ['baseline_LT','rotated_LT','auto_LT','delta_LT']},'fraction_improved':float(np.mean([x['delta_LT']>0 for x in rs]))})
with (OUT/'FIG6_18_PAIR_MEANS.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(pr[0]));w.writeheader();w.writerows(pr)
agg=[]
for model in ['graphsage','pca']:
 for donor in [1,2,3]:
  vals=[p for p in pr if p['model']==model and p['donor']==donor]
  item={'model':model,'donor':donor,'n_pairs':len(vals),'n_rotated_records':sum(x['n_rotated'] for x in vals)}
  for k in ['baseline_LT','rotated_LT','auto_LT','delta_LT']:
   v=[p[k] for p in vals];item[k+'_mean']=float(np.mean(v));item[k+'_sd_pairs']=float(np.std(v,ddof=1))
  agg.append(item)
for model in ['graphsage','pca']:
 vals=[p for p in pr if p['model']==model]
 print('OVERALL',model, 'rotated',np.mean([p['rotated_LT'] for p in vals]),'auto',np.mean([p['auto_LT'] for p in vals]),'delta',np.mean([p['delta_LT'] for p in vals]),'improved pairs',sum(p['delta_LT']>0 for p in vals),flush=True)
with (OUT/'FIG6_DONOR_MEANS.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(agg[0]));w.writeheader();w.writerows(agg)
print('DONOR means',agg)
