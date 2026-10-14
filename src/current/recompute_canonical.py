#!/usr/bin/env python3
"""Re-evaluate archived independently generated DLPFC latents with minimum-change optimal orthogonal completion.
Read-only relative to source archives. Checks historical default SVD replay first.
"""
import os, sys, re, json, csv, importlib.util, hashlib
from pathlib import Path
from collections import defaultdict
import numpy as np
ROOT=Path(os.environ.get('BIB_FULL9_ROOT','./BIB_PAPER1_INDEPENDENT_FULL_9PAIRS')).expanduser().resolve()
OUT=Path(os.environ.get('BIB_CANONICAL_OUT','./canonical')).expanduser().resolve();OUT.mkdir(parents=True,exist_ok=True)
spec=importlib.util.spec_from_file_location('historical_independent',ROOT/'paper1_independent_jsonfixed.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
orig_harmonize=m.harmonize
cached={}
def load(sid,seed,model):
    key=(sid,seed,model)
    if key not in cached:
        with np.load(ROOT/f'{sid}_seed{seed}_{model}_independent.npz',allow_pickle=False) as x:
            cached[key]={'H':x['H'], 'coords':x['coords'],'labels':x['labels']}
    return cached[key]

def canon_harmonize(*args,**kwargs):
    kwargs['completion']='min_change_sensitivity_only'
    return orig_harmonize(*args,**kwargs)

records=[]; comparisons=[]; all_min_q=[]
pattern=re.compile(r'pair_(\d+)_(\d+)_seed(\d+)_(pca|graphsage)\.json')
paths=sorted(ROOT.glob('pair_*.json'))
assert len(paths)==54, len(paths)
for ii,fp in enumerate(paths):
    mt=pattern.fullmatch(fp.name);assert mt,fp
    s,t,seed,model=mt.groups();seed=int(seed)
    with open(fp) as f:old=json.load(f)
    ds=load(s,seed,model);dt=load(t,seed,model)
    m.harmonize=orig_harmonize
    exact=m.posthoc_evaluation(ds['H'],dt['H'],ds,dt,K=7,tau=0.1)
    for key in ['baseline','Auto_conditional','Auto_procrustes_only','Auto_always_ZCA','Oracle_conditional']:
        original=old['scores'][key]
        replay=exact[key]
        oo=original['LT'] if key=='baseline' else original['metrics']['LT']
        rr=replay['LT'] if key=='baseline' else replay['metrics']['LT']
        if abs(oo-rr)>1e-10: comparisons.append((fp.name,key,oo,rr))
    m.harmonize=canon_harmonize
    canon=m.posthoc_evaluation(ds['H'],dt['H'],ds,dt,K=7,tau=0.1)
    row={'pair':s+'->'+t,'source':s,'target':t,'seed':seed,'model':model,'donor':1 if s in ['151507','151508','151509'] else 2 if s in ['151669','151670','151671'] else 3,'n_evaluable':exact['baseline']['n_evaluable']}
    row['baseline']=canon['baseline']['LT']
    for key in ['Auto_conditional','Auto_procrustes_only','Auto_always_ZCA','Oracle_conditional']:
        row[key]=canon[key]['metrics']['LT'];row[key+'_balacc']=canon[key]['metrics']['balanced_accuracy'];row[key+'_old']=old['scores'][key]['metrics']['LT'];row[key+'_change']=row[key]-row[key+'_old']
        row[key+'_rank']=canon[key]['diagnostics']['rank']
    row['baseline_balacc']=canon['baseline']['balanced_accuracy']
    row['delta_auto']=row['Auto_conditional']-row['baseline']
    records.append(row)
    if (ii+1)%10==0:print('Processed',ii+1,'/',len(paths),flush=True)

assert not comparisons, f'Default pipeline does not replicate prior score; first errors: {comparisons[:5]}'
keys=list(records[0])
with open(OUT/'CANONICAL_ALL_54.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(records)
by=defaultdict(list)
for r in records:
    by[(r['model'],r['donor'])].append(r)
summary=[]
for (model,donor),rs in sorted(by.items()):
    item={'model':model,'donor':donor,'n_unique_pairs':len(set(x['pair'] for x in rs)),'n_seed_records':len(rs)}
    for k in ['baseline','Auto_conditional','Auto_procrustes_only','Auto_always_ZCA','Oracle_conditional','baseline_balacc','Auto_conditional_balacc','delta_auto']:
        item[k]=float(np.mean([x[k] for x in rs]))
    item['n_improved']=sum(x['delta_auto']>1e-12 for x in rs)
    summary.append(item)
with open(OUT/'CANONICAL_DONOR_MEANS.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
for model in ['graphsage','pca']:
    rs=[r for r in records if r['model']==model]
    print(model, 'base',np.mean([r['baseline'] for r in rs]), 'auto',np.mean([r['Auto_conditional'] for r in rs]),'procr',np.mean([r['Auto_procrustes_only'] for r in rs]),'oracle',np.mean([r['Oracle_conditional'] for r in rs]),'improved',sum(r['delta_auto']>0 for r in rs), 'n',len(rs))
for r in summary: print('DONOR',r)
info={'n_records':len(records),'default_replay_mismatch':len(comparisons),'completion':'min_change_sensitivity_only','rank_cutoff_rel':1e-6,'source_snapshot':str(ROOT),'source_py_sha256':hashlib.sha256((ROOT/'paper1_independent_jsonfixed.py').read_bytes()).hexdigest(),'summary':summary}
(OUT/'RUN_MANIFEST.json').write_text(json.dumps(info,indent=2))
print('DONE',OUT)
