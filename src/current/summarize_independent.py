#!/usr/bin/env python3
"""Descriptive paired summary. Seeds are technical replicates, not donors."""
from pathlib import Path
import json,csv,sys,statistics
root=Path(sys.argv[1]) if len(sys.argv)>1 else Path('.')
donor_map={'151507':'Donor1','151508':'Donor1','151509':'Donor1','151669':'Donor2','151670':'Donor2','151671':'Donor2','151673':'Donor3','151674':'Donor3','151675':'Donor3'}
rows=[]
for path in sorted(root.glob('pair_*_seed*_*.json')):
 d=json.loads(path.read_text());sc=d['scores'];pre=sc['baseline'];r={'donor':donor_map.get(d['source'],'Unknown'),'pair':d['source']+'->'+d['target'],'model':d['model'],'seed':d['seed'],'n_target_labeled':pre['n_evaluable'],'baseline_lt':pre['LT'],'baseline_balacc':pre['balanced_accuracy']}
 for label,key in [('auto','Auto_conditional'),('no_zca','Auto_procrustes_only'),('always_zca','Auto_always_ZCA'),('oracle','Oracle_conditional')]:
  r[label+'_lt']=sc[key]['metrics']['LT'];r[label+'_balacc']=sc[key]['metrics']['balanced_accuracy']
 r['delta_auto']=r['auto_lt']-r['baseline_lt'];r['delta_no_zca']=r['no_zca_lt']-r['baseline_lt']
 r['covariance_distance']=sc['Auto_conditional']['diagnostics']['covariance_distance']
 r['zca_triggered']=sc['Auto_conditional']['diagnostics']['whitened']
 rows.append(r)
if not rows:raise SystemExit('No pair result files')
with (root/'PAIR_SEED_SUMMARY.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
metrics=['baseline_lt','baseline_balacc','auto_lt','auto_balacc','no_zca_lt','no_zca_balacc','always_zca_lt','oracle_lt','delta_auto','delta_no_zca']
groups={}
for r in rows:groups.setdefault((r['donor'],r['pair'],r['model']),[]).append(r)
paired=[]
for (donor,pair,model),a in sorted(groups.items()):
 d={'donor':donor,'pair':pair,'model':model,'n_seeds':len(a),'seeds':','.join(str(x['seed']) for x in a)}
 for k in metrics:d[k+'_seed_mean']=statistics.mean(x[k] for x in a)
 paired.append(d)
with (root/'PAIRED_MEANS.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(paired[0]));w.writeheader();w.writerows(paired)
summary={}
for donor in ['Donor1','Donor2','Donor3']:
 for model in ['graphsage','pca']:
  records=[r for r in paired if r['donor']==donor and r['model']==model]
  if records:
   summary[f'{donor}_{model}']={
     'number_of_pairs':len(records),
     'number_of_seed_pair_records':len([r for r in rows if r['donor']==donor and r['model']==model]),
     'baseline_LT_pair_mean':statistics.mean(r['baseline_lt_seed_mean'] for r in records),
     'Auto_LT_pair_mean':statistics.mean(r['auto_lt_seed_mean'] for r in records),
     'Auto_delta_pair_mean':statistics.mean(r['delta_auto_seed_mean'] for r in records),
     'paired_changes':[{'pair':r['pair'],'delta':r['delta_auto_seed_mean']} for r in records],
     'CAVEAT':'Not independent biological replicates; donor n=1 per entry.'
   }
(root/'SUMMARY_DESCRIPTIVE.json').write_text(json.dumps(summary,indent=2))
print('Read',len(rows),'pair-seed-model records and',len(paired),'pair-model summaries')
print(json.dumps(summary,indent=2))
if len(rows)!=54:
 print('CAUTION: full design expects 9 pairs x 3 seeds x 2 models = 54 pair-seed-model records; current records:',len(rows),file=sys.stderr)
