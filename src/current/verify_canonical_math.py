import numpy as np,glob,os,importlib.util,json
from pathlib import Path
B=Path(os.environ.get('BIB_FULL9_ROOT','./BIB_PAPER1_INDEPENDENT_FULL_9PAIRS')).expanduser().resolve()
spec=importlib.util.spec_from_file_location('original',B/'paper1_independent_jsonfixed.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
max_obj_gap=max_orth=max_change_violation=0.0
min_Q_singular=1.0
ranks={}
count=0
for f in sorted(B.glob('pair_*.json')):
 x=json.load(open(f));s=x['source'];t=x['target'];seed=x['seed'];model=x['model']
 a=np.load(B/f'{s}_seed{seed}_{model}_independent.npz');b=np.load(B/f'{t}_seed{seed}_{model}_independent.npz');hs=a['H'];ht=b['H'];axis=m.axis_from_source(a['coords']);gs,gt,_=m.auto_groups(a['coords'],b['coords'],axis,K=7)
 if m.covariance_distance(hs,ht)>.1:S=m.zca(hs);T=m.zca(ht)
 else:S=hs.astype(float)-hs.mean(0);T=ht.astype(float)-ht.mean(0)
 common=[int(j) for j in sorted(set(gs)&set(gt)) if (gs==j).sum()>=20 and (gt==j).sum()>=20]
 cs=np.array([S[gs==j].mean(0) for j in common]);ct=np.array([T[gt==j].mean(0) for j in common])
 M=cs.T@ct;U,sv,Vt=np.linalg.svd(M,full_matrices=True)
 r=int(np.sum(sv>1e-6*sv[0])) if sv[0]>0 else 0
 u0=U[:,r:];v0=Vt[r:,:].T;Q=v0.T@u0; q=np.linalg.svd(Q,compute_uv=False)
 min_Q_singular=min(min_Q_singular,float(q[-1]))
 Rc,rc=m.canonical_rotation(cs,ct);Rd=U@Vt
 errorth=np.linalg.norm(Rc.T@Rc-np.eye(Rc.shape[0]),'fro')
 obj=float(np.linalg.norm(cs@Rc-ct,'fro')**2-np.linalg.norm(cs@Rd-ct,'fro')**2)
 ch=float(np.linalg.norm(Rc-np.eye(Rc.shape[0]),'fro')-np.linalg.norm(Rd-np.eye(Rd.shape[0]),'fro'))
 max_obj_gap=max(max_obj_gap,abs(obj));max_orth=max(max_orth,errorth);max_change_violation=max(max_change_violation,ch)
 ranks[r]=ranks.get(r,0)+1;count+=1
result={'n_records':count,'rank_distribution':ranks,'max_abs_procrustes_obj_gap':max_obj_gap,'max_orthogonality_error':max_orth,'max_minimum_change_violation':max_change_violation,'minimum_Q_singular_value':min_Q_singular}
assert count==54 and max_obj_gap<1e-8 and max_orth<1e-10 and max_change_violation<1e-10 and min_Q_singular>1e-8
(Path(os.environ.get('BIB_CANONICAL_OUT','./canonical')).expanduser().resolve()/'CANONICAL_MATH_TEST.json').write_text(json.dumps(result,indent=2))
print(result)
