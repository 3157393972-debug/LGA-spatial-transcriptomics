#!/usr/bin/env python3
"""BIB Paper 1: prospective independent-embedding experiment.

Purpose: directly address editor's concern. Every section is preprocessed and
trained separately; NO MNN, NO pair-trained contrastive objective, NO shared
encoder weights, NO label use for encoder fitting or Auto correspondence.

This is a NEW experimental arm; it does NOT reproduce historical joint GraphSAGE
training. All numeric outcomes must be obtained by executing on raw DLPFC h5ad.

Requirements on AutoDL: scanpy, scikit-learn, scipy, torch, numpy, pandas.
Run from anywhere; historical author code not overwritten.
"""
from __future__ import annotations
import argparse, gc, hashlib, json, os, platform, random, sys
from pathlib import Path
import numpy as np
from scipy import sparse as sp
from scipy.stats import rankdata
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import normalize

LAYERS=['L1','L2','L3','L4','L5','L6','WM']
LABEL_MAP={f'Layer_{j}':f'L{j}' for j in range(1,7)}

def normalized_labels(y):
    out=[]
    for x in y:
        if x is None or str(x).lower() in ('nan','none','na','unlabeled',''):
            out.append('UNLABELLED')
        else:
            out.append(LABEL_MAP.get(str(x),str(x)))
    return np.asarray(out)

def set_seed(seed):
    import torch
    random.seed(seed);np.random.seed(seed);torch.manual_seed(seed)
    if torch.cuda.is_available():torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic=True
    torch.backends.cudnn.benchmark=False

def load_slice(h5ad_path,hvg=2000,n_pca=50):
    import scanpy as sc
    adata=sc.read_h5ad(str(h5ad_path))
    # IMPORTANT: do NOT discard unannotated cells/spots using ground truth.
    sc.pp.filter_genes(adata,min_cells=10)
    sc.pp.filter_cells(adata,min_genes=100)
    sc.pp.highly_variable_genes(adata,n_top_genes=hvg,flavor='seurat_v3')
    adata=adata[:,adata.var.highly_variable].copy()
    sc.pp.normalize_total(adata,target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata)
    sc.tl.pca(adata,n_comps=n_pca,random_state=0)
    X=adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
    pca=np.asarray(adata.obsm['X_pca'],dtype=np.float32)
    coords=np.asarray(adata.obsm['spatial'],dtype=np.float64)
    truth=np.asarray(adata.obs['layer_guess_reordered'] if 'layer_guess_reordered' in adata.obs else np.repeat('UNLABELLED',adata.n_obs))
    out={'X':X.astype(np.float32),'pca':pca,'coords':coords,'labels':normalized_labels(truth),
         'barcodes':np.asarray(adata.obs_names,dtype=str),'genes':np.asarray(adata.var_names,dtype=str)}
    del adata;gc.collect()
    return out

def sparse_historical_spatial_graph(coords,features,k=10,rank_block=500):
    """Replicates historical graph edge rule without forming dense NxN matrices.

    A is the symmetric matrix 0.5*(directed_weighted + transpose), where
    expression cosine is *Spearman* for directed edges with row<col only.
    This unusual triangular rule is inherited for a controlled comparison.
    """
    D=kneighbors_graph(coords,n_neighbors=k,mode='distance',include_self=False).tocsr()
    upper=sp.triu(D,k=1).tocoo()
    r=upper.row.astype(int);c=upper.col.astype(int);dist=upper.data
    positive=D.data[D.data>0]
    sigma=float(np.median(positive)) if len(positive) else 1.0
    sigma=max(sigma,1e-9)
    reverse=np.asarray(D[c,r]).ravel() if len(r) else np.array([])
    X=np.asarray(features,dtype=np.float64)
    ranks=rankdata(X,axis=1)
    ranks-=ranks.mean(axis=1,keepdims=True)
    norms=np.linalg.norm(ranks,axis=1)
    sim=np.zeros(len(r),dtype=np.float64)
    for first in range(0,len(r),rank_block):
        last=min(first+rank_block,len(r))
        u=r[first:last];v=c[first:last]
        numerator=np.einsum('ij,ij->i',ranks[u],ranks[v],optimize=True)
        denom=norms[u]*norms[v]
        cor=np.divide(numerator,denom,out=np.zeros_like(numerator),where=denom>0)
        sim[first:last]=np.clip(cor,-1,1)*0.5+0.5
    # For rows i<j: include i->j plus optional j->i if directed.
    a=.5*sim*(np.exp(-dist**2/(2*sigma**2)) + np.where(reverse>0,np.exp(-reverse**2/(2*sigma**2)),0.0))
    n=len(coords)
    A=sp.csr_matrix((np.r_[a,a],(np.r_[r,c],np.r_[c,r])),shape=(n,n))
    A.eliminate_zeros()
    degree=np.asarray(A.sum(axis=1)).ravel()
    m=float(degree.sum()*0.5)
    A_plus=A+sp.eye(n,format='csr')
    deg2=np.asarray(A_plus.sum(axis=1)).ravel()
    inv=np.divide(1,np.sqrt(deg2),out=np.zeros_like(deg2),where=deg2>0)
    A_norm=(sp.diags(inv)@A_plus@sp.diags(inv)).tocsr()
    return A,A_norm,degree,m

def to_sparse_tensor(A,device):
    import torch
    C=A.tocoo()
    ij=np.stack([C.row,C.col])
    return torch.sparse_coo_tensor(torch.from_numpy(ij).long(),torch.from_numpy(C.data.astype(np.float32)),size=C.shape,device=device).coalesce()

def train_local_graphsage(pca,coords,features,seed=42,epochs=500,device=None,verbose=False):
    import torch
    from torch import nn
    import torch.nn.functional as F
    set_seed(seed)  # BEFORE initializing either encoder.
    device=device or ('cuda' if torch.cuda.is_available() else 'cpu')
    A,An,deg,m=sparse_historical_spatial_graph(coords,features,k=10)
    A_t=to_sparse_tensor(A,device);An_t=to_sparse_tensor(An,device)
    deg_t=torch.tensor(deg,dtype=torch.float32,device=device)
    X=torch.tensor(pca,dtype=torch.float32,device=device)
    class Encoder(nn.Module):
        def __init__(self,input_dim):
            super().__init__()
            self.one=nn.Linear(input_dim*2,64,bias=False)
            self.two=nn.Linear(64*2,32,bias=False)
        def forward(self,X):
            h=F.leaky_relu(self.one(torch.cat([X,torch.sparse.mm(An_t,X)],dim=1)),0.2)
            return self.two(torch.cat([h,torch.sparse.mm(An_t,h)],dim=1))
    encoder=Encoder(X.shape[1]).to(device)
    mu=nn.Parameter(torch.empty(7,32,device=device))
    nn.init.xavier_uniform_(mu)
    optim=torch.optim.Adam(list(encoder.parameters())+[mu],lr=1e-3,weight_decay=1e-4)
    sched=torch.optim.lr_scheduler.CosineAnnealingLR(optim,T_max=epochs)
    history=[]
    for step in range(epochs):
        optim.zero_grad()
        H=encoder(X)
        R=F.softmax(F.normalize(H,dim=1)@F.normalize(mu,dim=1).T/0.18,dim=1)
        AR=torch.sparse.mm(A_t,R)
        # Same mathematical within-section modularity objective as dense B.
        trace=(R*AR).sum() - ((deg_t@R)**2).sum()/(2*m+1e-12)
        domain=-(trace/(m+1e-12))
        entropy=-(R*torch.log(R+1e-8)).sum(dim=1).mean()
        conf=R.max(dim=1).values.detach()
        mask=(conf>0.7).float().unsqueeze(1)
        smooth=R*mask
        lap=(torch.sparse.mm(A_t,smooth)-smooth).square().sum()
        # No anchor, no inter-slice community contrast, no cross-centroid term.
        loss=0.3*domain+0.03*entropy+0.01*lap
        loss.backward()
        nn.utils.clip_grad_norm_(list(encoder.parameters())+[mu],max_norm=1.0)
        optim.step();sched.step()
        if (step+1)%100==0 or step==epochs-1:
            record={'epoch':step+1,'loss':float(loss.item()),'domain':float(domain.item()),'entropy':float(entropy.item()),'laplacian':float(lap.item())}
            history.append(record)
            if verbose:print(record,flush=True)
    encoder.eval()
    with torch.no_grad():H=encoder(X).detach().cpu().numpy().astype(np.float32)
    del A_t,An_t,deg_t,X,encoder,mu,optim,sched
    gc.collect()
    if device=='cuda':torch.cuda.empty_cache()
    return H,history

def zca(H,eps=1e-8):
    centered=np.asarray(H,dtype=np.float64)-np.mean(H,axis=0)
    cov=(centered.T@centered)/max(1,len(H)-1)
    w,U=np.linalg.eigh(cov)
    return ((centered@U)*np.reciprocal(np.sqrt(np.maximum(w,0)+eps)))@U.T

def covariance_distance(X,Y):
    C1=np.cov(X,rowvar=False);C2=np.cov(Y,rowvar=False)
    return float(np.linalg.norm(C1-C2,'fro')/max(np.linalg.norm(C1,'fro'),np.linalg.norm(C2,'fro'),1e-12))

def axis_from_source(coords):
    from scipy.spatial import ConvexHull
    vertices=coords[ConvexHull(coords).vertices]
    center=coords.mean(axis=0)
    v=vertices[np.argmax(np.linalg.norm(vertices-center,axis=1))]-center
    return v/(np.linalg.norm(v)+1e-8)

def auto_groups(cs,ct,axis,K=7,scheme='pooled_quantile'):
    a=np.asarray(cs)@axis;b=np.asarray(ct)@axis
    # Normalize per section first, so differing coordinate ranges do not
    # by themselves assign all source spots to earlier groups.
    an=(a-a.min())/(np.ptp(a)+1e-8);bn=(b-b.min())/(np.ptp(b)+1e-8)
    if scheme=='pooled_quantile':edges=np.quantile(np.r_[an,bn],np.arange(1,K)/K)
    elif scheme=='equal_width':edges=np.linspace(0,1,K+1)[1:-1]
    else:raise ValueError(scheme)
    return np.searchsorted(edges,an,side='right'),np.searchsorted(edges,bn,side='right'),edges

def canonical_rotation(Cs,Ct):
    M=Cs.T@Ct
    U,s,Vt=np.linalg.svd(M,full_matrices=True)
    r=int(np.sum(s>1e-6*s[0])) if s[0]>0 else 0
    if r==len(U):return U@Vt,r
    U0=U[:,r:];V0t=Vt[r:,:]
    A,_,Bt=np.linalg.svd(V0t@U0,full_matrices=False)
    W=Bt.T@A.T
    return U[:,:r]@Vt[:r,:]+U0@W@V0t,r

def harmonize(Hs,Ht,gs,gt,whitening='conditional',tau=0.1,min_group_spots=20,completion='svd_default'):
    d=covariance_distance(Hs,Ht)
    active=whitening=='always' or (whitening=='conditional' and d>tau)
    S=zca(Hs) if active else Hs.astype(float)-Hs.mean(0)
    T=zca(Ht) if active else Ht.astype(float)-Ht.mean(0)
    common=[int(j) for j in sorted(set(gs)&set(gt)) if (gs==j).sum()>=min_group_spots and (gt==j).sum()>=min_group_spots]
    if len(common)<2:raise ValueError('insufficient common groups')
    Cs=np.asarray([S[gs==j].mean(0) for j in common]);Ct=np.asarray([T[gt==j].mean(0) for j in common])
    if completion=='svd_default':
        U,s,Vt=np.linalg.svd(Cs.T@Ct,full_matrices=True)
        rotation=U@Vt
        rank=int(np.sum(s>1e-6*s[0])) if s[0]>0 else 0
    elif completion=='min_change_sensitivity_only':
        rotation,rank=canonical_rotation(Cs,Ct)
    else:raise ValueError(completion)
    return S@rotation,T,{'covariance_distance':d,'whitened':bool(active),'common_groups':common,'rank':rank,'free_dim':S.shape[1]-rank,'completion':completion}

def full_target_lt(Hs,Ht,labels_s,labels_t):
    mask_s=np.isin(labels_s,LAYERS);mask_t=np.isin(labels_t,LAYERS)
    if not mask_s.any() or not mask_t.any():return {'LT':None,'balanced_accuracy':None,'per_layer_recall':{},'n_evaluable':0}
    S=normalize(Hs[mask_s],axis=1);T=normalize(Ht[mask_t],axis=1)
    y_s=labels_s[mask_s];y_t=labels_t[mask_t]
    pred=[]
    for st in range(0,len(T),350):
        nearest=(T[st:st+350]@S.T).argmax(axis=1)
        pred.extend(y_s[nearest].tolist())
    pred=np.asarray(pred)
    corr=pred==y_t
    recalls={k:float(np.mean(corr[y_t==k])) for k in LAYERS if (y_t==k).any()}
    return {'LT':float(corr.mean()),'balanced_accuracy':float(np.mean(list(recalls.values()))),
            'per_layer_recall':recalls,'n_evaluable':int(len(y_t)), 'target_label_fraction':float(mask_t.mean())}

def within_slice_probe(H,labels,seed=42):
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import train_test_split
    mask=np.isin(labels,LAYERS)
    X=H[mask];y=labels[mask]
    counts={k:int((y==k).sum()) for k in np.unique(y)}
    if len(counts)<2 or min(counts.values())<2:
        return {'probe_accuracy':None,'reason':'insufficient labelled classes','class_counts':counts}
    idx_train,idx_test=train_test_split(np.arange(len(y)),test_size=0.2,stratify=y,random_state=seed)
    clf=LogisticRegression(C=1.0,max_iter=2000,random_state=seed)
    clf.fit(X[idx_train],y[idx_train])
    return {'probe_accuracy':float(clf.score(X[idx_test],y[idx_test])),
            'n_labelled':int(len(y)),'class_counts':counts}

def posthoc_evaluation(hs,ht,ds,dt,K=7,tau=0.1,min_group_spots=20):
    axis=axis_from_source(ds['coords'])
    gs,gt,edges=auto_groups(ds['coords'],dt['coords'],axis,K,'pooled_quantile')
    scores={'baseline':full_target_lt(hs,ht,ds['labels'],dt['labels'])}
    cases={
      'Auto_conditional':(gs,gt,'conditional'),
      'Auto_procrustes_only':(gs,gt,'never'),
      'Auto_always_ZCA':(gs,gt,'always')}
    ys=ds['labels'];yt=dt['labels']
    label_map_to_int={k:j for j,k in enumerate(LAYERS)}
    oracle_gs=np.asarray([label_map_to_int.get(y,-1) for y in ys]);oracle_gt=np.asarray([label_map_to_int.get(y,-1) for y in yt])
    cases['Oracle_conditional']=(oracle_gs,oracle_gt,'conditional')
    for name,(a,b,w) in cases.items():
        try:
            # For Oracle, exclude -1 groups (unknown labels) from centroid fitting.
            if name.startswith('Oracle'):
                xa=a.copy();xb=b.copy()
                # use sentinel trick so negative ID cannot be considered common
                common_pos=sorted(set(xa[xa>=0])&set(xb[xb>=0]))
                xa[xa<0]=-100;xb[xb<0]=-200
                a,b=xa,xb
            S,T,diagnostics=harmonize(hs,ht,a,b,w,tau,min_group_spots)
            scores[name]={'metrics':full_target_lt(S,T,ys,yt),'diagnostics':diagnostics}
        except Exception as exc:scores[name]={'error':str(exc)}
    scores['spatial_prior']={'K':K,'axis':axis.tolist(),'pooled_depth_quantile_edges':edges.tolist(),
                             'source_group_counts':np.bincount(gs,minlength=K).tolist(),
                             'target_group_counts':np.bincount(gt,minlength=K).tolist()}
    return scores

def local_seed(seed,section_id):
    # Stable and role-independent: same section has the same training seed
    # whether it is source or target for a different section pair.
    h=int(hashlib.sha256(str(section_id).encode()).hexdigest()[:8],16)%1000000
    return int(seed+h)

def cli(args):
    import torch
    root=Path(args.data_dir);out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    descriptions={}
    for sid in (args.source,args.target):
        path=root/f'dlpfc_{sid}.h5ad'
        if not path.exists():raise FileNotFoundError(str(path))
        print('Preprocessing single section',sid,flush=True)
        descriptions[sid]=load_slice(path,hvg=args.hvg,n_pca=50)
    ds=descriptions[args.source];dt=descriptions[args.target]
    versions={'python':platform.python_version(),'numpy':np.__version__,'torch':torch.__version__}
    for seed in args.seeds:
        for model in args.models:
            embeddings={}
            for sid in (args.source,args.target):
                d=descriptions[sid]
                # Do not share optimizer, encoder state, PCA or source/target losses.
                cachepath=out/f'{sid}_seed{seed}_{model}_independent.npz'
                cached=bool(args.cache_latents and cachepath.exists())
                if cached:
                    with np.load(cachepath,allow_pickle=False) as stored:
                        assert stored['barcodes'].tolist()==d['barcodes'].tolist(), 'Cache barcode mismatch'
                        assert stored['genes'].tolist()==d['genes'].tolist(), 'Cache selected gene mismatch'
                        assert str(stored['pca_sha256'])==hashlib.sha256(d['pca'].tobytes()).hexdigest(), 'Cache PCA preprocessing mismatch'
                        assert int(stored['epochs'])==(0 if model=='pca' else args.epochs), 'Cache epoch mismatch'
                        H=stored['H'].copy()
                    loss=[]
                    print('Reusing independent section cache',cachepath,flush=True)
                elif model=='pca':H=d['pca'][:,:32].copy();loss=[]
                else:
                    print('Training independently',sid,model,'base seed',seed,'local seed',local_seed(seed,sid),flush=True)
                    H,loss=train_local_graphsage(d['pca'],d['coords'],d['X'],seed=local_seed(seed,sid),epochs=args.epochs,verbose=args.verbose)
                embeddings[sid]=H
                if (args.save_latents or args.cache_latents) and not cached:
                    np.savez_compressed(cachepath,H=H,coords=d['coords'],labels=d['labels'],
                                        barcodes=d['barcodes'],genes=d['genes'],pca_sha256=hashlib.sha256(d['pca'].tobytes()).hexdigest(),n_pca=50,epochs=0 if model=='pca' else args.epochs)
                if not cached:
                    (out/f'{sid}_seed{seed}_{model}_loss.json').write_text(json.dumps(loss,indent=2))
            score=posthoc_evaluation(embeddings[args.source],embeddings[args.target],ds,dt,K=args.K,tau=args.tau)
            result={'protocol':'strictly_independent_no_MNN_no_cross_slice_objective_v1',
                'model':model,'source':args.source,'target':args.target,'seed':seed,'section_initialization_seeds':{sid:local_seed(seed,sid) for sid in (args.source,args.target)},'epochs':0 if model=='pca' else args.epochs,
                'preprocessing':'independent_HVG_and_independent_PCA_by_section',
                'cross_section_training':False,'ground_truth_use':'metrics_and_Oracle_only',
                'versions':versions,'within_section_probe':{sid:within_slice_probe(embeddings[sid],descriptions[sid]['labels'],seed) for sid in (args.source,args.target)},'scores':score}
            dest=out/f'pair_{args.source}_{args.target}_seed{seed}_{model}.json'
            dest.write_text(json.dumps(result,indent=2,ensure_ascii=False))
            print('WROTE',dest,flush=True)
    print('COMPLETE; no result fields are claimed until this script runs on source data',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--data-dir',default=os.environ.get('LGA_DATA_DIR', './data'))
    p.add_argument('--output',default='./BIB_PAPER1_INDEPENDENT_RESULTS')
    p.add_argument('--source',default='151507');p.add_argument('--target',default='151508')
    p.add_argument('--seeds',type=int,nargs='+',default=[42,123,456])
    p.add_argument('--models',nargs='+',choices=['pca','graphsage'],default=['pca','graphsage'])
    p.add_argument('--hvg',type=int,default=2000);p.add_argument('--epochs',type=int,default=500)
    p.add_argument('--K',type=int,default=7);p.add_argument('--tau',type=float,default=.1)
    p.add_argument('--save-latents',action='store_true');p.add_argument('--cache-latents',action='store_true');p.add_argument('--verbose',action='store_true')
    cli(p.parse_args())
