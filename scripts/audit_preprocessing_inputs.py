#!/usr/bin/env python3
"""Read-only QC of the 14 LGA Paper-1 research input H5ADs.

Purpose: independently check preprocessing *compatibility*, not train anything.
Uses only numpy and h5py, with bounded reads for low-memory environments.
Run on the machine with H5AD inputs; does not transmit raw rows/barcodes.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import h5py


def decode(a):
    return [v.decode('utf-8','replace') if isinstance(v, (bytes, np.bytes_)) else str(v) for v in a]


def digest_names(ds, batch=4096):
    """Stable hash sensitive to ordered (not merely set-wise) gene identifiers."""
    h=hashlib.sha256()
    for p in range(0, len(ds), batch):
        for item in decode(ds[p:p+batch]):
            b=item.encode('utf-8')
            h.update(len(b).to_bytes(4,'little'))
            h.update(b)
    return h.hexdigest()


def fractions(x, eps=1e-4):
    if not len(x):
        return {"sampled_nnz": 0, "negative": 0, "nonfinite": 0,
                "noninteger": 0, "zero": 0, "min_value": None, "max_value": None}
    x=np.asarray(x,dtype=np.float64)
    finite=np.isfinite(x)
    good=x[finite]
    return {"sampled_nnz": int(len(x)), "negative": int((good < -eps).sum()),
            "nonfinite": int((~finite).sum()), "noninteger": int((np.abs(good-np.rint(good)) > eps).sum()),
            "zero": int((np.abs(good)<eps).sum()), "min_value": float(np.min(good)) if len(good) else None,
            "max_value": float(np.max(good)) if len(good) else None}


def sample_sparse_values(data, windows=16, size=2048, full_scan=False):
    n=len(data)
    if n==0:return fractions(np.empty(0))
    if full_scan:
        acc={"sampled_nnz":0,"negative":0,"nonfinite":0,"noninteger":0,"zero":0,
             "min_value":None,"max_value":None}
        for start in range(0,n,131072):
            item=fractions(data[start:start+131072])
            for k in ["sampled_nnz","negative","nonfinite","noninteger","zero"]:
                acc[k]+=item[k]
            for k in ["min_value","max_value"]:
                if item[k] is not None:
                    acc[k]=item[k] if acc[k] is None else (min(acc[k],item[k]) if k=="min_value" else max(acc[k],item[k]))
        return acc
    width=min(size,n)
    starts=np.unique(np.linspace(0,max(0,n-width),num=windows,dtype=np.int64))
    chunks=[np.asarray(data[int(s):int(s)+width]) for s in starts]
    return fractions(np.concatenate(chunks))


def inspect(file_path, rec, full_scan=False):
    with h5py.File(file_path,'r') as f:
        X=f['X']; nrows=int(rec['n_spots']); ncols=int(rec['n_genes'])
        assert isinstance(X,h5py.Group) and all(k in X for k in ('data','indices','indptr')), f"{file_path}: not expected CSR"
        indptr=np.asarray(X['indptr'][:],dtype=np.int64)
        assert len(indptr)==nrows+1 and int(indptr[-1])==len(X['data']), f"{file_path}: CSR index mismatch"
        assert int(X.attrs['shape'][0])==nrows and int(X.attrs['shape'][1])==ncols, f"{file_path}: dimension mismatch"
        values=sample_sparse_values(X['data'],full_scan=full_scan)
        n=min(96,nrows)
        row_idx=np.unique(np.linspace(0,nrows-1,num=n,dtype=int))
        sums=[]; nonzeros=[]
        for i in row_idx:
            lo,hi=int(indptr[i]),int(indptr[i+1])
            sums.append(float(np.asarray(X['data'][lo:hi],dtype=np.float64).sum()))
            nonzeros.append(hi-lo)
        coords=np.asarray(f['obsm']['spatial'][:],dtype=np.float64)
        assert coords.shape==(nrows,2), f"{file_path}: spatial dimension mismatch"
        names=f['var'][f['var'].attrs['_index']]
        assert len(names)==ncols, f"{file_path}: var length mismatch"
        label=f['obs'].get('layer_guess_reordered')
        labels=None
        if rec['id'].isdigit():
            assert isinstance(label,h5py.Group) and 'codes' in label and 'categories' in label, f"{file_path}: missing layer labels"
            cats=decode(label['categories'][:]);codes=np.asarray(label['codes'][:],dtype=int)
            assert len(codes)==nrows, f"{file_path}: label/spot mismatch"
            labels={str(c):int((codes==i).sum()) for i,c in enumerate(cats) if int((codes==i).sum())>0}
            labels['UNLABELLED']=int((codes<0).sum())
        result={"id":rec['id'],"matrix_shape":[nrows,ncols],"x_stored_dtype":str(X['data'].dtype),
           "x_stored_nonzero_elements":int(len(X['data'])),"sampled_matrix_values":values,
           "sampled_row_sums":{"n":len(sums),"min":round(min(sums),4),"median":round(float(np.median(sums)),4),
                                   "max":round(max(sums),4)},
           "sampled_row_nnz":{"min":min(nonzeros),"median":round(float(np.median(nonzeros)),1),"max":max(nonzeros)},
           "all_gene_names_order_sha256":digest_names(names),
           "spatial":{"all_finite":bool(np.isfinite(coords).all()),"unique_rows":int(len(np.unique(coords,axis=0))),
                      "x_range":[float(coords[:,0].min()),float(coords[:,0].max())],
                      "y_range":[float(coords[:,1].min()),float(coords[:,1].max())]},
           "layer_counts_including_unlabelled":labels,
           "unused_extra_embeddings":{k:list(f['obsm'][k].shape) for k in f['obsm'] if k!='spatial'}}
        return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-root',required=True,type=Path,help='Root containing j/ and mouse_brain/')
    p.add_argument('--manifest',type=Path,default=Path(__file__).resolve().parents[1]/'data'/'paper1_input_manifest.json')
    p.add_argument('--output',type=Path,default=Path('paper1_preprocessing_qc.json'))
    p.add_argument('--full-matrix',action='store_true',help='Scan ALL nonzero X values in 128K chunks (more I/O; no GPU)')
    args=p.parse_args()
    manifest=json.loads(args.manifest.read_text(encoding='utf-8'))
    report={"protocol":"Paper-1 first-pass preprocessing QC","scope":"read-only, no model training",
            "matrix_qc":"ALL nonzero elements" if args.full_matrix else "deterministic evenly spaced sample of X nonzero elements",
            "files":[],"conclusions":{}}
    for rec in manifest['files']:
        path=args.data_root/rec['relative_input_path']
        if not path.is_file():
            raise SystemExit(f'Missing input: {path}')
        result=inspect(path,rec,args.full_matrix)
        report['files'].append(result)
        s=result['sampled_matrix_values'];f=s['noninteger']/max(1,s['sampled_nnz'])
        print(f"{rec['id']}: dtype={result['x_stored_dtype']}, tested={s['sampled_nnz']:,} X values, "
              f"noninteger={100*f:.3f}%, median row sum={result['sampled_row_sums']['median']}, "
              f"missing labels={(result['layer_counts_including_unlabelled'] or {}).get('UNLABELLED','n/a')}",flush=True)
    hum=report['files'][:12];mouse=report['files'][12:]
    gene={r['id']:r['all_gene_names_order_sha256'] for r in hum}
    report['conclusions']['human_all_same_ordered_gene_list']=len(set(gene.values()))==1
    report['conclusions']['human_gene_list_digest_by_section']=gene
    report['conclusions']['donor_1_and_2_have_same_upstream_conversion_as_donor_3']=False
    report['conclusions']['human_missing_layer_labels']=sum(x['layer_counts_including_unlabelled']['UNLABELLED'] for x in hum)
    report['conclusions']['input_matrices_confirmed_raw_counts']=False
    report['conclusions']['raw_counts_comment']='Nonnegative integral X values are compatible with raw UMI counts but cannot by themselves prove upstream provenance or expression equivalence.'
    report['conclusions']['mouse_group_labels_only_added_during_experiment']=all(x['layer_counts_including_unlabelled'] is None for x in mouse)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('WROTE',args.output)
    print('DLPFC identical ordered genes:',report['conclusions']['human_all_same_ordered_gene_list'])
    print('Human missing annotation:',report['conclusions']['human_missing_layer_labels'])
    if any(x['sampled_matrix_values']['negative'] or x['sampled_matrix_values']['nonfinite'] for x in hum):
        print('WARNING: human X has negative or nonfinite values; check upstream pipeline')
    if any(x['sampled_matrix_values']['noninteger'] for x in hum):
        print('WARNING: at least one human X has noninteger values; seurat_v3 requires untransformed raw counts.')

if __name__=='__main__':main()