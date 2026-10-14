#!/usr/bin/env python3
"""Read-only comparison of processed Donor-3 H5AD X against original 10x counts.

Requires original vendor filtered_feature_bc_matrix.h5 files to be present.
Does not assume that the external processed H5AD was created from these raw files.
Only checks numeric count equality on available, matched barcodes and genes.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import h5py
import numpy as np

IDS=('151673','151674','151675','151676')

def dec(seq):
    return [s.decode('utf-8','replace') if isinstance(s,(bytes,np.bytes_)) else str(s) for s in seq]

def compare(h5ad,official,section,n_spots=96):
    with h5py.File(h5ad,'r') as a,h5py.File(official,'r') as b:
        g=b['matrix']
        official_genes=dec(g['features']['name'][:]);local_genes=dec(a['var'][a['var'].attrs['_index']][:])
        if len(official_genes)!=len(local_genes):
            return {'status':'incompatible_gene_count','local_n_genes':len(local_genes),'official_n_genes':len(official_genes)}
        if official_genes==local_genes:
            official_to_local=np.arange(len(local_genes));gene_mapping='identical_order'
        elif len(set(official_genes))==len(official_genes) and len(set(local_genes))==len(local_genes) and set(official_genes)==set(local_genes):
            index={x:i for i,x in enumerate(local_genes)}
            official_to_local=np.array([index[x] for x in official_genes]);gene_mapping='unique_names_mapped'
        else:
            return {'status':'ambiguous_gene_mapping','local_n_genes':len(local_genes),'official_n_genes':len(official_genes),'identical_gene_order':False}
        barcode=dec(g['barcodes'][:]);official_index={x:i for i,x in enumerate(barcode)}
        localbarcodes=dec(a['obs'][a['obs'].attrs['_index']][:])
        suffix='-'+section
        matched=[(j,official_index[bc[:-len(suffix)]]) for j,bc in enumerate(localbarcodes)
                 if bc.endswith(suffix) and bc[:-len(suffix)] in official_index]
        if len(matched)!=len(localbarcodes):
            return {'status':'unmatched_local_barcodes','matched':len(matched),'local':len(localbarcodes)}
        if n_spots>0 and len(matched)>n_spots:
            matched=[matched[i] for i in np.unique(np.linspace(0,len(matched)-1,n_spots,dtype=int))]
        x=a['X'];li=x['indptr'];lv=x['data'];lx=x['indices'];oi=g['indptr'];ov=g['data'];ox=g['indices']
        if len(oi)!=len(barcode)+1 or len(li)!=len(localbarcodes)+1:raise RuntimeError('Sparse matrix pointer inconsistent')
        bad=0;diff_max=0.0;noninteger=0;first=None
        for ix,(lidx,oidx) in enumerate(matched):
            lstart,lend=int(li[lidx]),int(li[lidx+1]);ostart,oend=int(oi[oidx]),int(oi[oidx+1])
            l_dense=np.zeros(len(local_genes),dtype=np.float64)
            o_dense=np.zeros(len(local_genes),dtype=np.float64)
            lvals=np.asarray(lv[lstart:lend],dtype=np.float64)
            ovals=np.asarray(ov[ostart:oend],dtype=np.float64)
            l_dense[np.asarray(lx[lstart:lend],dtype=np.intp)]=lvals
            o_dense[official_to_local[np.asarray(ox[ostart:oend],dtype=np.intp)]]=ovals
            noninteger+=int(np.count_nonzero(np.abs(lvals-np.rint(lvals))>1e-4))
            mx=float(np.max(np.abs(l_dense-o_dense)))
            diff_max=max(diff_max,mx)
            if mx>1e-5:
                bad+=1
                if first is None:first={'sample_index_no_barcode':ix,'max_absolute_count_difference':mx,'local_umi_sum':float(np.sum(l_dense)),'official_umi_sum':float(np.sum(o_dense))}
        return {'status':'compared','local_barcodes_matched':len(localbarcodes),'official_spots':len(barcode),
                'local_spots_checked':len(matched),'exact_with_tolerance':len(matched)-bad,'different_spots':bad,
                'max_absolute_count_difference':diff_max,'sampled_noninteger_nnz':noninteger,
                'gene_mapping':gene_mapping,'first_difference_without_barcode':first,
                'interpretation':'Only checked count matrix X, not spatial coordinates or layer label provenance.'}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-root',type=Path,required=True)
    p.add_argument('--official-dir',type=Path,required=True,help='Directory containing <section>_filtered_feature_bc_matrix.h5')
    p.add_argument('--n-spots',type=int,default=96,help='Evenly sampled local spots per section; 0 to check all')
    p.add_argument('--output',type=Path,default=Path('donor3_raw_count_comparison.json'))
    args=p.parse_args()
    records=[]
    for sec in IDS:
        path=args.data_root/'j'/f'dlpfc_{sec}.h5ad'
        official=args.official_dir/f'{sec}_filtered_feature_bc_matrix.h5'
        if not path.is_file():
            records.append({'id':sec,'status':'processed_input_missing'});continue
        if not official.is_file():
            records.append({'id':sec,'status':'official_precursor_not_available_locally'});continue
        rec={'id':sec,**compare(path,official,sec,args.n_spots)}
        records.append(rec)
        print(sec,rec['status'],rec.get('exact_with_tolerance','N/A'),'/',rec.get('local_spots_checked','N/A'),'matched numeric rows',flush=True)
    output={'purpose':'Distinguish raw count X from transformed/re-normalized X in externally preprocessed H5ADs',
            'note':'This compares only matching barcodes/features against a local 10x matrix; external processed release provenance remains unidentified.',
            'results':records}
    args.output.write_text(json.dumps(output,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('WROTE',args.output)
if __name__=='__main__':main()