import json,csv,collections,numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'serif','font.serif':['Liberation Serif','DejaVu Serif'],'font.size':11,'axes.labelsize':12,'axes.titlesize':11,'xtick.labelsize':10,'ytick.labelsize':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white','axes.linewidth':1.0,'ps.fonttype':42,'pdf.fonttype':42})
root=Path(__file__).resolve().parents[2];dst=Path(__import__('os').environ.get('LGA_FIGURE_OUT',str(root/'figures_generated')));dst.mkdir(parents=True,exist_ok=True)
D=json.load(open(root/'results/figure_inputs/e1_preprocessing.json'))
A={x['n_hvg']:x for x in D if x['mode']=='independent' and x['n_pca']==50}
S={x['n_hvg']:x for x in D if x['mode']=='shared' and x['n_pca']==50}
# Figure 5A: unpaired covariance discrepancy, correctly defined; only paired 1000/2000/3000 conditions
ks=[1000,2000,3000]; vals_shared=[S[k]['cd'] for k in ks];vals_ind=[A[k]['cd'] for k in ks]
fig,ax=plt.subplots(figsize=(5.4,3.6));x=np.arange(len(ks));w=.32
ax.bar(x-w/2,vals_shared,width=w,color='#D6DADE',edgecolor='black',linewidth=.7,label='Shared HVGs')
ax.bar(x+w/2,vals_ind,width=w,color='#4C78A8',edgecolor='black',linewidth=.7,label='Section-specific HVGs')
ax.set_xticks(x,labels=['1,000','2,000','3,000']);ax.set_xlabel('Number of HVGs');ax.set_ylabel('Normalized covariance distance');ax.set_ylim(0,.76);ax.legend(frameon=False,loc='upper left',fontsize=10)
fig.tight_layout(pad=1.0);fig.savefig(dst/'fig5_a.eps',format='eps',dpi=300);fig.savefig(dst/'fig5_a.pdf',dpi=300);plt.close(fig)
# Figure 5B: independent encoder parameterization, 12 pairs, five seeds, seeded paired means
def _encoder_rows(path):
    data=json.load(open(path))
    if isinstance(data,list): return data
    if isinstance(data,dict) and isinstance(data.get('results'),list): return data['results']
    raise ValueError(f'Unexpected encoder-sharing input format: {path}')
raw=_encoder_rows(root/'results/figure_inputs/encoder_sharing.json')
G=collections.defaultdict(list)
for r in raw:G[(r['donor'],r['pair'])].append(r)
order=sorted(G,key=lambda k:(int(k[0][-1]),k[1]))
y=np.arange(len(order));share=[np.mean([t['shared']['CD'] for t in G[k]]) for k in order];ind=[np.mean([t['independent']['CD'] for t in G[k]]) for k in order]
fig,ax=plt.subplots(figsize=(5.4,3.6))
for yy,a,b in zip(y,share,ind):ax.plot([a,b],[yy,yy],color='#9BA5AE',lw=1.15,zorder=1)
ax.scatter(share,y,s=25,c='#737E87',edgecolor='black',linewidth=.35,label='Shared encoder',zorder=2)
ax.scatter(ind,y,s=26,c='#4C78A8',edgecolor='black',linewidth=.35,label='Separate encoders',zorder=3)
ax.invert_yaxis();ax.set_yticks(y,labels=[f'D{k[0][-1]} {({"151507-151508":"AB","151508-151509":"BC","151507-151509":"AC","151509-151510":"CD","151669-151670":"AB","151670-151671":"BC","151669-151671":"AC","151671-151672":"CD","151673-151674":"AB","151674-151675":"BC","151673-151675":"AC","151675-151676":"CD"})[k[1]]}' for k in order]);ax.set_xlabel('Covariance Frobenius distance (a.u.)');ax.set_xlim(-.01,max(ind)*1.16)
ax.legend(frameon=False,loc='upper right',fontsize=9)
fig.tight_layout(pad=.8);fig.savefig(dst/'fig5_b.eps',format='eps',dpi=300);fig.savefig(dst/'fig5_b.pdf',dpi=300);plt.close(fig)
# Figure 6: controlled rotation on independent PCA+GraphSAGE (single common protocol)
pr=list(csv.DictReader(open(root/'results/verified/FIG6_18_PAIR_MEANS.csv')))
fig,axes=plt.subplots(1,2,figsize=(10.2,3.7),sharey=True)
for ax,model,title in zip(axes,['graphsage','pca'],['GraphSAGE (independent)','PCA (section-wise)']):
 rows=[r for r in pr if r['model']==model]
 by={d:[r for r in rows if int(r['donor'])==d] for d in [1,2,3]}
 x=np.arange(3);w=.34
 vr=[np.mean([float(r['rotated_LT']) for r in by[d]]) for d in [1,2,3]]
 va=[np.mean([float(r['auto_LT']) for r in by[d]]) for d in [1,2,3]]
 er=[np.std([float(r['rotated_LT']) for r in by[d]],ddof=1) for d in [1,2,3]]
 ea=[np.std([float(r['auto_LT']) for r in by[d]],ddof=1) for d in [1,2,3]]
 ax.bar(x-w/2,vr,w,color='#CFD4D7',edgecolor='black',linewidth=.75,yerr=er,error_kw={'lw':1,'capsize':3,'capthick':1},label='After rotation')
 ax.bar(x+w/2,va,w,color='#4C78A8',edgecolor='black',linewidth=.75,yerr=ea,error_kw={'lw':1,'capsize':3,'capthick':1},label='After LGA-Auto')
 for i,d in enumerate([1,2,3]):
  for q,r in enumerate(by[d]):
   ax.scatter(i-w/2, float(r['rotated_LT']),s=12,c='#5F6468',edgecolor='white',linewidth=.2,zorder=4)
   ax.scatter(i+w/2, float(r['auto_LT']),s=12,c='#16466B',edgecolor='white',linewidth=.2,zorder=4)
 ax.set_title(title);ax.set_xticks(x,labels=['Donor 1','Donor 2','Donor 3']);ax.set_ylim(0,.78);ax.set_xlabel('DLPFC donor')
 axes[0].set_ylabel('Target-to-source label-transfer accuracy')
 axes[0].legend(frameon=False,loc='upper left',fontsize=10)
fig.tight_layout(pad=1.2,w_pad=1.8);fig.savefig(dst/'fig6.eps',format='eps',dpi=300);fig.savefig(dst/'fig6.pdf',dpi=300);plt.close(fig)
print('Fig5 E1 available',[(k,S[k]['cd'],A[k]['cd']) for k in ks]);print('Fig5B raw covariance Frobenius',np.mean(share),np.mean(ind),'pairs',len(G));print('Wrote figures',list(dst.iterdir()))
