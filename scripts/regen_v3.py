import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
S5='/media/budheswar-dehury/DATADRIVE0/Fe-S/step5_scan'
df=pd.read_csv(f'{S5}/all_hits_final_v3.tsv',sep='\t',header=None,
   names=['genome','protein','hmm','score','family','ctype','status','provenance'])
func=df[df['status']=='functional'].copy()
print("functional:",len(func),"families:",func['family'].nunique())
print("by tier:",func['provenance'].value_counts().to_dict())
mat=func.pivot_table(index='family',columns='genome',values='protein',aggfunc='count',fill_value=0).astype(int)
mat.to_csv(f'{S5}/fes_family_matrix_v3.tsv',sep='\t')
mat.to_csv(f'{S5}/csv_family_genome_matrix_v3.csv')
ngen=mat.shape[1]
ab=func.groupby('family').agg(n_proteins=('protein','size'),cluster_type=('ctype',lambda x:x.mode().iloc[0]),
   tier=('provenance',lambda x:x.mode().iloc[0])).reset_index().sort_values('n_proteins',ascending=False)
ab.to_csv(f'{S5}/csv_family_abundance_v3.csv',index=False)
binary=(mat>0).astype(int)
prev=binary.sum(axis=1)
pa=pd.DataFrame({'family':prev.index,'n_genomes_present':prev.values,'prevalence_pct':(prev.values/ngen*100).round(1)})
cc=int(np.ceil(0.95*ngen)); rc=int(np.ceil(0.10*ngen))
pa['category']=pd.cut(pa['n_genomes_present'],bins=[-1,rc-1,cc-1,ngen],labels=['rare','accessory','core'])
pa=pa.sort_values('n_genomes_present',ascending=False)
pa.to_csv(f'{S5}/csv_family_presence_absence_v3.csv',index=False)
binary.to_csv(f'{S5}/csv_family_binary_matrix_v3.csv')
# abundance bar (color by tier: core_validated vs extended_pfam)
counts=ab.set_index('family')['n_proteins'].sort_values()
tier=ab.set_index('family')['tier']
colors=['#08519c' if tier[f]=='core_validated' else '#d94801' for f in counts.index]
fig,axx=plt.subplots(figsize=(8,11))
axx.barh(counts.index,counts.values,color=colors)
for i,(f,v) in enumerate(counts.items()): axx.text(v+max(counts.values)*0.01,i,str(v),va='center',fontsize=8)
axx.set_xlabel('Number of proteins'); axx.set_title(f'Fe-S family abundance (v3: {len(func)} proteins, 44 families)')
from matplotlib.patches import Patch
axx.legend(handles=[Patch(color='#08519c',label='core (validated)'),Patch(color='#d94801',label='extended (Pfam)')],loc='lower right')
plt.tight_layout(); plt.savefig(f'{S5}/fig_family_abundance_v3.png',dpi=300,bbox_inches='tight')
plt.savefig(f'{S5}/fig_family_abundance_v3.pdf',bbox_inches='tight'); plt.close()
# heatmaps
g=sns.clustermap(np.log1p(mat),cmap='YlGnBu',figsize=(16,11),xticklabels=False,yticklabels=True,cbar_kws={'label':'log(1+count)'})
g.fig.suptitle('Fe-S family distribution (v3, 44 families)',y=1.01)
g.savefig(f'{S5}/fig_family_heatmap_counts_v3.png',dpi=300,bbox_inches='tight'); plt.close()
b2=binary.loc[prev.sort_values(ascending=False).index]
fig,ax=plt.subplots(figsize=(16,11))
sns.heatmap(b2,cmap=['#f0f0f0','#08519c'],cbar=False,ax=ax,xticklabels=False,yticklabels=True)
ax.set_xlabel(f'{ngen} genomes'); ax.set_title('Fe-S family presence/absence (v3)')
for i,(f,p) in enumerate(prev.sort_values(ascending=False).items()): ax.text(b2.shape[1]+2,i+0.5,f'{p}/{ngen}',va='center',fontsize=7)
plt.tight_layout(); plt.savefig(f'{S5}/fig_family_presence_absence_v3.png',dpi=300,bbox_inches='tight'); plt.close()
print("\n=== core/accessory/rare (v3) ===")
print(pa['category'].value_counts().to_dict())
print("\n=== NEW families prevalence ===")
for f in ['IspH_LytB','Dihydroxyacid_dehydratase_IlvD','Quinolinate_synthase_NadA','ThiC','CDGSH_NEET','SufU_assembly','NifU_assembly','FhuF_siderophore_reductase']:
    if f in prev.index: print(f"  {f}: {prev[f]}/{ngen} ({100*prev[f]//ngen}%)")
