import pandas as pd
S5='/media/budheswar-dehury/DATADRIVE0/Fe-S/step5_scan'

df = pd.read_csv(f'{S5}/all_hits_final.tsv', sep='\t', header=None,
                 names=['genome','protein','hmm','score','family','ctype','status'])
func = df[df['status']=='functional']

# 1. ABUNDANCE csv: family, n_proteins, cluster_type
ab = func.groupby('family').agg(n_proteins=('protein','size'),
                                cluster_type=('ctype', lambda x: x.mode().iloc[0])).reset_index()
ab = ab.sort_values('n_proteins', ascending=False)
ab.to_csv(f'{S5}/csv_family_abundance.csv', index=False)

# 2. MATRIX csv (family x genome counts) — already have fes_family_matrix_final.tsv; also write as csv
mat = pd.read_csv(f'{S5}/fes_family_matrix_final.tsv', sep='\t', index_col=0)
mat.to_csv(f'{S5}/csv_family_genome_matrix.csv')

# 3. PRESENCE/ABSENCE csv: family, n_genomes_present, prevalence_pct, category
binary=(mat>0).astype(int)
prev=binary.sum(axis=1)
pa=pd.DataFrame({'family':prev.index,'n_genomes_present':prev.values,
                 'prevalence_pct':(prev.values/277*100).round(1)})
pa['category']=pd.cut(pa['n_genomes_present'],bins=[-1,27,262,277],
                      labels=['rare','accessory','core'])
pa=pa.sort_values('n_genomes_present',ascending=False)
pa.to_csv(f'{S5}/csv_family_presence_absence.csv',index=False)

# 4. binary matrix csv (for heatmap plotting)
binary.to_csv(f'{S5}/csv_family_binary_matrix.csv')

print("WROTE 4 CSVs:")
print(" csv_family_abundance.csv       (family, n_proteins, cluster_type)")
print(" csv_family_genome_matrix.csv   (family x 277 genomes, counts)")
print(" csv_family_presence_absence.csv(family, n_genomes, %, category)")
print(" csv_family_binary_matrix.csv   (family x 277, 0/1)")
print("\nABUNDANCE preview:")
print(ab.to_string(index=False))
print("\nPRESENCE preview:")
print(pa.to_string(index=False))
