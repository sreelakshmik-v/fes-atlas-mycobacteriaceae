import os, glob
from collections import defaultdict, Counter
BASE='/media/budheswar-dehury/DATADRIVE0/Fe-S'

# ---- 1. load product names for the new proteins ----
new_rows=[]  # (genome, locus, scan_family, score)
for line in open(f'{BASE}/step11b_expand/genuinely_new.tsv'):
    g,l,fam,sc=line.rstrip('\n').split('\t')
    if fam=='RmlD_sub_bind': continue          # false family - drop
    new_rows.append((g,l,fam,sc))

need=defaultdict(set)
for g,l,fam,sc in new_rows: need[g].add(l)
prod={}
for g,loci in need.items():
    tsv=f'{BASE}/bakta_out/{g}/{g}.tsv'
    if not os.path.exists(tsv): continue
    for line in open(tsv):
        if line.startswith('#'): continue
        c=line.split('\t')
        if len(c)>7 and c[5] in loci: prod[(g,c[5])]=c[7]

# ---- 2. classify each new protein -> final family + cluster type ----
NEW_FAM={
 'ILVD_EDD_N':('Dihydroxyacid_dehydratase_IlvD','4Fe-4S'),
 'LYTB':('IspH_LytB','4Fe-4S'),
 'NadA':('Quinolinate_synthase_NadA','4Fe-4S'),
 'ThiC_Rad_SAM':('ThiC','4Fe-4S'),
 'zf-CDGSH':('CDGSH_NEET','2Fe-2S'),
 'FhuF_C':('FhuF_siderophore_reductase','2Fe-2S'),
 'PhnJ':('PhnJ_CP_lyase','4Fe-4S'),
 'AFOR_C':('Aldehyde_Fd_oxidoreductase','4Fe-4S'),
 'Mob_synth_C':('Molybdopterin_MoaA_like','4Fe-4S'),
}
AUGMENT={  # new-scan family -> existing family name (cluster type from existing)
 'Oxidored_q6':('NADH_dehydrogenase','4Fe-4S'),
 'NADH_4Fe-4S':('NADH_dehydrogenase','4Fe-4S'),
 'NIR_SIR':('Nitrite_sulfite_reductase','4Fe-4S'),
 'Ring_hydroxyl_A':('Rieske','2Fe-2S'),
}
out_rows=[]   # (genome, locus, hmm, score, family, ctype, status)
dropped=Counter(); added=Counter()
for g,l,fam,sc in new_rows:
    p=prod.get((g,l),'').lower()
    # drop pseudogenes
    if '(pseudo)' in p: dropped['pseudogene']+=1; continue
    # junk/pseudogene-only families
    if fam in ('UPF0004','Fer2_3','Fer4_8'): dropped[fam]+=1; continue
    # --- Aconitase split by product ---
    if fam=='Aconitase':
        if 'isopropylmalate' in p:
            final,ct='IPM_isomerase_LeuC','4Fe-4S'
        else:
            final,ct='Aconitase','4Fe-4S'
    # --- SDH_beta misclassification fix ---
    elif fam=='SDH_beta':
        if 'serine' in p: final,ct='Serine_dehydratase','4Fe-4S'
        else: final,ct='SDH_FRD','4Fe-4S'
    # --- EndIII: split Nth vs MutY/adenine ---
    elif fam=='EndIII_4Fe-2S':
        if 'adenine' in p or 'a/g' in p: final,ct='DNA_glycosylase_MutY','4Fe-4S'
        else: final,ct='DNA_glycosylase_Nth','4Fe-4S'
    # --- NifU: separate the 53 mislabeled Rieske ---
    elif fam in ('NifU','NifU_N'):
        if 'rieske' in p: final,ct='Rieske','2Fe-2S'
        elif 'sufu' in p: final,ct='SufU_assembly','assembly'
        else: final,ct='NifU_assembly','assembly'
    elif fam in NEW_FAM:
        final,ct=NEW_FAM[fam]
    elif fam in AUGMENT:
        final,ct=AUGMENT[fam]
    else:
        dropped['unmapped_'+fam]+=1; continue
    out_rows.append((g,l,f'PFAM_{fam}',sc,final,ct,'functional'))
    added[final]+=1

print("=== ADDED proteins by final family ===")
for fam,n in added.most_common(): print(f"{n:5d}  {fam}")
print(f"\nTOTAL ADDED: {sum(added.values())}")
print("\n=== DROPPED ===")
for k,n in dropped.most_common(): print(f"{n:5d}  {k}")

# ---- 3. merge with existing v2 -> v3 ----
existing=[l.rstrip('\n') for l in open(f'{BASE}/step5_scan/all_hits_final_v2.tsv')]
with open(f'{BASE}/step5_scan/all_hits_final_v3.tsv','w') as o:
    for line in existing: o.write(line+'\n')
    for g,l,hmm,sc,fam,ct,st in out_rows:
        o.write(f"{g}\t{l}\t{hmm}\t{sc}\t{fam}\t{ct}\t{st}\n")

# ---- 4. new totals ----
funct=sum(1 for r in existing[1:] if r.split('\t')[6:7]==['functional']) if len(existing)>1 else 0
# recount properly
fcount=0; famset=Counter()
for line in open(f'{BASE}/step5_scan/all_hits_final_v3.tsv'):
    p=line.rstrip('\n').split('\t')
    if len(p)>=7 and p[6]=='functional':
        fcount+=1; famset[p[4]]+=1
print(f"\n=== NEW REPERTOIRE (v3) ===")
print(f"total functional Fe-S proteins: {fcount}")
print(f"total families: {len(famset)}")
print(f"(was 24,480 / 34 families)")
