import os, glob
BASE='/media/budheswar-dehury/DATADRIVE0/Fe-S'
S5=f'{BASE}/step5_scan'
BAKTA=f'{BASE}/bakta_out'
fes = {}
with open(f'{S5}/all_hits_final_v2.tsv') as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p)>=7 and p[6]=='functional':
            fes.setdefault(p[0], {})[p[1]] = (p[4], p[5])
print(f"genomes: {len(fes)}; proteins: {sum(len(v) for v in fes.values())}")
out = open(f'{BASE}/step11_function/fes_annotated.tsv','w')
out.write("genome\tlocus\tfamily\tcluster_type\tgene\tproduct\tKEGG_KO\tEC\tCOG\tGO\n")
nw=0; nko=0; miss=0
for genome, loci in fes.items():
    tsv = f'{BAKTA}/{genome}/{genome}.tsv'
    if not os.path.exists(tsv):
        cand = [c for c in glob.glob(f'{BAKTA}/{genome}*/*.tsv') if 'inference' not in c and 'hypothetical' not in c]
        if not cand:
            miss+=1; continue
        tsv = cand[0]
    with open(tsv) as f:
        for line in f:
            if line.startswith('#') or line.startswith('Sequence'): continue
            c = line.rstrip('\n').split('\t')
            if len(c)<6: continue
            locus = c[5]
            if locus in loci:
                fam, ctype = loci[locus]
                gene = c[6] if len(c)>6 else ''
                product = c[7] if len(c)>7 else ''
                dbx = (c[8] if len(c)>8 else '').replace(' ','')
                parts = dbx.split(',')
                ko = ';'.join(x.split(':',1)[1] for x in parts if x.startswith('KEGG:'))
                ec = ';'.join(x.split(':',1)[1] for x in parts if x.startswith('EC:'))
                cog= ';'.join(x.split(':',1)[1] for x in parts if x.startswith('COG:'))
                go = ';'.join(x.split(':',1)[1] for x in parts if x.startswith('GO:'))
                out.write(f"{genome}\t{locus}\t{fam}\t{ctype}\t{gene}\t{product}\t{ko}\t{ec}\t{cog}\t{go}\n")
                nw+=1
                if ko: nko+=1
out.close()
print(f"wrote {nw} annotated; {nko} have KEGG KO ({100*nko//max(nw,1)}%); missing-tsv genomes: {miss}")
