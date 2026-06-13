import sys
from collections import Counter
from math import log2
aln=sys.argv[1]
# read alignment
seqs={}; name=None
for line in open(aln):
    if line.startswith('>'): name=line[1:].split()[0]; seqs[name]=[]
    else: seqs[name].append(line.strip())
seqs={k:''.join(v) for k,v in seqs.items()}
mat=list(seqs.values())
L=len(mat[0]); N=len(mat)
print(f"# {aln}: {N} seqs x {L} cols")
print("col\tconservation\tmost_common_aa\tpct")
scores=[]
for i in range(L):
    col=[s[i] for s in mat if s[i]!='-']
    if not col: scores.append((i+1,0,'-',0)); continue
    cnt=Counter(col)
    # conservation = 1 - normalized Shannon entropy
    tot=len(col); ent=-sum((c/tot)*log2(c/tot) for c in cnt.values())
    maxent=log2(min(20,len(cnt))) if len(cnt)>1 else 1
    cons=1-(ent/log2(20)) if len(cnt)>1 else 1.0
    top,topn=cnt.most_common(1)[0]
    scores.append((i+1,round(cons,3),top,round(100*topn/tot)))
# write all + report most conserved
out=open(aln.replace('.aln','_conservation.tsv'),'w')
out.write("col\tconservation\ttop_aa\tpct_gapless\n")
for s in scores: out.write(f"{s[0]}\t{s[1]}\t{s[2]}\t{s[3]}\n")
out.close()
# highlight cysteines (Fe-S ligands) and their conservation
cys=[s for s in scores if s[2]=='C' and s[3]>=50]
print(f"\n# Conserved cysteine columns (Fe-S ligand candidates): {len(cys)}")
for s in sorted(cys,key=lambda x:-x[1])[:15]:
    print(f"  col {s[0]}: conservation={s[1]}, C in {s[3]}% of seqs")
# top conserved overall
print(f"\n# Top 10 most conserved positions:")
for s in sorted(scores,key=lambda x:-x[1])[:10]:
    print(f"  col {s[0]}: {s[2]} (cons={s[1]}, {s[3]}%)")
