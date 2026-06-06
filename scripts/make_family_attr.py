import csv, sys

# usage: python3 make_family_attr.py <cyto_nodes_export.csv> <output_attr.tsv>
nodes_csv = sys.argv[1]
out_path  = sys.argv[2]

S8='/media/budheswar-dehury/DATADRIVE0/Fe-S/step8_coordination'

# master map: genome|protein -> fine_family + coordination status (built once, reused for every family)
fine={}; status={}
with open(f'{S8}/coordination_refined.tsv') as f:
    for line in f:
        p=line.rstrip('\n').split('\t')
        if len(p)>=9:
            key=f"{p[0]}|{p[1]}"
            fine[key]=p[2]      # family (fine)
            status[key]=p[8]    # coordination status

out=open(out_path,'w')
out.write("shared name\tfine_family\tcoord_status\tn_cys\tpredicted_cluster\n")
# also load n_cys + cluster for richer coloring
ncys={}; clu={}
with open(f'{S8}/coordination_refined.tsv') as f:
    for line in f:
        p=line.rstrip('\n').split('\t')
        if len(p)>=9:
            key=f"{p[0]}|{p[1]}"; ncys[key]=p[3]; clu[key]=p[5]

matched=0; total=0
with open(nodes_csv, newline='') as f:
    for row in csv.DictReader(f):
        total+=1
        key=row['Description'].split()[0]
        zid=row['shared name']
        if key in fine:
            out.write(f"{zid}\t{fine[key]}\t{status[key]}\t{ncys[key]}\t{clu[key]}\n"); matched+=1
        else:
            out.write(f"{zid}\tNA\tNA\tNA\tNA\n")
out.close()
print(f"{nodes_csv}: {total} nodes, {matched} matched -> {out_path}")
