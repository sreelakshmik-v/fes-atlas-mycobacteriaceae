import re

# (cluster_type, min_cys_expected, primary_motif, relaxed_motif, note)
# min_cys_expected reflects REAL coordination biology (3-Cys families set to 3)
FAMILY_MODEL = {
 'Radical_SAM':            ('4Fe-4S',3, r'C.{3}C.{2}C', r'C.{1,6}C.{1,5}C','CxxxCxxC + SAM'),
 'Radical_SAM_MoaA':       ('4Fe-4S',3, r'C.{3}C.{2}C', r'C.{1,6}C.{1,5}C','radical-SAM'),
 'Radical_SAM_MiaB':       ('4Fe-4S',3, r'C.{3}C.{2}C', r'C.{1,6}C.{1,5}C','radical-SAM +2nd cluster'),
 'Radical_SAM_LipA':       ('4Fe-4S',3, r'C.{3}C.{2}C', r'C.{1,6}C.{1,5}C','radical-SAM +2nd cluster'),
 'Radical_SAM_BioB':       ('4Fe-4S',3, r'C.{3}C.{2}C', r'C.{1,6}C.{1,5}C','radical-SAM +2Fe2S'),
 'Radical_SAM_RlmN':       ('4Fe-4S',3, r'C.{3}C.{2}C', r'C.{1,6}C.{1,5}C','radical-SAM'),
 'Radical_SAM_HemN':       ('4Fe-4S',3, r'C.{3}C.{2}C', r'C.{1,6}C.{1,5}C','radical-SAM'),
 'Radical_SAM_PqqE':       ('4Fe-4S',3, r'C.{3}C.{2}C', r'C.{1,6}C.{1,5}C','radical-SAM'),
 'FbiC_FO_synthase':       ('4Fe-4S',3, r'C.{3}C.{2}C', r'C.{1,6}C.{1,5}C','FO synthase'),
 'Ferredoxin_4Fe4S':       ('4Fe-4S',3, r'C.{2}C.{2}C.{2,40}C', r'C.{1,15}C.{1,15}C','bacterial Fdx (3-4 Cys, may have Asp 4th)'),
 'Ferredoxin_2Fe2S':       ('2Fe-2S',4, r'C.{4}C.{2}C.{20,45}C', r'C.{1,12}C.{1,15}C.{10,60}C','plant-type 2Fe-2S Fdx'),
 'Rieske':                 ('2Fe-2S',2, r'C.H.{10,55}C.{2}H', r'C.H.{5,70}C.{1,4}H','2Cys+2His'),
 'Rieske_QcrA':            ('2Fe-2S',2, r'C.H.{10,55}C.{2}H', r'C.H.{5,70}C.{1,4}H','Rieske QcrA'),
 'Rieske_steroid_oxygenase':('2Fe-2S',2,r'C.H.{10,55}C.{2}H', r'C.H.{5,70}C.{1,4}H','Rieske steroid'),
 'Rieske_nitrite_reductase':('2Fe-2S',2,r'C.H.{10,55}C.{2}H', r'C.H.{5,70}C.{1,4}H','Rieske nitrite'),
 'WhiB':                   ('4Fe-4S',4, r'C.{1,3}C.{10,40}C.{1,3}C', r'C.{1,8}C.{5,50}C.{1,8}C','WhiB 4-Cys'),
 'SDH_FRD':                ('multi',4,  r'C.{2,5}C.{2,5}C', r'C.{1,12}C.{1,12}C','SDH/FRD 2Fe2S+3Fe4S+4Fe4S'),
 'NADH_dehydrogenase':     ('multi',4,  r'C.{2,5}C.{2,5}C', r'C.{1,12}C.{1,12}C','Complex I Fe-S'),
 'Molybdoenzyme':          ('4Fe-4S',4, r'C.{2}C.{2}C.{2,40}C', r'C.{1,12}C.{1,12}C.{1,60}C','molybdoenzyme'),
 'Aconitase':              ('4Fe-4S',3, r'C.{20,200}C.{5,40}C', r'C.{10,250}C.{2,60}C','3-Cys + open site'),
 'IPM_isomerase_LeuC':     ('4Fe-4S',3, r'C.{20,200}C.{5,40}C', r'C.{10,250}C.{2,60}C','3-Cys + open site'),
 'DNA_glycosylase':        ('4Fe-4S',4, r'C.{2,8}C.{2,8}C.{2,8}C', r'C.{1,15}C.{1,15}C.{1,15}C','glycosylase 4-Cys'),
 'DNA_glycosylase_MutY':   ('4Fe-4S',4, r'C.{2,8}C.{2,8}C.{2,8}C', r'C.{1,15}C.{1,15}C.{1,15}C','MutY 4-Cys'),
 'DNA_glycosylase_Nth':    ('4Fe-4S',4, r'C.{2,8}C.{2,8}C.{2,8}C', r'C.{1,15}C.{1,15}C.{1,15}C','Nth 4-Cys'),
 'SUF_biogenesis':         ('assembly',0,r'C.{1,60}C', r'C.{1,100}C','transient assembly cluster'),
 'Glutamate_synthase':     ('3Fe-4S',3, r'C.{2,20}C.{2,20}C', r'C.{1,40}C.{1,40}C','3Fe-4S, 3-Cys'),
 'Nitrite_sulfite_reductase':('4Fe-4S',4,r'C.{2,10}C.{30,80}C.{2,5}C', r'C.{1,20}C.{10,100}C.{1,12}C','siroheme-coupled 4Fe-4S'),
 'IspG_HMBPP_synthase':    ('4Fe-4S',3, r'C.{2,30}C.{2,30}C', r'C.{1,60}C.{1,60}C','IspG/GcpE 3-Cys'),
 'Ferredoxin_reductase_FprB':('FAD',0,  r'C.{1,200}C', r'C.{1,300}C','flavoenzyme, Fe-S via partner'),
 'FeS_reductase_CCG':      ('4Fe-4S',2, r'C.{1,10}C.{1,10}C', r'C.{1,30}C.{1,30}C','CCG cysteine-rich'),
 'SerA_dehydrogenase':     ('4Fe-4S',3, r'C.{2,40}C.{2,40}C', r'C.{1,80}C.{1,80}C','SerA'),
 'Serine_dehydratase':     ('4Fe-4S',3, r'C.{2,40}C.{2,40}C', r'C.{1,80}C.{1,80}C','serine dehydratase'),
 'CobG_precorrin':         ('4Fe-4S',4, r'C.{2,30}C.{2,30}C.{2,30}C', r'C.{1,60}C.{1,60}C.{1,60}C','CobG'),
 'FadH_dienoyl_CoA':       ('4Fe-4S',4, r'C.{2,30}C.{2,30}C.{2,30}C', r'C.{1,60}C.{1,60}C.{1,60}C','FadH'),
 'DUF779_FeS_candidate':   ('unknown',2,r'C.{1,40}C', r'C.{1,80}C','DUF779 candidate'),
}

def analyze(seq, family):
    seq = seq.upper().replace('*','')
    cys=[i+1 for i,a in enumerate(seq) if a=='C']
    his=[i+1 for i,a in enumerate(seq) if a=='H']
    n_cys=len(cys)
    model=FAMILY_MODEL.get(family)
    if model is None:
        return (n_cys,'NA','unassigned',f"Cys@{','.join(map(str,cys[:6]))}",'low')
    ctype,min_cys,primary,relaxed,note=model
    has_p=bool(re.search(primary,seq))
    has_r=bool(re.search(relaxed,seq))
    enough=(min_cys==0) or (n_cys>=min_cys)
    if has_p and enough:
        conf='high'; matched='primary'
    elif has_r and enough:
        conf='high'; matched='relaxed'
    elif enough and ctype in ('assembly','FAD'):
        conf='medium'; matched='special'
    elif enough:
        conf='medium'; matched='cys_count_ok'
    else:
        conf='low'; matched='insufficient_cys'   # genuine structural target
    coord=f"Cys@{','.join(map(str,cys[:6]))}"
    if 'Rieske' in family: coord+=f"|His@{','.join(map(str,his[:4]))}"
    return (n_cys,matched,ctype,coord,conf)

def iter_fasta(path):
    h,s=None,[]
    for line in open(path):
        line=line.rstrip()
        if line.startswith('>'):
            if h: yield h,''.join(s)
            h=line[1:];s=[]
        else: s.append(line)
    if h: yield h,''.join(s)

fam_of={}
for ln in open('fes_protein_list.tsv'):
    g,p,f=ln.rstrip('\n').split('\t'); fam_of[f"{g}|{p}"]=(g,p,f)

out=open('coordination_per_protein.tsv','w')
out.write("genome\tprotein\tfamily\tn_cys\tmotif_matched\tpredicted_cluster\tcoordinating_residues\tconfidence\n")
n=0
for hdr,seq in iter_fasta('fes_all_proteins.faa'):
    key=hdr.split()[0]
    if key not in fam_of: continue
    g,p,f=fam_of[key]
    nc,mm,ct,co,cf=analyze(seq,f)
    out.write(f"{g}\t{p}\t{f}\t{nc}\t{mm}\t{ct}\t{co}\t{cf}\n"); n+=1
out.close()
print("analyzed",n)
