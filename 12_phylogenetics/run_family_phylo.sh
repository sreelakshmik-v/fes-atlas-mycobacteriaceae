#!/bin/bash
BASE=/media/budheswar-dehury/DATADRIVE0/Fe-S
cd $BASE/step13_phylo
FAM="$1"
MINSEQ="${2:-30}"
python3 extract_family.py "$FAM"
seqkit rmdup -s $FAM.faa -o ${FAM}_uniq.faa 2>/dev/null
N=$(grep -c '^>' ${FAM}_uniq.faa)
echo "$FAM unique: $N"
if [ "$N" -lt "$MINSEQ" ]; then echo "$FAM: SKIP (<$MINSEQ)"; exit 0; fi
IN=${FAM}_uniq.faa
if [ "$N" -gt 500 ]; then
  mmseqs easy-cluster ${FAM}_uniq.faa ${FAM}_f tmp_$FAM --min-seq-id 0.70 -c 0.8 >/dev/null 2>&1
  IN=${FAM}_f_rep_seq.fasta
  echo "$FAM filtered: $(grep -c '^>' $IN) reps"
fi
mafft --localpair --maxiterate 1000 --thread 16 $IN > ${FAM}.aln 2>/dev/null
trimal -in ${FAM}.aln -out ${FAM}_trim.aln -gt 0.5 2>/dev/null
COLS=$(awk '/^>/{if(s)print length(s);s="";next}{s=s$0}END{print length(s)}' ${FAM}_trim.aln | sort -nu | tail -1)
echo "$FAM aligned: $(grep -c '^>' ${FAM}_trim.aln) seqs x $COLS cols"
iqtree2 -s ${FAM}_trim.aln -B 3000 -m MFP -nstop 300 -pers 0.3 -T AUTO --prefix ${FAM}_tree -redo >/dev/null 2>&1
[ -f ${FAM}_tree.treefile ] && echo "$FAM TREE DONE" || echo "$FAM tree FAILED"
