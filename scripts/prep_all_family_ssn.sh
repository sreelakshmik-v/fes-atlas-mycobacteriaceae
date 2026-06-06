#!/bin/bash
cd /media/budheswar-dehury/DATADRIVE0/Fe-S/step9_ssn
mkdir -p per_family_ssn
S8=/media/budheswar-dehury/DATADRIVE0/Fe-S/step8_coordination

# superfamily map (rolls up sub-families), excluding likely-nonFeS
awk -F'\t' '$9!="likely_NONFES_demote"{
  fine=$3; sup=fine
  if(fine ~ /^Radical_SAM/) sup="Radical_SAM"
  else if(fine ~ /^Rieske/) sup="Rieske"
  else if(fine ~ /^Ferredoxin_4/) sup="Ferredoxin_4Fe4S"
  else if(fine ~ /^Ferredoxin_2/) sup="Ferredoxin_2Fe2S"
  else if(fine ~ /^DNA_glycosylase/) sup="DNA_glycosylase"
  else if(fine ~ /^Nitrite/) sup="Nitrite_sulfite_reductase"
  print $1"|"$2"\t"sup
}' $S8/coordination_refined.tsv > superfamily_map.tsv

# get all superfamilies and their sizes
cut -f2 superfamily_map.tsv | sort | uniq -c | sort -rn > superfamily_sizes.txt
echo "=== all superfamily sizes ==="
cat superfamily_sizes.txt
echo ""
echo "=== building SSN inputs for families with >=50 proteins ==="
while read n fam; do
  if [ "$n" -ge 50 ]; then
    awk -F'\t' -v f="$fam" '$2==f{print $1}' superfamily_map.tsv > /tmp/fids.txt
    seqkit grep -f /tmp/fids.txt ssn_input_all.faa > per_family_ssn/${fam}.faa 2>/dev/null
    cd-hit -i per_family_ssn/${fam}.faa -o per_family_ssn/${fam}_r85.faa -c 0.85 -n 5 -M 3000 -T 8 >/dev/null 2>&1
    reps=$(grep -c '^>' per_family_ssn/${fam}_r85.faa)
    echo "$fam: $n proteins -> $reps reps (SSN-worthy)"
  else
    echo "$fam: $n proteins -> SKIP (too small for meaningful SSN; shown in full-repertoire SSN)"
  fi
done < superfamily_sizes.txt
echo "DONE"
