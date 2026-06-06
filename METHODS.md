# The Fe–S Cluster Protein Atlas of *Mycobacteriaceae*
## Complete Step-by-Step Methods and Results Log (Steps 1–10)

**Project:** Comprehensive identification, classification, coordination analysis, network analysis, and structural validation of iron–sulfur (Fe–S) cluster proteins across the family *Mycobacteriaceae* (genus *Mycobacterium* sensu lato).

**Research question:** *What is the complete set of Fe–S cluster-containing proteins in mycobacteria, and how do they contribute to metabolism, redox homeostasis, and virulence?* (Half 1 = identification/classification; Half 2 = functional contribution.)

**Compute environment:** local Linux server (28 cores, 15 GB RAM; NVIDIA RTX A400, 4 GB VRAM); heavy data on `/media/budheswar-dehury/DATADRIVE0/Fe-S`. Conda envs: `bakta`, `fes_tools` (mmseqs2, mafft, trimal, hmmer, seqkit, pandas, cd-hit, foldseek). TMalign at `/usr/bin/TMalign`.

**Methodological template:** adapted from the Alba-domain SSN study (Sci Rep s41598-024-79937-4): PSI-BLAST/HMM → SSN (EFI-EST) → HHpred → structure (TMalign) → phylogeny. Reference biology: Vallières et al. *Metallomics* 2024 (Fe–S review); Wehrspan et al. 2022 (AlphaFold Fe–S coordination geometry); MetalPredator.

---

## SCOPE DECISIONS (locked at project start)

1. **Taxonomic scope:** family *Mycobacteriaceae* (five genera per GTDB/LPSN: *Mycobacterium*, *Mycobacteroides*, *Mycolicibacterium*, *Mycolicibacter*, *Mycolicibacillus*).
2. **Genome set:** GTDB species-representative genomes (R232, 15 Apr 2026), one per 95% ANI species cluster — avoids *M. tuberculosis* oversampling.
3. **Annotation:** uniform re-annotation with Bakta (not as-deposited NCBI annotations).
4. **Definition of "Fe–S protein":** proteins that coordinate an Fe–S cluster ([2Fe-2S], [3Fe-4S], [4Fe-4S]) as cofactor — clients + carriers; SUF/ISC biogenesis kept and tagged (systems-level framing).
5. **Two-phase question:** identify completely (Steps 1–10), then map function (Step 11+).

---

# STEP 1 — Curated, metadata-rich seed dataset

**Goal:** a locked, provenance-tracked seed FASTA with per-protein metadata. (Previous attempt failed by clustering before locking seeds and by losing metadata — corrected here.)

### 1.1 — Source
- Curated Fe–S protein supplementary (Vallières-lineage), 9 sheets: 1 SUF-machinery wide matrix + 8 species client sheets — Mtb (103), Mab (136), Mav (107), Mmar (108), Msm (136), Mul (75), Mbv (70), Mle (50). Identification method column = "UniProt + AlphaFold" / Estellon et al. (structurally vetted).

### 1.2 — Extract + validate accession list
```bash
libreoffice --headless --convert-to csv --outdir . seeds.xlsx
awk -F',' 'NR>1 {print $2}' seeds.csv | sed 's/[[:space:]]//g' | grep -v '^$' > ids_raw.txt   # ~775
sort -u ids_raw.txt > ids_unique.txt
# UniProt accession format check:
grep -vE '^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$' ids_unique.txt
```

### 1.3 — Fetch sequences (UniProt ID Mapping API)
```bash
JOB=$(curl -s --request POST 'https://rest.uniprot.org/idmapping/run' \
  --form 'from=UniProtKB_AC-ID' --form 'to=UniProtKB' \
  --form "ids=$(paste -sd, ids_unique.txt)" | grep -o '"jobId":"[^"]*"' | cut -d'"' -f4)
curl -s "https://rest.uniprot.org/idmapping/uniprotkb/results/stream/$JOB?format=fasta&compressed=false" -o seeds_raw.fasta
grep -c '^>' seeds_raw.fasta
```

### 1.4 — Reconcile requested vs returned (QC step previously skipped)
```bash
grep '^>' seeds_raw.fasta | cut -d'|' -f2 | sort -u > ids_returned.txt
comm -23 ids_unique.txt ids_returned.txt > ids_missing.txt   # resolve each: merged/demerged/obsolete
```

### 1.5 — Family-coverage confirmation (biological check)
```bash
grep '^>' seeds_raw.fasta > headers.txt
for fam in ferredoxin Rieske "radical SAM" aconitase dehydratase NEET CDGSH glycosylase WhiB SufB SufD IscA NfuA; do
  printf "%-15s %d\n" "$fam" "$(grep -ic "$fam" headers.txt)"; done
grep -iE 'whib' headers.txt   # WhiB1-7 must be present (central to redox/virulence half)
```
**Result:** all major Fe–S families present; WhiB fully covered; aconitase present under "aconitate hydratase" (9); 7 suf/isc/nfu entries present (kept, tagged biogenesis). Composition biologically coherent (Rieske, sulfite reductase, radical-SAM RlmN ×5, WhiB1–7, glycosylases). A few "Uncharacterized/Putative" entries flagged for scrutiny at structural stage.

### 1.6 — Rename headers to stable IDs
```bash
grep '^>' seeds_raw.fasta | sed -E 's/^>(sp|tr)\|([A-Z0-9]+)\|([A-Z0-9_]+) (.*) OS=(.*) OX=([0-9]+).*GN=([^ ]+).*/\2\t\1\t\4\t\5\t\6\t\7/' > meta_from_headers.tsv
# (all 774 parsed cleanly, zero GN=-less failures)
seqkit replace -p '^\w+\|(\w+)\|.*$' -r '{kv}' -k header_map_clean.tsv seeds_raw.fasta > seeds_renamed.fasta
```

### 1.7 — Deduplicate at 100% identity only (NOT CD-HIT 90%)
```bash
seqkit rmdup -s seeds_renamed.fasta -o seeds_locked.fasta -D seq_duplicates.txt
seqkit stats seeds_locked.fasta
```
**Result:** 774 input → **717 unique sequences** (57 collapses = byte-identical Mtb/Mbv orthologs + a few Mul/Mmar pairs; recorded in `seq_duplicates.txt`, not lost).

### 1.8 — Lock + document
```bash
chmod 444 seeds_locked.fasta header_map_clean.tsv seq_duplicates.txt
# master_metadata.tsv = 774 rows, cols: uniprot_acc, sprot_trembl, protein_name, species, species_taxid, gene_name
```
**Deferred (documented in README, trigger = before Step 11):** `category` (client/biogenesis/carrier) and `source_sheet`.

**STEP 1 DELIVERABLES:** `seeds_locked.fasta` (717), `master_metadata.tsv` (774 rows), `seq_duplicates.txt`, `SEED_README.txt`, `SEED_STATS.txt`.

---

# STEP 2 — Genome collection + uniform annotation

### 2.1 — GTDB representative set
- GTDB Release R232 (15 Apr 2026). Filtered `f__Mycobacteriaceae` species representatives → genome accession list.

### 2.2 — Download genomes (NCBI `datasets`)
```bash
datasets download genome accession --inputfile got_accessions.txt --include genome
# nucleotide FASTAs for each representative
```

### 2.3 — Quality filter (CheckM2)
- Retained genomes with completeness ≥95%, contamination ≤5%.

### 2.4 — Uniform annotation (Bakta)
```bash
bakta --db .../bakta_db/db --genus Mycobacterium --threads 16 <genome>.fna --output <out>
# full Bakta DB v6.0 (32 GB)
```
**Result:** **277 genomes** uniformly annotated. Verified real proteome count = 277 (excluding `*.hypotheticals.faa`); CDS range 3,186–7,862, mean 5,490.

**Documented limitation:** *M. leprae* (GCF_000195855.1, TN strain) = 3,989 CDS but only 25 Bakta-pseudogenes vs ~1,300 curated pseudogenes — Bakta/Pyrodigal calls intact ORFs over degraded regions.

**STEP 2 DELIVERABLES:** `bakta_out/` (277 dirs), `proteomes/` (277 .faa), `master_metadata.tsv`.

---

# STEP 3 — Family-aware HMM construction

### 3.1 — Cluster seeds (mmseqs2)
```bash
mmseqs easy-cluster seeds_locked.fasta seeds_clust tmp --min-seq-id 0.30 -c 0.70
```
**Result:** 210 clusters. Size distribution: 80 singletons, 42 doubletons, 19 triplets, 17 quads (=158 "thin"); 52 clusters ≥5 members. Largest = 26 members (MMAR_2932 → Ferredoxin_4Fe4S/Fer4_13).

### 3.2 — Build custom HMMs for the 52 large clusters
```bash
# per cluster:
mafft --auto cluster.faa > cluster.aln
trimal -in cluster.aln -out cluster.trim -gt 0.2
hmmbuild cluster.hmm cluster.trim
```
- GA (gathering) thresholds calibrated by scanning each HMM back against the seed set (lowest GA = 64.7 bits).
**Result:** 52 custom seed-built HMMs.

### 3.3 — Cover the 158 thin clusters with Pfam HMMs
```bash
hmmfetch Pfam-A.hmm <family>   # for each
```
- 17 Pfam Fe–S families fetched: Fer2, Fer4, Fer4_2/6/7/8/9/10/12/17, Rieske, Rieske_2, Radical_SAM, Molybdop_Fe4S4, Whib, SUFBD_core, UDG.

### 3.4 — Orphan triage
- 6 un-clustered orphans (`orphans_flagged.tsv`): MAB_4421 → SDH (covered by SDH HMM); MSMEG_1081 → Rieske (covered); **BCG_2352, MMAR_0590, Rv0246, Rv0311 = candidate-novel** (4 orphans carried forward for structural adjudication in Step 10).

### 3.5 — Assemble final library
**Result:** `fes_final.hmm` = **69 HMMs** (52 custom + 17 Pfam), all GA-tagged. Smoke-test on Mtb proteome = 68 unique proteins recovered.

**STEP 3 DELIVERABLES:** `fes_final.hmm` (+ .h3*), `hmm_family_map.tsv` (69 lines: HMM→family→provisional cluster_type), `hmms/`, alignments, `orphans_flagged.tsv`.

---

# STEP 4 — Proteome scan

### 4.1 — Scan all 277 proteomes
```bash
for prot in proteomes/*.faa; do
  hmmsearch --cpu 16 --cut_ga \
    --tblout scan_tbl/$(basename $prot .faa).tbl \
    --domtblout scan_dom/$(basename $prot .faa).dom \
    fes_final.hmm $prot
done
```

### 4.2 — Best-hit resolution
- 44,260 raw rows → **25,996 unique Fe–S proteins** (one best HMM hit per protein).
**Result:** mean 93.8 Fe–S proteins/genome; range **30** (*M. leprae* + *M. lepromatosis*, leprosy bacilli, joint-lowest) to **141** (*M. moriokaense*, environmental, Rieske-rich — validates biology).

**STEP 4 DELIVERABLES:** `all_hits_besthit.tsv`, `scan_tbl/`, `scan_dom/`.

---

# STEP 5–7 — Family classification, false-positive demotion, pseudogene flag, LOCK

### 5.1 — Map HMM → family (35 genuine families)
- `hmm_family_map.tsv` resolved each HMM (locus-tag-named HMMs resolved via UniProt names + Pfam hmmscan).

### 5.2 — Demote non-Fe–S families (independent KEGG evidence)
- **ParA_ATPase** (655 hits, KEGG K03496) — not Fe–S → demoted.
- **PrpD** methylcitrate dehydratase (10 hits, K01720) — not Fe–S → demoted.
- DUF779 (6 hits, K09959, relates Rv0526) — kept as candidate.
```bash
# removed 665 (ParA + PrpD) -> all_hits_validated.tsv
```

### 5.3 — Pseudogene flag (from Bakta gff3)
- 3,666 pseudogene tags across genomes → only **19 pseudogene-derived Fe–S hits** dataset-wide (Bakta calls few pseudogenes).

### 5.4 — LOCK
**Result:** `all_hits_final.tsv` = **25,312 functional Fe–S proteins, 35 families**, matrix 35 × 277.

### 5.5 — Family abundance (top families)
| Family | n | Family | n |
|---|---|---|---|
| Ferredoxin_4Fe4S | 4,681 | WhiB | 2,200 |
| Rieske | 4,266 | Ferredoxin_2Fe2S | 2,159 |
| Molybdoenzyme | 2,683 | DNA_glycosylase | 832 |
| Radical_SAM | 2,511 | NADH_dehydrogenase | 722 |

### 5.6 — Core/accessory/rare (prevalence across 277)
- **14 CORE families** (≥263/277, i.e. ≥95%): DNA_glycosylase, Glutamate_synthase, Ferredoxin_4Fe4S, WhiB, SDH_FRD, Radical_SAM, SUF_biogenesis (all 277); Rieske, Ferredoxin_2Fe2S (276); Molybdoenzyme (275); NADH_dehydrogenase (274); FprB (273); IspG (272); Radical_SAM_MoaA (266).
- **19 accessory**, **4 rare** (CobG 19, FadH 12, Serine_dehydratase 10, DUF779 6).
- *M. leprae* = 28 functional Fe–S proteins.

**STEP 5–7 DELIVERABLES:** `all_hits_final.tsv`, `fes_family_matrix_final.tsv`, `genome_family_counts_final.tsv`, `pseudogene_locustags.tsv`, CSV exports, 4 figures (abundance bar, clustered heatmap, binary presence/absence).

---

# STEP 8 — Coordination analysis (sequence layer)

### 8.1 — Family-aware Cys/His motif scan
- Script `coordination_scan3.py`: per-family expected coordinating-Cys count reflecting REAL biology (3-Cys families set to 3, not 4): radical-SAM (3 Cys + SAM), ferredoxin 4Fe4S (3 Cys ± Asp 4th), aconitase/IPM-isomerase (3 Cys + open substrate site), Rieske (2 Cys + 2 His), glutamate synthase [3Fe-4S] (3 Cys); genuine 4-Cys families = WhiB, molybdoenzyme, glycosylase.
- Two motif tiers per family (primary canonical spacing + relaxed variant spacing).

### 8.2 — Confidence assignment
**Result (final, after correcting per-family Cys minima):**
- **high = 18,050 (71%)** | medium = 6,473 (26%) | **low = 789 (3%)**
- Motif: 14,215 canonical (primary) + 3,835 variant (relaxed) = 18,050 with recognizable clustered motif (~56% textbook motif).

### 8.3 — Cluster-type distribution
- **[4Fe-4S] = 15,476 (61%)**, **[2Fe-2S] = 7,476 (30%)**, multi 1,129, assembly 556, [3Fe-4S] 396, FAD 273, unknown 6.

### 8.4 — Coordination refinement (flag likely non-Fe–S)
- The 789 low-confidence proteins broke down by family: **DNA_glycosylase 549**, SerA 135, WhiB 80, then small tails.
- DNA_glycosylase split: 477 FeS_supported (≥3 Cys), 221 ambiguous (2 Cys), 134 likely-nonFeS (≤1 Cys).
- SerA: 135/136 likely-nonFeS (≤2 Cys); full-length 529-aa phosphoglycerate dehydrogenase, KEGG non-Fe–S.
**Result (`coordination_refined.tsv`):** 24,717 FeS_supported, 326 ambiguous_structural, 269 likely_NONFES_demote → structural-target list = 601 (326 + 269 + 6 DUF779).

### 8.5 — KEY METHODOLOGICAL CAVEAT (recorded)
Cysteine COUNT is a screen, not proof of Fe–S binding (binding depends on 3D Sγ–Sγ distance). "Too few Cys" = sound demotion floor; "enough Cys" = only coordination-CONSISTENT, not confirmed. The FeS_supported set is "sequence/coordination-supported," not "confirmed." Structural validation (Step 10) is the arbiter.

**STEP 8 DELIVERABLES:** `coordination_per_protein.tsv`, `coordination_refined.tsv`, `structural_targets.tsv`, `coordination_scan3.py`.

---

# STEP 9 — Sequence Similarity Networks (SSN) — Alba-paper style

### 9.1 — Full-repertoire SSN input
```bash
# exclude likely-nonFeS (269); 25,312 - 269 = 25,043 proteins
awk -F'\t' '$9!="likely_NONFES_demote"{print $1"|"$2"\t"$3}' coordination_refined.tsv > ssn_protein_family.tsv
cd-hit -i ssn_input_all.faa -o ssn_rep70.faa -c 0.70 -n 5 -M 4000 -T 8
```
**Result:** 25,043 → **2,053 representatives** at 70% identity (CD-HIT 55% gave 980; 70% chosen for richer node count, closer to Alba's ~6,000).

### 9.2 — EFI-EST + Cytoscape
- Uploaded `ssn_rep70.faa` to EFI-EST (FASTA option, headers preserved). EFI-EST assigns `zzz####` node IDs; original `genome|locus_tag` preserved in Description column.
- Alignment score threshold = 4 (≈ Alba paper). Downloaded XGMML → Cytoscape → Prefuse Force-Directed layout.
- Node attribute (family + super_family) built by parsing Description (`make_family_attr.py`); imported keying on `shared name`.

### 9.3 — Result
- The 35 families separate cleanly into distinct clusters = **independent validation of HMM classification** (sequence clustering recovers functional families).
- Super-family node counts: Rieske 494, Ferredoxin_4Fe4S 474, Ferredoxin_2Fe2S 374, WhiB 286, Molybdoenzyme 169, Radical_SAM 115, DNA_glycosylase 96; small families = peripheral specks (expected).

### 9.4 — Per-family SSNs
- Reusable prep (`prep_all_family_ssn.sh`): families with ≥50 representatives at 85% are SSN-worthy. **Worth submitting:** Rieske (1,560 reps), Ferredoxin_4Fe4S (1,321), Ferredoxin_2Fe2S (1,076), Molybdoenzyme (827), WhiB (595), Radical_SAM (585), DNA_glycosylase (388), FprB (101), CCG (83), NADH (68), SUF (54).
- **Too conserved for SSN** (few reps despite many proteins — under purifying selection): SDH_FRD (9), Glutamate_synthase (21), Nitrite_reductase (15), Aconitase (3), IspG (2), IPM_isomerase (2). [Biologically meaningful: conserved-essential vs diverse-accessory split.]
- **Rieske per-family SSN done:** sub-types (steroid_oxygenase, nitrite_reductase, QcrA) separate spatially from core Rieske = functional diversification visible. Disconnection between sub-clusters = the signal (distinct sub-families), not a bug.
- `make_family_attr.py` is reusable per family (4 colorable attributes: fine_family, coord_status, n_cys, predicted_cluster). Remaining per-family SSNs to be batched.

**STEP 9 DELIVERABLES:** `ssn_rep70.faa`, `rep70_family_attr2.tsv`, `cyto_*` files, `superfamily_map.tsv`, `per_family_ssn/`, `make_family_attr.py`, full-repertoire + Rieske SSN figures.

---

# STEP 10 — Structural validation (Depth A: validate + novelty)

**Tools:** HHpred (MPI Bioinformatics Toolkit, fold ID) | AlphaFold3 (structures; local GPU too small at 4 GB) | TMalign (local, fold superposition) | AlphaFill (cofactor transplant) | foldseek | BioPython MMCIFParser (Sγ–Sγ geometry; [4Fe-4S]~6.3 Å, [2Fe-2S]~3.5 Å per Wehrspan).

## 10.1 — DUF779 (strongest novelty result) → IscA/HesB Fe–S assembly fold
- 6 proteins, one homogeneous family (158–166 aa; 2 mis-labeled "Acetaldehyde dehydrogenase" by Bakta). One HHpred representative suffices (`LNBHIF_00484`, H37Rv ortholog, relates Rv0526; seed = BCG_0499/A0A0H3M2A1).
- **HHpred:** 98.8% to HesB/IscA Fe–S-biogenesis fold — PF05610 (DUF779) 100%, then 2P2E "Putative Fe-S biosynthesis protein", SCOP HesB-like (d1nwba, d2apna), IscA/YfhF (d1r94a), ECOD Fe-S_biosyn, PF01521 "Iron-sulphur cluster biosynthesis" — all ~98.8%.
- **TMalign (AF3 model vs references):** vs IscA 1R94 = TM **0.64**; vs 2P2E = TM **0.56** (both >0.5 = same fold; ~12–16% seq id; DUF779 = IscA core + ~70 extra residues).
- **Cys geometry (AF3 monomer):** 5 Cys, only one close pair (Cys35–Cys36 = 3.8 Å); rest 17–38 Å. Consistent with IscA-family **inter-subunit (dimer-interface)** coordination — monomer cannot show the full cage.
- **AlphaFill:** no cofactor transplanted (empty table; IscA structures usually apo) — binding cannot be confirmed via this route.
- **VERDICT:** IscA/HesB-like fold confirmed ×2 (HHpred + TMalign); putative Fe–S assembly/scaffold protein. Binding inferred, not proven. KEPT.
- **CAVEAT:** known fold (not a novel fold). **Literature check pending** (DUF779↔Fe–S may already be published) before claiming novelty.

## 10.2 — DNA_glycosylase reclassification (MAJOR correction)
**Trigger:** coordination flagged 549/837 DNA_glycosylase as low-Cys; user questioned validity of Cys-count filtering.
- Defining HMM = **Pfam UDG** (function-based model = uracil-DNA glycosylase, non-Fe–S by definition).
- **HHpred on 4 reps** (2, 2, 9, 10 Cys — controls BPHPIM_02234/KMCKDB_02864 + ambiguous ENJDEI_05570/HNDBBG_03927): ALL hit **UDG fold (SCOP c.18.1)**, non-Fe–S templates (1UDG, 1UI0, 1VK2, 4WPK). Even 9–10 Cys reps were non-Fe–S (scattered cysteines, no FCL) — **directly vindicates the Cys-count concern.**
- **FCL motif scan** on the 283 ≥4-Cys members (regex `C.{3,8}C.{1,4}C.{2,8}C`): **0 of 283** have the iron-sulfur cluster loop motif.
- **MutY + Nth reps** (IHOGBB_05096, CGIDCJ_00434; both 10 Cys, C-terminal FCL): HHpred → **HhH-GPD fold (SCOP a.96)**, EndoIII/MutY templates **with bound [4Fe-4S]** (1ORN, 2ABK, 3N5N "HET: SF4").
- **THREE converging lines** → broad `DNA_glycosylase` (UDG, 832 functional) is **non-Fe–S uracil glycosylase → DEMOTED**; genuine Fe–S glycosylases = **MutY (187) + Nth (172)** (HhH-GPD [4Fe-4S]), retained.
- **Note on rigor:** initially over-reached toward demoting on 4 examples; corrected by requiring the HMM-identity + FCL-scan evidence before family-level demotion.

## 10.3 — Orphan adjudication (4 candidate-novel singletons)
| Orphan | aa | HHpred fold | Verdict |
|---|---|---|---|
| BCG_2352 | 128 | Molybdoenzyme/formate-dehydrogenase [4Fe-4S] (YjgC; 2IV2 FE4S4) | **Genuine Fe–S** |
| MMAR_0590 | 102 | UDG fold; 5 Cys, **no FCL** | non-Fe–S |
| Rv0246 | 436 | MFS membrane transporter (100%) | non-Fe–S |
| Rv0311 | 409 | DUF6125 / TRAPP / H-NOX | non-Fe–S |
- **1 of 4 confirmed Fe–S; 3 reclassified.** Honest outcome — flagging was appropriately cautious.

## 10.4 — SerA_dehydrogenase
- 135/136 ≤2 Cys; full-length 529-aa phosphoglycerate dehydrogenase; KEGG non-Fe–S enzyme. Flagged likely non-Fe–S (count + length + annotation conclusive; structural HHpred optional).

## 10.5 — Repertoire correction propagated
```bash
# demote UDG: all_hits_final.tsv -> all_hits_final_v2.tsv (status nonFeS_UDG_demoted)
# 837 DNA_glycosylase rows = 832 functional + 5 pseudogene
# 25,312 - 832 = 24,480 functional
# regenerate matrix, 4 CSVs, 3 figures (all *_corrected)
```
**CORRECTED REPERTOIRE:** **24,480 functional Fe–S proteins | 34 families | 13 core** (DNA_glycosylase UDG dropped from core; MutY 67.5% + Nth 62.1% are accessory).
- 13 core families: Glutamate_synthase, Ferredoxin_4Fe4S, WhiB, SDH_FRD, Radical_SAM, SUF_biogenesis (277); Rieske, Ferredoxin_2Fe2S (276); Molybdoenzyme (275); NADH_dehydrogenase (274); FprB (273); IspG (272); Radical_SAM_MoaA (266).

## 10.6 — Validation scoping (why not all 35 families)
- Families built from curated UniProt Fe–S seeds (≈25 families) are literature-grounded; definitional-Fe–S Pfam families (Fer2/Fer4/Rieske/Radical_SAM/Molybdop_Fe4S4/Whib) have Fe–S binding intrinsic to the domain. **UDG was the unique risk:** the only family defined by a function-based Pfam model (not a cofactor-based one), which over-recruited non-Fe–S members. Risk-targeted validation, not blanket.

**STEP 10 DELIVERABLES:** `step10_structural/` (HHpred .hhr files, AF3 .cif/.pdb models, TMalign outputs, reference PDBs 1R94/2P2E, `flagged_validation/`, `orphans4.faa`, FCL scans), `all_hits_final_v2.tsv`, `*_corrected` matrix/CSVs/figures, `STEP10_structural_validation_README.md`.

---

# CURRENT STATUS & REMAINING WORK

**DONE (Half 1 — identification/classification, complete & structurally validated):**
- Steps 1–10: seeds → genomes → annotation → HMMs → scan → classification → demotion → coordination → SSN → structural validation.
- Final repertoire: **24,480 Fe–S proteins, 34 families, 13 core**, across 277 genomes.

**POLISH (close Depth A clean):**
- DUF779 literature check (novelty framing).
- Propagate UDG flag into `coordination_refined.tsv` for full consistency.
- Optional positive-control HHpred (1 Fer4 + 1 Radical_SAM) + SerA HHpred.
- Recolor full-repertoire SSN by corrected status (corroborates UDG correction).

**REMAINING ANALYSES (roadmap):**
- **Step 11 — Functional mapping (Half 2 of research question):** KEGG/KO + GO (kofamscan) → metabolism (TCA, respiration, cofactor biosynthesis); redox homeostasis (ferredoxins, WhiB regulators, Rieske); virulence (WhiB3/WhiB6/WhiB7, ESX, SUF operon) + Mtb essentiality (DeJesus TnSeq). Requires completing deferred `category` metadata column.
- **Step 12 — Per-family SSNs:** batch the 11 qualifying families through EFI-EST.
- **Step 13 — Phylogenetics:** per-family mafft→trimal→IQ-TREE; gain/loss across 277; dN/dS.
- **Step 14 — Depth-B comprehensive structural analysis:** representative folds per family, TMalign all-vs-all clustering (Alba Fig 3 equivalent).
- **Step 15 — Synthesis / thesis writing.**

---

*Document generated 2026-06-03. All counts traceable to `all_hits_final_v2.tsv` and `coordination_refined.tsv`.*
