# The Fe–S Cluster Protein Atlas of *Mycobacterium*
## Complete step-by-step methods log (Steps 1–15)

**Project:** Genome-wide identification, classification, coordination analysis, network analysis,
structural comparison, conservation, phylogenetics and essentiality mapping of iron–sulfur (Fe–S)
cluster proteins across the genus *Mycobacterium*.

**Manuscript:** Sreelakshmi K V, Sarode VV, Shenoy G G, Dehury B. *Unified by chemistry, not fold:
A genus-wide atlas of the mycobacterial iron–sulfur proteome and its antitubercular target space.*

**Research question:** What is the complete set of Fe–S cluster-containing proteins in
*Mycobacterium*, how are those proteins related structurally and evolutionarily, and what does the
repertoire imply for antitubercular target selection?

**Compute environment:** local Linux server (28 cores, 15 GB RAM; NVIDIA RTX A400, 4 GB VRAM); data
on `/media/budheswar-dehury/DATADRIVE0/Fe-S`. Conda envs: `bakta`, `fes_tools` (mmseqs2, mafft,
trimal, hmmer, seqkit, pandas, cd-hit, foldseek). TM-align at `/usr/bin/TMalign`.

**Reference biology:** Vallières et al. *Metallomics* 2024 (Fe–S review); Wehrspan et al. 2022
(AlphaFold Fe–S coordination geometry); MetalPredator.

---

## SCOPE DECISIONS (locked at project start)

1. **Taxonomic scope:** the genus *Mycobacterium* as circumscribed by GTDB. All 277
   species-representative genomes analysed fall in `g__Mycobacterium` under GTDB R232. Under NCBI
   taxonomy the same set is distributed among *Mycobacterium*, *Mycolicibacterium* and
   *Mycobacteroides*; GTDB circumscription is used throughout and all counts refer to it.
2. **Genome set:** GTDB species-representative genomes (R232, 15 Apr 2026), one per 95 % ANI species
   cluster — avoids *M. tuberculosis* oversampling. Prevalence is therefore expressed as the
   percentage of species, not of genomes sampled.
3. **Annotation:** uniform re-annotation with Bakta, not as-deposited NCBI annotations.
4. **Definition of an Fe–S protein:** a protein that coordinates an Fe–S cluster ([2Fe–2S],
   [3Fe–4S], [4Fe–4S]) as cofactor — clients plus carriers; SUF/ISC biogenesis components retained
   and tagged.
5. **Two discovery tiers:** a seed-based core (Steps 1–10) and a Pfam–GO expansion (Step 11), kept
   separable throughout by a `provenance` column.

---

# STEP 1 — Curated, metadata-rich seed dataset

**Goal:** a locked, provenance-tracked seed FASTA with per-protein metadata.

### 1.1 Source
Curated Fe–S protein supplementary (Vallières lineage), 9 sheets: 1 SUF-machinery matrix plus 8
species client sheets — Mtb (103), Mab (136), Mav (107), Mmar (108), Msm (136), Mul (75), Mbv (70),
Mle (50). Identification method = "UniProt + AlphaFold" / Estellon et al. (structurally vetted).

### 1.2 Extract and validate accession list
```bash
libreoffice --headless --convert-to csv --outdir . seeds.xlsx
awk -F',' 'NR>1 {print $2}' seeds.csv | sed 's/[[:space:]]//g' | grep -v '^$' > ids_raw.txt
sort -u ids_raw.txt > ids_unique.txt
grep -vE '^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$' ids_unique.txt
```

### 1.3 Fetch sequences (UniProt ID Mapping API)
```bash
JOB=$(curl -s --request POST 'https://rest.uniprot.org/idmapping/run' \
  --form 'from=UniProtKB_AC-ID' --form 'to=UniProtKB' \
  --form "ids=$(paste -sd, ids_unique.txt)" | grep -o '"jobId":"[^"]*"' | cut -d'"' -f4)
curl -s "https://rest.uniprot.org/idmapping/uniprotkb/results/stream/$JOB?format=fasta&compressed=false" -o seeds_raw.fasta
```

### 1.4 Reconcile requested versus returned
```bash
grep '^>' seeds_raw.fasta | cut -d'|' -f2 | sort -u > ids_returned.txt
comm -23 ids_unique.txt ids_returned.txt > ids_missing.txt
```
Each missing accession resolved individually as merged, demerged or obsolete.

### 1.5 Family-coverage confirmation
All major Fe–S families present; WhiB1–7 fully covered; aconitase present under "aconitate
hydratase" (9); 7 suf/isc/nfu entries retained and tagged as biogenesis.

### 1.6 Rename headers to stable IDs
```bash
seqkit replace -p '^\w+\|(\w+)\|.*$' -r '{kv}' -k header_map_clean.tsv seeds_raw.fasta > seeds_renamed.fasta
```

### 1.7 Deduplicate at 100 % identity only
```bash
seqkit rmdup -s seeds_renamed.fasta -o seeds_locked.fasta -D seq_duplicates.txt
```
774 input → **717 unique sequences** (57 collapses, byte-identical Mtb/Mbv orthologues and a few
Mul/Mmar pairs; recorded, not discarded).

**Deliverables:** `seeds_locked.fasta` (717), `master_metadata.tsv` (774 rows), `seq_duplicates.txt`.

---

# STEP 2 — Genome collection and uniform annotation

### 2.1 GTDB representative set
GTDB release R232 (15 Apr 2026); species representatives retrieved and filtered to the genus.

### 2.2 Download genomes
```bash
datasets download genome accession --inputfile got_accessions.txt --include genome
```

### 2.3 Quality filter (CheckM2)
Genomes retained at completeness ≥ 95 % and contamination ≤ 5 %.

### 2.4 Uniform annotation (Bakta)
```bash
bakta --db .../bakta_db/db --genus Mycobacterium --threads 16 <genome>.fna --output <out>
```
**277 genomes** uniformly annotated (Bakta DB v6.0). CDS range 3,186–7,862, mean 5,490.

**Documented limitation:** *M. leprae* (GCF_000195855.1) yields 3,989 CDS but only 25 Bakta
pseudogenes against roughly 1,300 curated pseudogenes — Pyrodigal calls intact ORFs over degraded
regions. Pseudogene-derived calls are therefore under-flagged in the leprosy bacilli.

---

# STEP 3 — Family-aware HMM construction

### 3.1 Cluster seeds
```bash
mmseqs easy-cluster seeds_locked.fasta seeds_clust tmp --min-seq-id 0.30 -c 0.70
```
210 clusters: 80 singletons, 42 doubletons, 19 triplets, 17 quads (158 thin); 52 clusters ≥ 5
members. Largest 26 members.

### 3.2 Build custom HMMs for the 52 large clusters
```bash
mafft --auto cluster.faa > cluster.aln
trimal -in cluster.aln -out cluster.trim -gt 0.2
hmmbuild cluster.hmm cluster.trim
```
Gathering thresholds calibrated by scanning each HMM back against the seed set (lowest GA 64.7 bits).

### 3.3 Cover the 158 thin clusters with Pfam HMMs
17 Pfam Fe–S families fetched: Fer2, Fer4, Fer4_2/6/7/8/9/10/12/17, Rieske, Rieske_2, Radical_SAM,
Molybdop_Fe4S4, Whib, SUFBD_core, UDG.

### 3.4 Orphan triage
6 unclustered orphans; 2 covered by existing HMMs, 4 carried forward for structural adjudication
(Step 10.3).

### 3.5 Final library
`fes_final.hmm` = **69 HMMs** (52 custom + 17 Pfam), all GA-tagged.

---

# STEP 4 — Proteome scan

```bash
for prot in proteomes/*.faa; do
  hmmsearch --cpu 16 --cut_ga \
    --tblout scan_tbl/$(basename $prot .faa).tbl \
    --domtblout scan_dom/$(basename $prot .faa).dom \
    fes_final.hmm $prot
done
```
44,260 raw rows → **25,996 unique proteins** after best-hit resolution (one best HMM per protein).

---

# STEP 5–7 — Family classification, false-positive demotion, pseudogene flag, lock

### 5.1 Map HMM to family
`hmm_family_map.tsv` resolves each HMM; locus-tag-named HMMs resolved via UniProt names and Pfam
hmmscan.

### 5.2 Demote non-Fe–S families on independent KEGG evidence
- **ParA_ATPase** (655 hits, K03496) — not Fe–S, demoted.
- **PrpD** methylcitrate dehydratase (10 hits, K01720) — not Fe–S, demoted.
- **DUF779** (6 hits, K09959) — retained as a candidate pending structural adjudication (Step 10.1).

### 5.3 Pseudogene flag
3,666 pseudogene tags across genomes, yielding 19 pseudogene-derived hits dataset-wide.

### 5.4 Lock
`all_hits_final.tsv` = 25,312 proteins, 35 families.

---

# STEP 8 — Coordination analysis, sequence layer

### 8.1 Family-aware Cys/His motif scan
`coordination_scan3.py` applies a per-family expected coordinating-cysteine count reflecting real
biology rather than a uniform four: radical SAM (3 Cys + SAM), ferredoxin [4Fe–4S] (3 Cys ± Asp),
aconitase and IPM isomerase (3 Cys + open substrate site), Rieske (2 Cys + 2 His), glutamate
synthase [3Fe–4S] (3 Cys); genuine 4-Cys families are WhiB, molybdoenzyme and glycosylase. Two motif
tiers per family (canonical spacing, relaxed variant).

### 8.2 Confidence assignment
high 18,050 (71 %), medium 6,473 (26 %), low 789 (3 %).

### 8.3 Cluster-type distribution at this stage
[4Fe–4S] 15,476, [2Fe–2S] 7,476, multi 1,129, assembly 556, [3Fe–4S] 396, FAD 273, unknown 6.

### 8.4 Refinement
`coordination_refined.tsv`: 24,717 Fe–S supported, 326 ambiguous, 269 likely non-Fe–S;
structural-target list of 601.

### 8.5 Methodological caveat (load-bearing)
Cysteine **count** is a screen, not proof of Fe–S binding, which depends on three-dimensional
Sγ–Sγ distance. Too few cysteines is a sound demotion floor; enough cysteines is only
coordination-consistent, never confirmation. Step 10 and Step 13 are the arbiters. This caveat is the
reason coordination is re-tested structurally in Step 13 rather than taken from this step.

---

# STEP 9 — Sequence similarity networks

### 9.1 Representative selection
```bash
cd-hit -i ssn_input_all.faa -o ssn_rep70.faa -c 0.70 -n 5 -M 4000 -T 8
```
Representatives at 70 % identity. Rebuilt on the v3 atlas in Step 11 (see 11.4), giving the 2,234
representatives reported in the manuscript.

### 9.2 EFI-EST and Cytoscape
`ssn_rep70.faa` submitted to EFI-EST; node IDs assigned as `zzz####` with `genome|locus_tag`
preserved in the Description column. XGMML downloaded and laid out in Cytoscape (Prefuse
force-directed). Node attributes (family, super-family, coordination status, cysteine count,
predicted cluster) built by `make_family_attr.py` and keyed on `shared name`.

### 9.3 Threshold selection
Networks were generated across a range of alignment-score thresholds and the threshold adopted for
the reported network was selected as the lowest value at which the major families separated into
discrete components without fragmenting internally. The adopted threshold is recorded in the
Cytoscape session files.

### 9.4 Result
Families separate into discrete clusters, providing independent validation of the HMM
classification. Per-family networks were generated for families with ≥ 50 representatives;
families too conserved to yield enough representatives (SDH_FRD, glutamate synthase, nitrite
reductase, aconitase, IspG, IPM isomerase) are reported as such, the conserved-essential versus
diverse-accessory split being itself informative.

---

# STEP 10 — Structural validation of the seed-derived core

**Tools:** HHpred (fold identification), AlphaFold3 (models), TM-align (superposition), AlphaFill
(cofactor transplant), foldseek, BioPython MMCIFParser for Sγ–Sγ geometry.

### 10.1 DUF779 → IscA/HesB assembly fold
6 proteins, one homogeneous family (158–166 aa). HHpred 98.8 % to the HesB/IscA Fe–S biogenesis
fold (PF05610, 2P2E, SCOP HesB-like, IscA/YfhF, ECOD Fe-S_biosyn, PF01521). TM-align of the
AlphaFold3 model against IscA 1R94 gives TM 0.64 and against 2P2E TM 0.56, both above the 0.5
same-fold threshold at 12–16 % sequence identity. The monomer shows one close cysteine pair
(3.8 Å) with the remainder 17–38 Å apart, consistent with the inter-subunit coordination
characteristic of the IscA family, which a monomer cannot display. AlphaFill transplanted no
cofactor, as is usual for apo IscA structures. **Verdict:** IscA/HesB-like fold confirmed by two
independent routes; a putative Fe–S assembly scaffold. Binding is inferred, not proven, and the
family is reported throughout as `DUF779_FeS_candidate`.

### 10.2 DNA glycosylase reclassification
The defining HMM is Pfam UDG, a function-based rather than cofactor-based model. HHpred on four
representatives spanning 2–10 cysteines returned the UDG fold (SCOP c.18.1) against non-Fe–S
templates in every case, including the 9- and 10-cysteine representatives — scattered cysteines
with no iron-sulfur cluster loop. A regex scan for the FCL motif (`C.{3,8}C.{1,4}C.{2,8}C`) across
the 283 members carrying ≥ 4 cysteines returned **0 of 283**. MutY and Nth representatives, by
contrast, returned the HhH-GPD fold (SCOP a.96) against EndoIII and MutY templates with bound
[4Fe–4S]. On these three converging lines the broad UDG `DNA_glycosylase` family was demoted and
the genuine Fe–S glycosylases MutY and Nth retained as separate families.

### 10.3 Orphan adjudication
| Orphan | aa | HHpred fold | Verdict |
|---|---|---|---|
| BCG_2352 | 128 | molybdoenzyme / formate dehydrogenase [4Fe–4S] (2IV2) | genuine Fe–S |
| MMAR_0590 | 102 | UDG fold, 5 Cys, no FCL | non-Fe–S |
| Rv0246 | 436 | MFS membrane transporter | non-Fe–S |
| Rv0311 | 409 | DUF6125 / TRAPP / H-NOX | non-Fe–S |

One of four confirmed; three reclassified.

### 10.4 SerA_dehydrogenase — retained, flagged
The sequence-layer scan flagged 135 of 136 SerA proteins as carrying ≤ 2 cysteines, and KEGG
annotates SerA as a non-Fe–S phosphoglycerate dehydrogenase. The family is nonetheless **retained**
in the atlas, for three reasons. First, the family derives from the curated seed set, whose entries
were structurally vetted upstream, so demotion would discard curated evidence on a sequence-count
screen that Step 8.5 explicitly identifies as insufficient grounds for family-level demotion.
Second, the UDG case (10.2) established the evidentiary bar for demoting a whole family — HMM
identity, an independent motif scan and fold-level evidence had to converge — and SerA meets only
the first of those. Third, retention is conservative in the direction that matters: the family is
reported throughout as coordination-unconfirmed, it is excluded from every claim that rests on
confirmed coordination, and its status is carried in the `confidence` column of
`06_coordination/coordination_refined_v2.tsv` so that any downstream user can exclude it. Readers
wishing to exclude SerA can do so without recomputing the atlas; the effect is to reduce the atlas
by one family and by the SerA protein count.

### 10.5 Scope of validation
Families built from curated UniProt Fe–S seeds are literature-grounded, and definitional Fe–S Pfam
families (Fer2, Fer4, Rieske, Radical_SAM, Molybdop_Fe4S4, Whib) have Fe–S binding intrinsic to the
domain. UDG was the unique risk: the only family defined by a function-based model, which
over-recruited non-Fe–S members. Validation at this step was therefore risk-targeted rather than
blanket; all 44 families are tested structurally in Step 13.

**Core tier after Step 10:** **24,480 proteins, 34 families** (`all_hits_final_v2.tsv`).

---

# STEP 11 — Pfam–GO expansion and the v3 atlas

**Goal:** recover divergent or under-annotated Fe–S proteins invisible to the seed-derived HMMs.

### 11.1 Define the Fe–S Pfam set
Pfam families were selected by gene ontology annotation via pfam2go, retaining those carrying
Fe–S-binding GO terms (`GO:0051536` iron–sulfur cluster binding and its descendants
`GO:0051537`, `GO:0051538`, `GO:0051539`). The resulting HMM set is
`03_hmms/extended_pfam/feS_new_families.hmm`, with the family mapping in
`03_hmms/extended_pfam/all_feS_pfam.tsv`.

### 11.2 Scan and merge
The 277 proteomes were rescanned with the expanded set and hits merged with the core tier,
retaining the better-scoring assignment per protein and recording the discovery route in a
`provenance` column (`core_validated` / `extended_pfam`).

### 11.3 Validate the expansion
Expanded proteins were required to pass the same coordination screen as the core, and their
placement was checked in the sequence-similarity network: expanded proteins integrate within the
same components as seed-derived proteins rather than forming satellite clusters.

### 11.4 v3 atlas
**28,789 proteins across 44 families** — 24,480 (85.0 %) `core_validated` and 4,309 (15.0 %)
`extended_pfam`. Ten families derive all their members from the expansion, including IspH (LytB),
NadA, ThiC, IlvD, SufU, CDGSH/NEET and NifU. The SSN was rebuilt on this atlas, giving 2,234
representatives.

### 11.5 Prevalence tiers
Because one genome is retained per species, prevalence is the percentage of species encoding a
family: **21 core** (≥ 95 % of species), **18 accessory** (10–95 %), **5 rare** (< 10 %); nine
families are present in every species sampled.

**Deliverables:** `04_results/all_hits_final_v3.tsv`, `04_results/fes_family_matrix_v3.tsv`,
`05_sequences/fes_functional_28789.faa`, `data/TableS4_fes_atlas_28789_proteins.csv`.

---

# STEP 12 — Functional and pathway classification

Functional annotation of all Fe–S proteins was performed with **eggNOG-mapper v3.3.0** against the
**eggNOG v7** orthology database. Clusters of orthologous groups (COG) categories, KEGG pathway and
module assignments and EC numbers were transferred from the best-scoring orthologue. COG categories
were assigned to 99.8 % of proteins and KEGG pathways, excluding generic global-metabolism maps, to
52.5 %.

Each family was additionally assigned to one of **19 curated functional categories** (electron
transfer, energy metabolism, redox regulation, cofactor biosynthesis, amino-acid biosynthesis, DNA
repair and so on) through `09_functional/family_function_map.tsv`. Category assignment is at family
level and is independent of the eggNOG annotation, so the two provide mutually independent views of
the same proteome.

**Deliverables:** `09_functional/fes_v3_eggnog.tsv`, `09_functional/fes_v3_categorized.tsv`,
`09_functional/csv_functional_categories_v3.csv`, `09_functional/family_function_map.tsv`.

---

# STEP 13 — Structural modelling, coordination geometry and fold clustering

### 13.1 Representative modelling
A representative sequence from each of the 44 families was modelled with **AlphaFold3**. Models were
generated **without the Fe–S cofactor** — that is, as apo structures — so that the coordination test
below is a property of the predicted protein geometry rather than of a modelled ligand. Model
quality is recorded as mean pLDDT and pTM in `13_structural/model_plddt.tsv` and
`13_structural/model_confidence.tsv`.

### 13.2 Coordination-geometry validation
Coordination was assessed geometrically rather than by sequence count, addressing the caveat
recorded in Step 8.5. A family was scored as presenting a bona fide coordination site when the
model showed **three or more cysteine sulfur atoms converging within 7.5 Å**, or, for Rieske-type
families, a **2Cys–2His** arrangement. The 7.5 Å window follows the Sγ–Sγ separations reported for
experimentally determined [4Fe–4S] and [2Fe–2S] sites (Wehrspan et al. 2022).

**Result:** 37 of 44 families (84 %) present a coordination site — 34 by cysteine triad and 3 by the
Rieske arrangement. Of the remaining seven, two are cluster-delivery scaffolds, which bind their
cofactor transiently and are not expected to present a static site, and five are
coordination-unconfirmed. SerA_dehydrogenase (10.4) is reported among the unconfirmed set.

### 13.3 All-against-all structural comparison
Pairwise structural similarity among the 44 representatives was quantified by all-against-all
alignment with **TM-align**, each pair normalised by the shorter chain and the resulting matrix
symmetrised. The matrix is `13_structural/tmalign_similarity_matrix.csv`.

### 13.4 Fold clustering
Distances were taken as 1 − TM-score and clustered by **complete linkage**. Complete linkage was
chosen deliberately: average linkage chained on this matrix, producing a single group of 34–39
families at every cut tested, because many Fe–S folds share small similar substructures that average
linkage treats as grounds for merging. Complete linkage requires every member of a group to exceed
the threshold against every other member, which is the semantics that "same fold" demands.

The dendrogram was cut at **TM-score 0.5**, the conventional threshold for shared fold assignment
(Xu & Zhang 2010), giving **8 fold groups** and leaving **25 families with no structural relative**.
Within-group similarity far exceeds between-group similarity (mean TM-score 0.647 against 0.277;
Wilcoxon signed-rank, p < 0.001).

**Deliverables:** `13_structural/`, `fig4_fold_groups.csv`, AlphaFold3 model coordinates.

---

# STEP 14 — Phylogenetics and residue conservation

### 14.1 Species tree
A concatenated alignment of conserved single-copy marker genes was generated across the 277
genomes and used to infer the species tree with **IQ-TREE 2**, the substitution model selected by
**ModelFinder** and branch support assessed by 1,000 ultrafast bootstrap replicates. The tree is
`12_phylogenetics/species_tree.treefile`.

### 14.2 Family trees
```bash
mafft --auto family.faa > family.aln
trimal -in family.aln -out family.trim -gt 0.2
iqtree2 -s family.trim -m MFP -B 1000 -T AUTO --prefix family
```
Per-family maximum-likelihood trees were inferred for the families with sufficient membership and
are in `12_phylogenetics/*.treefile`. The pipeline script is
`12_phylogenetics/run_family_phylo.sh`.

### 14.3 Gene-tree / species-tree congruence
Congruence was quantified as the **Spearman rank correlation between patristic distances** — for
every pair of genomes represented in a family tree, the path-length distance in that tree against
the path-length distance in the species tree. This measure is independent of rooting and tolerates
paralogy, both of which defeat clade-retention measures on unrooted trees carrying multiple copies
per genome. Where a genome contributes several paralogues, the minimum pairwise distance is used.

### 14.4 Residue conservation
Per-column conservation was computed from the trimmed family alignments, with the proportion of
gap-free sequences recorded alongside so that conservation at sparsely occupied columns is not
over-read. Columns corresponding to predicted cluster-coordinating cysteines were identified from
the Step 13.2 geometry and mapped onto the alignment.

### 14.5 Statistical treatment of conservation
For each family, the mean per-column conservation of predicted coordinating cysteine positions was
compared with that of all remaining alignment columns, giving **one paired difference per family**.
Significance within a family was assessed by Mann–Whitney U with **Benjamini–Hochberg** correction
across the 37 families with resolved coordination sites; consistency of the direction of effect
across families was assessed by a **two-sided sign test**.

This design replaces an earlier pooled comparison that treated every alignment column in every
family as an independent observation. That pooling is pseudoreplication — columns within a family
are not independent — and it produced an implausibly small p value. Only the per-family paired
analysis is reported.

**Result:** the difference is positive in all 37 families (median +0.202, range +0.069 to +0.518;
two-sided sign test, p = 1.5 × 10⁻¹¹) and individually significant in 31 of 37 after
Benjamini–Hochberg correction.

**Deliverables:** `11_conservation/conservation.py`, `11_conservation/conservation_summary.tsv`,
`12_phylogenetics/`.

---

# STEP 15 — *M. tuberculosis* essentiality, selectivity and phylogenomic distribution

### 15.1 Essentiality mapping
The Fe–S proteins of *M. tuberculosis* H37Rv were cross-referenced with the genome-wide
transposon-mutagenesis essentiality calls of DeJesus et al. (2017). Of the 72 Fe–S proteins
encoded by H37Rv, 57 mapped to Rv numbers and 55 carried definitive calls; 15 could not be mapped
and 2 were uncertain. Mapping is recorded in `10_tb_analysis/mtb_fes_rv_mapped.tsv`.

### 15.2 Enrichment test
Essential and growth-defect calls among Fe–S proteins were compared with the genome-wide reference
by **two-sided Fisher's exact test** with the Haldane–Anscombe correction applied to the odds ratio.
Both arms were evaluated on the same category definition and genes with uncertain calls were
excluded from both.

**Result:** 16 of 55 (29.1 %) against 596 of 3,943 genome-wide (15.1 %) — a 1.93-fold enrichment
(odds ratio 2.30, 95 % CI 1.28–4.15; p = 7.6 × 10⁻³). A sensitivity analysis across three
alternative category groupings gave fold enrichments of 1.84 to 2.02 (all p < 0.025).

### 15.3 Target selectivity
The Fe–S proteins of H37Rv were searched against the human reference proteome (UniProt
**UP000005640**) with **DIAMOND blastp** in `--very-sensitive` mode. A protein was considered to
lack a detectable human homologue when no alignment reached **E ≤ 1 × 10⁻³**; no identity or
coverage threshold was applied. The absence of secondary filters is deliberate — the question is
whether any detectable homologue exists, and an identity or coverage cut-off would convert weak but
real homologues into apparent absences, inflating the selectivity claim.

**Result:** 11 of the 16 essential or growth-defect proteins return no alignment at this threshold.

### 15.4 Allometry of Fe–S content
Fe–S complement was regressed against predicted proteome size across the 277 genomes by ordinary
least squares, with outliers identified by studentised residual, and separately by a
negative-binomial fit on log–log scale to estimate the scaling exponent.

**Result:** Spearman ρ = 0.77 (p = 1.3 × 10⁻⁵⁴); OLS R² = 0.64; 18.8 additional Fe–S proteins per
1,000 predicted proteins; negative-binomial scaling exponent 1.07, with Fe–S proteins constituting
1.89 ± 0.22 % of the proteome. *M. leprae* and *M. lepromatosis* are the only genomes beyond a
studentised residual of −3.

### 15.5 Phylogenomic distribution
Genomes were grouped by cutting the species tree at a depth chosen by bisection to yield the target
number of clades. Repertoire similarity between genomes was computed as **Jaccard distance** on
family presence/absence and ordinated by **principal coordinates analysis**. Family-level
enrichment between groups was tested by Fisher's exact test with Benjamini–Hochberg correction.

**Phylogenetic control.** Clinically relevant species are concentrated in a single clade, so any
enrichment found across all 277 genomes may reflect shared ancestry rather than clinical status. The
identical test was therefore repeated **within** that clade, where phylogeny is approximately held
constant. A family that remains significant within the clade is associated with clinical status; one
that loses significance was a clade marker. This control is implemented in `check_fig8_confound.py`
and its output is reported alongside the primary test.

**Deliverables:** `10_tb_analysis/`, `14_taxonomy/`, `fig8_genome_coords.csv`,
`fig8_within_clade_test.csv`.

---

## Software

| Tool | Purpose |
|---|---|
| Bakta | uniform genome annotation |
| CheckM2 | genome quality filtering |
| MMseqs2 | seed clustering |
| MAFFT | multiple sequence alignment |
| trimAl | alignment trimming |
| HMMER (hmmbuild, hmmsearch, hmmfetch, hmmscan) | profile HMM construction and search |
| CD-HIT | SSN representative selection |
| EFI-EST, Cytoscape | sequence-similarity networks |
| HHpred | remote fold identification |
| AlphaFold3 | structure prediction |
| TM-align | structural superposition |
| AlphaFill | cofactor transplantation |
| foldseek | structure search |
| IQ-TREE 2, ModelFinder | phylogenetic inference |
| eggNOG-mapper v3.3.0 (eggNOG v7) | functional annotation |
| DIAMOND | human homology search |
| SeqKit, Biopython, pandas, SciPy, scikit-learn | data handling and statistics |

Exact version strings and full command lines for each step are recorded in the step directories.

---

*All counts in this document refer to the v3 atlas and are traceable to
`04_results/all_hits_final_v3.tsv`, `09_functional/fes_v3_categorized.tsv` and
`06_coordination/coordination_refined_v2.tsv`.*
