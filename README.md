# Fe-S Cluster Protein Atlas of *Mycobacteriaceae*

Genome-wide identification, classification, coordination, network, structural, and
functional analysis of iron-sulfur (Fe-S) cluster proteins across 277 *Mycobacteriaceae*
species-representative genomes (GTDB R232).

## Headline result
**28,789 Fe-S proteins across 44 families** (21 core, present in >=95% of genomes).

The atlas has two confidence tiers (see the `provenance` column in 04_results/all_hits_final_v3.tsv):
- **core_validated (24,480 proteins, 34 families):** seed-derived families validated through
  coordination analysis, sequence-similarity networks, and structural fold assignment.
- **extended_pfam (4,309 proteins, 10 families):** comprehensive Pfam-based detection
  (GO-defined Fe-S domains via pfam2go), coordination-validated. Includes the drug targets
  IspH/LytB and IlvD, the CDGSH/NEET 2Fe-2S class, quinolinate synthase (NadA), ThiC,
  and ISC biogenesis machinery (SufU/NifU).

## Repository layout
| Folder | Contents |
|---|---|
| 01_seeds/ | 717 curated Fe-S seed proteins + metadata |
| 02_genomes/ | 277 genome accession list |
| 03_hmms/ | Seed-derived 69-HMM library + extended_pfam/ (GO-defined Fe-S Pfam set) |
| 04_results/ | v3 hit table (tiered), family x genome matrix, abundance/presence CSVs |
| 05_sequences/ | FASTA of all 28,789 predicted Fe-S proteins |
| 06_coordination/ | Per-protein coordination (cluster type, confidence) |
| 07_ssn/ | SSN representative sequences + Cytoscape attributes |
| 08_structural/ | Structural validation (DUF779 model, Step 10 summary) |
| 09_functional/ | Family->function map, categorized proteins, functional categories, species map |
| scripts/ | Analysis scripts |
| figures/ | Family abundance (tier-colored), heatmap, presence/absence |

## Key files
- 04_results/all_hits_final_v3.tsv — genome, locus, HMM, score, family, cluster_type, status, provenance
- 05_sequences/fes_functional_28789.faa — all predicted Fe-S protein sequences
- 09_functional/family_function_map.tsv — family -> metabolic/redox/virulence category

## Pipeline
Curated seeds (717) -> 277 GTDB genomes (Bakta) -> HMM library (mmseqs clusters + Pfam)
-> hmmsearch scan -> classify / demote non-Fe-S / lock -> coordination -> SSN -> structural
validation (core, 24,480) -> comprehensive Pfam expansion (extended, +4,309) -> functional
mapping (metabolism / redox / virulence) + TB-vs-non-TB comparison.

See METHODS.md for full detail. Template: Alba-domain SSN study (Sci Rep 2024, s41598-024-79937-4).

## Functional landscape (28,789 proteins)
Electron transfer + energy metabolism dominate (~49%); the Fe-S proteome is enriched in
TB-relevant systems: respiration (NADH/SDH/QcrA), the MEP isoprenoid pathway (IspG + IspH,
drug targets), cholesterol catabolism Rieske oxygenases, WhiB redox regulators (WhiB3/6/7),
F420 biosynthesis, and SUF/ISC Fe-S biogenesis.

## Status
Half 1 (identification / classification / validation): complete.
Half 2 (functional mapping): in progress.
