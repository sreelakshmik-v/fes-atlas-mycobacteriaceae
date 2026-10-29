# Fe–S Cluster Protein Atlas of *Mycobacterium*

Data and code accompanying the manuscript:

> Sreelakshmi K V, Sarode VV, Shenoy G G, Dehury B. *Unified by chemistry, not fold: A genus-wide
> atlas of the mycobacterial iron–sulfur proteome and its antitubercular target space.*

Genome-wide identification, classification, coordination, network, structural and functional
analysis of iron–sulfur (Fe–S) cluster proteins across 277 species-representative genomes of the
genus *Mycobacterium* (GTDB R232).

## Headline result

28,789 Fe–S proteins across 44 families (21 core, present in >=95 % of species).

Two confidence tiers (see the `provenance` column in `04_results/all_hits_final_v3.tsv`):

- **core_validated** (24,480 proteins, 34 families): seed-derived families validated through
  coordination analysis, sequence-similarity networks and structural fold assignment.
- **extended_pfam** (4,309 proteins, 10 families): GO-defined Fe–S Pfam domains via pfam2go,
  coordination-validated. Includes IspH/LytB, IlvD, CDGSH/NEET, NadA, ThiC and SufU/NifU.

## Repository layout

| Folder | Contents |
|---|---|
| `01_seeds/` | 717 curated Fe–S seed proteins + metadata |
| `02_genomes/` | 277 genome accession list |
| `03_hmms/` | Seed-derived HMM library + `extended_pfam/` |
| `04_results/` | v3 hit table (tiered), family x genome matrix, abundance CSVs |
| `05_sequences/` | FASTA of all 28,789 predicted Fe–S proteins |
| `06_coordination/` | Per-protein coordination (cluster type, confidence) |
| `07_ssn/` | SSN representative sequences + Cytoscape attributes |
| `08_structural/` | Structural validation summary |
| `09_functional/` | Family to function map, categorized proteins, species map |
| `10_tb_analysis/` | *M. tuberculosis* H37Rv essentiality intersection (DeJesus et al. 2017) |
| `11_conservation/` | Per-column conservation and ligand-cysteine analysis |
| `12_phylogenetics/` | Species tree (277 taxa) and per-family maximum-likelihood trees |
| `13_structural/` | AlphaFold3 model confidence and all-against-all TM-score matrix |
| `14_taxonomy/` | Species tree annotations (iTOL) |
| `data/` | Annotated atlas table (28,789 proteins), CSV and XLSX |
| `figures/manuscript/` | Final manuscript figures |
| `scripts/` | Pipeline and figure-generation scripts |

## Key files

- `data/TableS4_fes_atlas_28789_proteins.csv` — the full annotated atlas
- `04_results/all_hits_final_v3.tsv` — genome, locus, HMM, score, family, cluster type, provenance
- `13_structural/tmalign_similarity_matrix.csv` — all-against-all TM-scores for the 44 representatives
- `METHODS.md` — full methods, software versions and command lines

## Pipeline

Curated seeds (717) -> 277 GTDB genomes (Bakta) -> HMM library (mmseqs clusters + Pfam) ->
hmmsearch scan -> classify / demote non-Fe–S / lock -> coordination -> SSN -> structural
validation (core, 24,480) -> Pfam expansion (extended, +4,309) -> functional mapping ->
conservation, phylogenetics and *M. tuberculosis* essentiality analysis.

See `METHODS.md` for full detail.

## Licence

Code in `scripts/` is released under the MIT Licence; data files under CC-BY-4.0.
