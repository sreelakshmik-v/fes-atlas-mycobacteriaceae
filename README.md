# Fe-S Cluster Protein Atlas of *Mycobacteriaceae*

Genome-wide identification, classification, coordination, network, and structural analysis
of iron-sulfur (Fe-S) cluster proteins across 277 *Mycobacteriaceae* species-representative
genomes (GTDB R232).

## Headline result
**24,480 functional Fe-S proteins** in **34 families** (13 core, present in >=95% of genomes),
identified by 69 family-aware HMMs, classified by coordination motif, validated by
sequence-similarity networks and structural fold analysis (HHpred / TMalign / AlphaFold3).

## Repository layout
| Folder | Contents |
|---|---|
| 01_seeds/ | 717 curated Fe-S seed proteins + metadata |
| 02_genomes/ | 277 genome accession list (genomes via NCBI datasets) |
| 03_hmms/ | 69-HMM library + HMM-to-family map |
| 04_results/ | Final hit table, family x genome matrix, abundance/presence CSVs |
| 05_sequences/ | FASTA of all 24,480 predicted Fe-S proteins |
| 06_coordination/ | Per-protein coordination (cluster type, confidence) |
| 07_ssn/ | SSN representative sequences + Cytoscape attributes |
| 08_structural/ | Structural validation (DUF779 model, Step 10 summary) |
| scripts/ | Analysis scripts |
| figures/ | Abundance bar, clustered heatmap, presence/absence heatmap |

## Key files
- 04_results/all_hits_final_v2.tsv — genome, locus, HMM, score, family, cluster_type, status
- 05_sequences/fes_functional_24480.faa — all predicted Fe-S protein sequences
- 03_hmms/hmm_family_map.tsv — HMM -> family -> cluster-type mapping

## Pipeline
Curated seeds (717) -> 277 GTDB genomes (Bakta) -> 69 HMMs (mmseqs clusters + Pfam)
-> hmmsearch scan -> classify (35 families) / demote non-Fe-S / lock -> coordination analysis
-> SSN (EFI-EST/Cytoscape) -> structural validation -> corrected 24,480 / 34 families.

See METHODS.md for full step-by-step detail. Template adapted from the Alba-domain SSN
study (Sci Rep 2024, s41598-024-79937-4).

## Status
Half 1 (identification / classification / validation): complete.
Half 2 (functional mapping: metabolism / redox / virulence): in progress.
