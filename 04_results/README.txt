PRIMARY FILE: all_hits_final_v3.tsv (28,789 Fe-S proteins, 44 families)
  Columns: genome, locus, hmm, score, family, cluster_type, status, provenance
  The 'provenance' column marks each protein:
    core_validated  = 24,480 seed-derived, fully validated (coordination + SSN + structural)
    extended_pfam   = 4,309 added via comprehensive Pfam scan (coordination-validated)
  Filter status=='functional' for the Fe-S set; status=='nonFeS_UDG_demoted' are excluded UDG.

all_hits_final_v2.tsv = the 24,480 validated-core tier only (subset of v3's core_validated).
All CSVs/matrix here are v3 (44 families). v2 retained for the high-confidence-core reference.
