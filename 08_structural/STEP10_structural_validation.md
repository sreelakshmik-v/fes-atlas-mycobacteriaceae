# Step 10 — Structural Validation (Depth A)

## Corrected repertoire
- 24,480 functional Fe-S proteins (was 25,312); 34 families (was 35); 13 core (was 14)

## DUF779 -> IscA/HesB Fe-S assembly fold (kept)
- 6 proteins, one family (158-166 aa). HHpred 98.8% to HesB/IscA fold (PF01521, SCOP b.124.1; templates 1R94/2P2E/1NWB).
- TMalign vs IscA(1R94) TM=0.64; vs 2P2E TM=0.56 (both >0.5 = same fold).
- Cys geometry (monomer): one close pair (3.8A), rest dispersed; consistent with IscA inter-subunit (dimer) coordination.
- AlphaFill: no cofactor transplanted. Binding inferred, not proven.
- Known fold (not novel). HesB Fe-S proteins known in mycobacteria (M. smegmatis MSMEG_4272). Contribution = DUF779-family fold assignment. Verify vs InterPro/Pfam before publication.

## Broad DNA_glycosylase (UDG) -> non-Fe-S uracil glycosylase (demoted)
- Defining HMM = Pfam UDG (function-based, non-Fe-S). HHpred on 4 reps (2-10 Cys): all UDG fold (SCOP c.18.1), non-Fe-S.
- FCL motif scan: 0 of 283 >=4-Cys members have the iron-sulfur cluster loop.
- 832 functional proteins demoted (status nonFeS_UDG_demoted).
- Genuine Fe-S glycosylases retained: MutY (187), Nth (172) = HhH-GPD fold (SCOP a.96) with [4Fe-4S] (1ORN/2ABK/3N5N, HET:SF4).

## Orphans (4 candidate-novel)
- BCG_2352 (128aa): molybdoenzyme/FDH [4Fe-4S] fold -> Fe-S
- MMAR_0590 (102aa): UDG fold, no FCL -> non-Fe-S
- Rv0246 (436aa): MFS membrane transporter -> non-Fe-S
- Rv0311 (409aa): DUF6125/TRAPP -> non-Fe-S

## SerA_dehydrogenase: 135/136 <=2 Cys; 529-aa phosphoglycerate dehydrogenase; KEGG non-Fe-S. Flagged.

## Methodological note
Cysteine count is a screen, not proof (binding depends on 3D distance). "Too few Cys" = sound demotion floor; "enough Cys" = coordination-consistent only. UDG case: 9-10 Cys proteins still non-Fe-S (no FCL). Structural fold = arbiter.

## Tools: HHpred, AlphaFold3, TMalign, AlphaFill, foldseek, BioPython.
