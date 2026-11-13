# Phase 3: Shafer et al. (PMID 38231165) evidence review

**Reviewed:** 12 November 2026. **Materials:** the author manuscript (`nihms-1950807.pdf`), Supplement 1 (Extended Data Figures 1–4), Supplement 2 (Methods), the retained IEDB query output and the original six-row review CSV. These source files are not distributed with this repository.

## Evidence checked against the paper

- Main article, pages 3-4: the authors made overlapping 15-amino-acid peptides spanning mutations studied in metastatic breast cancer. T-cell lines were expanded using peptide-pulsed dendritic cells and tested against individual peptides by IFN-γ ELISpot. Responses to reactive mutant peptides were compared with their WT counterparts.
- Supplement 1, Extended Data Figure 1 page 1: TP53 R248Q and R248W were tested using healthy-donor samples, including the participant labelled donor 13. The figure does **not** test TP53 R249S and is not an HCC experiment.
- Supplement 1, Extended Data Figure 2 page 2: for TP53 R248Q/R248W in donor 13, the WT peptides numbered (3), (4) and (5) show 536, 3,610 and 0 IFN-γ spot-forming units (SFU) per 200,000 cells, respectively. This agrees with the positive, positive and negative qualitative IEDB outcomes, each represented twice among the six WT-peptide assay annotations.
### Figure-level WT/mutant comparison (not an assay-ID mapping)

Extended Data Figure 2 also compares WT and mutant peptides for donor 13. For **TP53 R248W**, the values are WT(3) 536 versus Neo(3) 317, WT(4) 3,610 versus Neo(4) 2,121, and WT(5) 0 versus Neo(5) 898. For **TP53 R248Q**, they are WT(3) 536 versus Neo(3) 1,538, WT(4) 3,610 versus Neo(4) 2,767, and WT(5) 0 versus Neo(5) 892 (IFN-γ SFU per 200,000 cells, as labelled in the figure). The responses to WT(3) and WT(4) show that these responses were not uniformly mutant-exclusive. None of the observations establishes an R249S response, CD4 specificity or a particular HLA-II restriction.

- Supplement 2, Methods page 1: mutational libraries use overlapping 15-mers. This explains why several WT peptide windows may occur for one mutation. However, the supplement does not provide a direct cross-reference between exact peptide sequence, IEDB assay ID, figure panel and donor. The proposed assignment of the 15-mers to WT(3), WT(4) and WT(5) relies on their position, window order and qualitative outcomes.
- Main article, page 4, section on specificity: the expanded T-cell lines used for the TP53 WT-peptide experiments were not shown to be CD4-defined. Figure 2 and subsequent HLA-I work concern isolated **ESR1** TCR clonotypes; those restriction findings cannot be assigned to the TP53 WT assays.
- Main article, pages 5-6, HLA restriction and immunopeptidomics: HLA-I restriction and immunopeptidomics findings concern ESR1 mutants, not the TP53 WT 15-mers. The reviewed material did not establish their HLA-II binding, CD4 recognition or endogenous presentation.

## Record-level audit

| IEDB ID(s) | WT peptide sequence | IEDB result | Matching donor 13 figure pattern | ID-level assignment |
|---|---|---|---|---|
| 22792072, 22792080 | CMGGMNRRPILTIIT | Positive, Positive | WT(3) = 536 SFU | Not confirmed |
| 22792073, 22792081 | MNRRPILTIITLEDS | Positive, Positive | WT(4) = 3610 SFU | Not confirmed |
| 22792074, 22792082 | RRPILTIITLEDSSG | Negative, Negative | WT(5) = 0 SFU | Not confirmed |

**Interpretation:** The sequence matches concern *WT comparator peptides* used in a different study of TP53 R248Q/R248W. They are not experiments on the R249S mutant. The six assay annotations represent three WT sequences from one publication; assigning those sequences to two hypothetical DRB1 settings in the present analysis does not create independent HLA-supported observations.

## Attribution that remains unresolved

The supplementary figures do not identify which individual IEDB assay ID belongs to R248Q or R248W, nor do they link IDs directly to donors. The figure pattern is consistent with donor 13, but this remains an inference until checked against the full IEDB assay export, curation records or an author-supplied cross-reference. No row is therefore treated as verified at individual assay-ID level. The IEDB qualitative results have not been changed.

## Phase 3 conclusion

- **Exact WT-sequence annotation:** yes (3 unique WT sequences).
- **HCC-derived mutant R249S T-cell assay:** no matching record returned in the completed exact sequence query.
- **Experimental HLA-II binding:** not established from this publication or current IEDB results.
- **TP53 WT-specific CD4 response and DRB1*01:01/04:01 restriction:** not established.
- **Endogenous presentation of these TP53 WT peptides:** not established.
- **Model independent validation:** not possible with these records.

The Phase 2 rankings remain unchanged. A future source audit would need the full IEDB assay export or another direct cross-reference linking each ID to the donor and TP53 mutation panel. Further literature and database searches could establish whether any direct R249S mutant-peptide evidence exists. Until then, the missing evidence should remain an **unresolved evidence gap**, not be reported as experimental nonresponse. The proposed experiments address that gap but have not been performed.
