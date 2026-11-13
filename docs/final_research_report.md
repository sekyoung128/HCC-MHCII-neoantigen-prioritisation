# HCC MHC-II Mutation-Derived Peptide Prioritisation and Evidence Review

**Research report — computational case study**

**Date:** 13 November 2026

**Status:** Computational work and published-evidence review completed. Four follow-up experiments were proposed, but no laboratory validation of the selected TP53 R249S or CTNNB1 S45P peptides was undertaken.

**Scope:** Two HCC-associated missense variants, two hypothetical HLA-DR settings, overlapping WT/MUT peptides, three selection rules, an assay-evidence review and a separate retrospective HCC binding comparison.

> **Interpretation:** This is a reproducible comparison of candidate-selection rules and published evidence. It does not establish naturally presented HCC neoantigens, mutant-specific CD4 responses, superior predictive performance or vaccine efficacy.

## Abstract

**Background.** Predicted peptide–MHC class II binding is not equivalent to T-cell recognition, allele-specific restriction, or natural antigen processing and presentation. Ranking rules based on absolute peptide-binding predictions, mutation-associated shifts and unresolved experimental questions may nominate different candidates even when applied to the same input data.

**Objective.** To examine how three defined selection rules prioritise peptides derived from TP53 R249S and CTNNB1 S45P, determine what experimental evidence can be attributed to the exact sequences, and identify questions for follow-up experiments.

**Methods.** WT and mutant peptides of 15, 18 and 21 amino acids were assessed with NetMHCIIpan 4.3k (PEPTIDE input, EL+BA output) for hypothetical DRB1*01:01 and DRB1*04:01 settings. The resulting 432 peptide-state–HLA records provided 216 matched WT/MUT-window–HLA comparisons. Three rules, fixed after pilot inspection, prioritised candidates from this common pool: A favoured low absolute mutant BA rank; B favoured improvements in mutant BA rank before applying absolute-rank tie-breaks; and C sought questions arising from BA/EL disagreement, predicted core changes or HLA sensitivity. IEDB and CEDAR exact-sequence results were checked against their assay sources. A separate set of 13 HCC peptides with published HOM2 binding classifications from Repáraz et al. (2022) was compared descriptively with NetMHCIIpan 4.3k predictions. Four conditional experimental cases were then proposed.

**Results.** Among the 216 matched comparisons, three mutant candidates had BA percentile ranks at or below 5 and four had EL ranks at or below 5; 23 had positive WT-minus-MUT BA-rank shifts. In the strategy-specific top tens, A concentrated on one `(gene, HLA, predicted mutant core)` group, whereas B and C covered four and five such groups, respectively. A and B shared no top-ten entries; A and C shared none; B and C shared one (Jaccard ~0.053). Six exact-sequence T-cell assay annotations were retrieved, all corresponding to three **WT control** peptide sequences in a distinct TP53 R248Q/R248W study; no exact mutant peptide-specific HLA-II binding, confirmed CD4 response or allele-level restriction was established for the nominated R249S/S45P candidates. In the separate retrospective cohort, 10 of 13 candidates were reported as HOM2 binding-score-positive. A mutant BA rank ≤5 occurred in four of 10 reported binders **and** one of three non-binders; improved predicted BA rank was seen in eight of 10 binders **and all three** non-binders.

**Conclusion.** The three rules produced different candidate lists, but these differences cannot be interpreted as differences in immunogenicity. Careful source attribution prevented unrelated WT-control assays from being counted as mutant-specific evidence. The small, preselected historical comparison cannot independently validate the predictor. The study therefore offers a reproducible selection analysis and proposals for testing its unresolved biological questions, rather than a validated HCC neoantigen discovery.

## 1. Research question and study design

The main question was how an absolute-binding rule, a mutation-aware rule and a rule designed around experimental questions would change the candidates selected from one HCC peptide–HLA prediction set. I also examined whether the exact peptide sequences had supporting evidence for MHC-II binding, CD4 recognition, allele-specific restriction or processing-dependent presentation, and what findings would warrant further testing.

The analysis was exploratory. Phase 2 rules were fixed before the Phase 2 rankings were calculated, but the pilot and extension results had already been inspected. This was not a preregistered or blinded evaluation. Although the dataset contains 216 window–HLA comparisons, they arise from **only two variants** and are not independent mutations or patient observations. Both HLA alleles were hypothetical rather than patient-typed.

## 2. Methods

### 2.1 Mutation and peptide definition

Two HCC-associated missense variants were selected for an initial case study: TP53 R249S and CTNNB1 S45P. Matched WT and mutant peptide windows were constructed, preserving exact sequence and residue position, with peptide lengths of 15, 18 or 21 residues. The same windows were considered for DRB1*01:01 and DRB1*04:01. All predictions used NetMHCIIpan-4.3k with PEPTIDE-format input, EL+BA outputs, and context encoding off in the underlying runs. The Phase 1 extension was explicitly post-pilot.

Absolute mutant BA/EL percentile ranks were kept separate from mutation-associated changes in those ranks:

`deltaBA = WT %Rank_BA − MUT %Rank_BA`

`deltaEL = WT %Rank_EL − MUT %Rank_EL`

Lower percentile ranks indicate more favourable predictions within each model endpoint. A positive delta means a lower mutant rank; it does not measure a change in affinity or establish binding or mutant-specific CD4 recognition.

### 2.2 Three selection strategies

**Strategy A — BA-first.** Candidates were ordered by increasing mutant `%Rank_BA`, using the fixed EL and candidate-ID tie-breaks. The rule favours stronger predicted absolute binding.

**Strategy B — mutation-aware.** Candidates with `deltaBA > 0` came first, followed by absolute mutant BA rank and the specified tie-breaks. This ordering is a defined decision rule, not a trained immunogenicity score.

**Strategy C — experiment-informative.** A fixed round-robin procedure drew from three groups: opposed BA/EL changes, changes in predicted WT/MUT core offset, and differences between HLA settings. The aim was to select questions worth testing, rather than maximise a predicted immunogenicity score.

All three rules were applied to the same 216 rows. I compared their top-5, top-10 and top-20 selections, Jaccard overlap and candidate/core composition. Grouping by predicted 9-mer core gives a useful indication of redundancy, but does not identify independently processed epitopes.

### 2.3 Experimental evidence attribution

The full peptide sequences were searched against IEDB and CEDAR using separate retrieval workflows. An assay hit was not automatically counted as HLA-II or CD4 evidence. For each relevant record, the review considered the exact peptide and WT/MUT status, study context, assay endpoint, CD4 identity, experimental HLA information, links to the source figures, WT cross-reactivity and whether antigen was supplied as peptide or processed from a full-length protein. Where a search returned no eligible record, this was recorded as a gap in the retrieved evidence, not a negative experimental result.

### 2.4 Separate historical experimental comparator

Thirteen mutation-derived peptides with figure-level experimental HLA-DR/HOM2 binding classifications were extracted from Repáraz et al. (2022). Matching WT/MUT sequences were submitted to NetMHCIIpan-4.3k for DRB1*01:01, generating 26 peptide predictions. The original TXT output, tab-separated XLS, frozen manifest and figure-derived binary labels were cross-checked. Published experimental binder status used the paper's competition binding-score threshold (`BS > 0.5`); it was **not** a measured IC50 threshold. Numeric figure bar heights were not presented as extracted exact measurements.

The historical panel concerned mutations other than TP53 R249S and CTNNB1 S45P. The source study had already used an earlier NetMHCIIpan version to select its peptides, and the available material did not consistently establish the exact experimental HLA-DR1 suballele. The comparison is therefore descriptive, not a prospective or independent validation.

### 2.5 Experiment-informed follow-up

Phase 4 v2 selected three locked Phase 2 rank-one representatives (A/B/C) and one explicitly post hoc CTNNB1 case (D). For each, the proposal distinguished four conditional evidence gates: biochemical peptide–HLA-II binding, functional peptide-pulsed CD4 recognition, allele-specific HLA-II restriction, and full-length antigen processing-dependent recognition. Suggested experiments were **not performed**, and numerical experimental decision thresholds were not invented in the absence of validated SOPs and available experimental resources.

## 3. Results

### 3.1 Prediction inventory

Phase 1 generated **432 peptide-state–HLA predictions**, corresponding to **216 paired WT/MUT-window–HLA comparisons** involving 216 unique exact WT+MUT peptide sequences. Quality control reported consistent model settings and complete matched-pair comparisons. Across 216 pairs:

| Model-derived observation | Count |
|---|---:|
| Mutant `%Rank_BA ≤ 1` | 0 |
| Mutant `%Rank_BA ≤ 5` | 3 |
| Mutant `%Rank_EL ≤ 1` | 0 |
| Mutant `%Rank_EL ≤ 5` | 4 |
| Mutant BA **and** EL rank ≤ 5 | 2 |
| `deltaBA > 0` | 23 |
| `deltaEL > 0` | 118 |
| Predicted WT/MUT core offset changed | 41 |
| Mutation within predicted MUT 9-aa core | 64 |

These counts refer to predicted window–HLA records. They should not be read as numbers of independent neoantigens, patients or confirmed binders.

### 3.2 Differences between ranking objectives

| Top-10 descriptive feature | A: BA-first | B: Mutation-aware | C: Experiment-informative |
|---|---:|---:|---:|
| Distinct `(gene, HLA, predicted MUT core)` groups | 1 | 4 | 5 |
| Distinct predicted MUT 9-mer core **sequences**, ignoring gene/HLA grouping | 1 | 3 | 3 |
| Candidates with improved mutant BA rank | 0 | 10 | 2 |
| Candidates with mutant BA rank ≤5 | 3 | 0 | 0 |
| TP53 cases | 10 | 10 | 9 |
| CTNNB1 cases | 0 | 0 | 1 |

The group counts above include gene and HLA context; they are **not** counts of unique 9-mer sequences. In particular, B and C each contain only three distinct predicted 9-mer core sequences, and C consists of nine TP53 and one CTNNB1 nominations. Top-10 overlaps were A–B: `0/20` union (Jaccard 0), A–C: `0/20` (0), and B–C: `1/19` (~0.053). Given their different aims, the rules were expected to select different candidates. A favoured a narrow group of cores, while B selected candidates with better mutant-relative-to-WT rankings despite weak absolute BA ranks. These are trade-offs in candidate selection, not evidence that any rule predicts CD4 responses more accurately.

### 3.3 Exact-sequence IEDB/CEDAR evidence

The exact-sequence inventory contained **216 unique peptide sequences**. Retrieved T-cell evidence amounted to **six IEDB assay annotations** for **three TP53 WT peptide sequences**, attributed to the TP53 R248Q/R248W context reported by Shafer et al. (2024), a metastatic breast-cancer investigation. Four annotations were qualitatively positive and two negative; they were not uniquely mapped to the precise donor/figure subpanel at individual assay-ID level. These results were **not** R249S mutant peptide responses, and CD4 identity, exact HLA-II restriction and natural presentation were not verified for the sequences of interest.

The supplementary CEDAR exact search returned three epitope records and the same six T-cell assay IDs, with no additional MHC assay records retrieved. CEDAR–IEDB overlap is **not independent evidence replication**. Within the completed query scope, no exact mutant peptide-specific HLA-II binding, confirmed CD4 recognition or allele-specific restriction was established. The search documents what was found within its scope; it cannot show that the peptides fail to bind or are non-immunogenic.

### 3.4 Retrospective Repáraz et al. peptide-binding comparison

The historical comparator comprised 13 mutation-derived peptides: **10** published HOM2 binding-score-positive and **3** non-binders. All 26 WT/MUT prediction rows matched across TXT and source XLS, and all 13 WT/MUT pairs differed at one residue.

| Mutant NetMHCIIpan-4.3k criterion | Among 10 experimental binders | Among 3 experimental non-binders |
|---|---:|---:|
| BA rank ≤5 | 4 | 1 |
| BA rank ≤10 | 7 | 2 |
| EL rank ≤5 | 5 | 0 |
| Improved BA rank vs matched WT | 8 | 3 |

Individual cases show why these endpoints should be kept apart: **SQL** was reported as a non-binder despite a mutant BA rank of 4.83; **RLV** was a binder with rank 16.79; and **AQR** was a binder with rank 10.15 despite an unfavourable WT-to-MUT BA change. The fact that BA-rank improvement was also present in *all three* experimental non-binders illustrates why a favourable mutation-associated rank change should not be treated as direct binding evidence.

These results do **not** support reliable estimates of AUROC, sensitivity or specificity, or an independent model benchmark: the source investigators selected candidates using NetMHCIIpan 3.2; only three experimental non-binders were included; possible training-set overlap was not excluded; reported experimental HLA-DR1 context was not guaranteed to be exact *01:01 for every case; and gene/variant/sequence provenance discrepancies remain flagged in the source evidence matrix. Specifically, the source material contains an **ELF2 D291E versus D291J** mutation-label discrepancy and **ARHGEF40 versus ARHGEF4** notation for the SPE case; the legacy SPE HLA prediction refers to **DRB1*01:02**, not the **DRB1*01:01** used for the new comparison. These discrepancies are audit flags, not resolved biological identities. The published biochemical endpoint was an HOM2 competition-binding **score**, with **BS > 0.5** used to classify binders. It was not an IC50 measurement, CD4 response or measure of tumour presentation. Separately calculated `%Rank_BA` cutoffs must not be conflated with the NetMHCIIpan web output's `BindLevel`, which can reflect a different output/threshold context.

### 3.5 Four conditional experimental cases

| Case | Selection origin | WT/MUT example and allele | Central uncertainty |
|---|---|---|---|
| A | Phase 2 BA-first rank 1 | TP53 R249S, 15-mer beginning at residue 248, DRB1*04:01 | Relatively favourable BA/EL ranks despite no mutant-over-WT predicted binding advantage; measured binding and CD4 selectivity remain distinct |
| B | Phase 2 mutation-aware rank 1 | TP53 R249S, 15-mer beginning at 242, DRB1*01:01 | `deltaBA = +0.92` with mutant BA rank 31.87; is any biochemical WT/MUT difference reproducible? |
| C | Phase 2 experiment-informative rank 1 | TP53 R249S, 15-mer beginning at 243, DRB1*04:01 | Opposite BA/EL rank-change directions and a predicted core-offset switch, despite poor absolute ranks |
| D | **Post hoc exploratory** | CTNNB1 S45P, 15-mer with mutation at residue 15/15, DRB1*04:01 | Favourable EL-related rank versus poor BA rank and potential sensitivity to peptide-window definition |

**A** is an illustrative binding-versus-recognition case; WT/MUT core predictions are identical (`ILTIITLED`) and exclude the substituted residue. **B** tests the interpretation of a small ranking improvement rather than presuming strong binding. **C** requires initial binding feasibility before any resource-intensive experimental binding-register attribution; a computational 9-mer core shift is not itself a demonstrated physical register shift. **D** is explicitly separate from the frozen ranking comparison; alternative mutation-centred peptide windows would create *new exploratory candidates*, not a revision of Phase 2.

The v2 plan sets out a proposed comparison, necessary controls, remaining uncertainties and decision points for each case. It also records which reagents and assay conditions would need to be established before laboratory work could begin. An anti-HLA-DR blocking effect alone would not establish restriction to a particular DRB1 suballele. Recognition of full-length antigen in an engineered APC model would not automatically establish presentation by patient HCC tumours.

## 4. Discussion

The central finding is that **the choice of ranking rule changes which peptides are prioritised**. A produced more favourable predicted absolute BA ranks but little predicted WT-to-MUT binding advantage and marked core redundancy. B consistently favoured prediction shifts but nominated weaker absolute BA candidates. C surfaced questions about endpoint and HLA sensitivity but was not designed to maximise predicted binders. Overlap statistics measure the separation of these selection objectives; they are not biological performance measures.

The source-level evidence audit imposed an additional constraint: exact-sequence database hits can arise from WT controls tested in a **different variant, tumour setting and functional context**. Counting those assays as R249S-specific CD4 evidence or as confirmed DRB1 restriction would misrepresent what was observed. Likewise, missing database matches are not negative assays.

The historical comparison provided published experimental binding categories for a different, preselected set of HCC peptides. Several cases did not agree with a straightforward reading of predicted BA rank. Because peptide selection already used NetMHCIIpan 3.2, and because exact HLA specificity and model training overlap are unresolved, it provides useful cautionary examples, **not out-of-sample proof that NetMHCIIpan-4.3k generalises** or that any A/B/C selection rule is more accurate.

The Phase 4 plan proposes experiments in sequence, with further work depending on the preceding result. It does not assume that every candidate should undergo all four assays. Experimental feasibility (defined HLA reagents, compatible APCs, functional CD4 cells, replication, and SOP-defined measurement criteria) is currently unverified. Positive biochemical binding alone cannot establish mutant-specific CD4 recognition; positive peptide-pulsed T-cell responses do not demonstrate naturally processed tumour presentation.

## 5. Limitations

1. **Small number of biological variants:** 216 comparisons arise from only TP53 R249S and CTNNB1 S45P, with extensive overlap among windows. They are not independently sampled mutations, patients or tumours.
2. **Hypothetical HLA and disease context:** DRB1*01:01 and DRB1*04:01 were model input scenarios; patient genotype, tumour HLA-II expression, full-length antigen expression and processing competence were not independently established.
3. **No new empirical validation of nominated peptides:** No direct R249S/S45P mutant peptide–HLA-II binding, validated CD4 response, allele-specific restriction or endogenous presentation has been demonstrated by the project.
4. **Post-pilot analysis:** Prediction-window expansion, method specification and the later four-case experimental plan followed inspection of preliminary outputs. Locking a subsequent script and hashing inputs supports reproducibility, not prospective validation.
5. **Endpoint mismatch:** BA rank, EL rank, biochemical binding scores, CD4 functional assays, HLA restriction and antigen processing measure distinct properties and cannot be substituted for one another.
6. **Database evidence is scope-limited:** Exact-sequence IEDB/CEDAR searches and source-level attribution do not establish comprehensive absence of experimental literature or negative immunogenicity.
7. **Retrospective comparator selection bias:** The 13 historical peptides were preselected with an earlier NetMHCIIpan model; there were only three experimental non-binders, uncertain exact HLA suballele correspondence, possible training-set overlap and noted annotation discrepancies.
8. **Experimental design remains conditional:** No wet-lab SOP, sample access, reagent availability, power calculation or validated numerical go/no-go threshold was established. The v2 plan must not be represented as an executed experimental protocol.

## 6. Conclusions and project decision

This study provides a reproducible account of how three selection rules prioritise HCC mutation-derived MHC-II peptides. It shows why absolute binding predictions and WT-to-MUT ranking changes should be considered separately, how the selected lists differ, and why published assay records need to be traced back to their experimental context. The remaining uncertainties inform the proposed follow-up experiments.

The appropriate next step is to present the computational work and evidence review with their limitations, while retaining the Phase 4 v2 plan as a proposal. The findings do **not** justify claims of validated HCC MHC-II neoantigens, mutant-selective CD4 responses, a superior selection rule, generalisable NetMHCIIpan performance or vaccine efficacy. Any further evaluation would need experimental peptide–HLA-II binding and/or functional CD4 outcomes linked to exact sequences and allele information, suitable controls, and an independent cohort designed to minimise information leakage.

## 7. Reproducibility and source audit

**Primary project outputs and provenance:**

| Phase | Core source file(s) | Purpose |
|---|---|---|
| 1 | `results/phase1_extension/analysis/wt_mut_comparison_216.csv`; Phase 1 extension `analysis_qc_and_results.txt` | Frozen matched predictions and QC |
| 2 | `phase2_rules.json`; `phase2_method_lock.md`; `results/phase2/phase2_all_216_ranked.csv`; `results/phase2/phase2_topK_overlap.csv` | Rule definitions, rankings and overlap |
| 3 | Archived IEDB/CEDAR retrieval logs and `phase3_WT_Tcell_assay_review_FINAL_v2.csv` (**not distributed in this repository**) | Source-level T-cell assay attribution |
| 3C | `results/phase3/historical_comparison/phase3_historical_comparison_13.csv`; `phase3_historical_QC_and_results.txt` | Separate historical comparison |
| 4 | `results/phase4/phase4_case_registry.csv` (**not distributed**); `results/phase4/phase4_experiment_plan_v2.md`; `results/phase4/phase4_experiment_decision_matrix_v2.csv` | Proposed cases and conditional experiment decisions |

**Recorded file hashes:**

- Phase 2 source matched-comparison CSV SHA-256: `7b1cebf1d4c282a1041c7c94bd15392ee99cc8f78de9fd8c7d716415d6415b6f`.
- Phase 2 frozen rules SHA-256: `5e0c9b369777cdb9e3cae4b46dc064e57e4e1315363dbfce0cac3eeb5c024f6c`.
- Phase 4 v2 experiment plan SHA-256 recorded in QC: `d2df1879594eb149c56423c3de75547da69fdd86514dbeed3259ea4c7f1b0054`.
- Phase 4 v2 decision matrix SHA-256 recorded in QC: `3e8cdc182512900954e99fe22f66e5df07f79a1665741aa9cc64b834a64d6e41`.

**Figures included in the repository:**

- Phase 2: `results/phase2/fig1_topk_jaccard.png`, `fig2_rank_A_vs_B.png`, `fig3_absolute_ba_vs_delta.png` and `fig4_top10_membership.png`.
- Phase 3C: `results/phase3/historical_comparison/phase3_fig1_mutant_ba_by_assay_label.png` and `phase3_fig2_mutant_ba_vs_delta.png`.

Phase 1 figures from the original workspace are not included in this snapshot.

**Literature reviewed:**

- Repáraz et al. (2022), *Neoantigens as potential vaccines in hepatocellular carcinoma*, *Journal for ImmunoTherapy of Cancer*, DOI: https://doi.org/10.1136/jitc-2021-003978; main article figures and supplementary evidence tables. The experimental HLA-DR/HOM2 binding classifications are attributed to this source.
- Shafer et al. (2024), PMID: https://pubmed.ncbi.nlm.nih.gov/38231165/; main text and supplementary figures assessed for TP53 WT-control assay attribution. The examined evidence was not interpreted as R249S-specific CD4/HLA-II validation.
- IEDB and CEDAR API results were interpreted alongside the original publications. The raw query logs remain in the separate research workspace and are not distributed here.
