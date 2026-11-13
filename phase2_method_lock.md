# Phase 2 — Post-pilot method lock (11 November 2026)

## Scientific question
How does the choice of selection rule affect which HCC-associated MHC-II peptide–HLA candidates are prioritised? The three rules focus on absolute binding prediction, WT-to-mutant rank changes and questions for experimental follow-up.

## Status and provenance
These rules were fixed after reviewing the Phase 1 pilot and extension summaries, but **before calculating the Phase 2 rankings**. The comparison is exploratory, not prospectively preregistered. Experimental outcomes were neither available to nor used by the ranking algorithms. This avoids using observed assay labels to choose among the candidates, but it does not rule out overlap with model training data, earlier scientific knowledge or decisions influenced by the pilot.

Input: exact byte-level SHA-256 `7b1cebf1d4c282a1041c7c94bd15392ee99cc8f78de9fd8c7d716415d6415b6f` for `wt_mut_comparison_216.csv` (216 matched WT/MUT peptide-window–HLA units). Gene/variant: CTNNB1 S45P; TP53 R249S. Alleles: DRB1_0101 and DRB1_0401 (hypothetical); peptide lengths 15, 18 and 21 aa. Only **two distinct mutations** are represented; overlapping windows and HLA settings must not be counted as independent biological cases.

## Definitions
For each window–HLA candidate, `deltaBA = Rank_BA_WT - Rank_BA_MUT`, `deltaEL = Rank_EL_WT - Rank_EL_MUT`. A lower percentile rank is more favourable within a model endpoint. A positive delta indicates a lower mutant rank, not measured binding or mutant-specific T-cell recognition.

## Three rules; identical 216-row pool
**A — BA-first:** ascending mutant `Rank_BA_MUT`; ties ascending mutant `Rank_EL_MUT`; further ties alphabetical `candidate_id = pair_id + "|" + MHC`.

**B — mutation-aware (BA):** prefer `deltaBA>0` before `deltaBA<=0`; within each group ascending mutant `Rank_BA_MUT`; ties descending `deltaBA`; then candidate_id. This ordering is a fixed decision rule, **not** a validated neoantigen score. A positive shift does not establish binding when the absolute mutant rank remains poor.

**C — experiment-informative:** three pre-specified queues and a deterministic round-robin, without replacement: (1) opposite-sign nonzero deltaBA and deltaEL, sorted decreasing `|deltaBA-deltaEL|`, then increasing `min(mutant BA rank, mutant EL rank)`, then candidate_id; (2) changed predicted WT/MUT core offset, sorted `mutation inside MUT core` first, increasing mutant BA rank, decreasing `|deltaBA-deltaEL|`, then candidate_id; (3) HLA sensitivity for the same window, with opposite-sign BA delta across HLA first, then decreasing absolute between-HLA mutant BA-rank difference, then increasing mutant BA rank, then candidate_id. Select one unused candidate from each queue cyclically, skipping already selected records. When no queue can supply more, fill from all remaining by ascending mutant EL rank then candidate_id. The selection reason is recorded. **C is not an immunogenicity or vaccine-suitability score.** BA and EL represent different prediction endpoints. Choosing candidates because those outputs differ does not establish that either is wrong.

## Comparisons and reporting
Use global top-5, top-10 and top-20 without imposing equal quotas per gene or HLA. Report exact top-K overlap, Jaccard intersection/union, all-candidate Spearman rank correlations (descriptive), per-strategy composition by gene/HLA, number of unique predicted mutant 9-aa cores, and proportion with positive deltaBA or predicted mutation-in-core. For core redundancy, report number of unique `(gene, HLA, Core_MUT)` groups among top-10; this is **only a heuristic**, not proof of independent biological epitopes. Report absolute mutant ranks and mutation-associated changes **separately**.

The analysis neither trains nor validates a predictive model. It cannot establish CD4 recognition, tumour presentation, allele-specific restriction, vaccine efficacy or independent predictive performance. Because A, B and C answer different selection questions, differences between their rankings should not be presented as biological superiority.

## Software and reproducibility
Run `python phase2_prioritise.py` from the project root after placing `phase2_rules.json` and the script there. Inputs should be in `results/phase1_extension/analysis/wt_mut_comparison_216.csv`. The script checks the exact input SHA-256 and exports `results/phase2/` outputs and hashes. If a rerun produces different CSV bytes, investigate the cause rather than changing the frozen inputs or rules to force a match.
