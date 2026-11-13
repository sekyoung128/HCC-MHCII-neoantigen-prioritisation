# Phase 4 — Experiment-informed hypotheses and conditional validation plan (v2)

> **Status and provenance:** This is a **post hoc, hypothesis-generating proposal**, not a wet-lab protocol or record of completed experiments. Original Phase 2 candidate ranks, selection rules, and Phase 3 evidence files remain unchanged. Do not present any result here as verified MHC-II presentation, HLA restriction, or CD4 immunogenicity. This v2 document revises the explanatory decision logic only; original generated Phase 4 files should be retained for audit.

## 1. Scientific objective and scope

For HCC-associated **TP53 R249S** and **CTNNB1 S45P**, separate five biological endpoints: (i) algorithmic peptide–HLA-II binding-related ranks, (ii) directly observed peptide binding, (iii) functional CD4 recognition, (iv) experimentally demonstrated allele-specific HLA-II restriction, and (v) processing-dependent presentation. The aim is to propose experiments that could distinguish these possibilities, **not** to claim validated neoantigens or that a selection strategy is superior.

### Frozen predicted inputs (unchanged)

| Case | Origin | WT → MUT sequence | HLA | BA WT → MUT | EL WT → MUT | Core WT → MUT |
|---|---|---|---|---:|---:|---|
| A | A — BA-first (Phase 2 rank 1) | `RRPILTIITLEDSSG` → `RSPILTIITLEDSSG` | DRB1_0401 | 3.08 → 4.35 | 2.46 → 2.95 | `ILTIITLED` → `ILTIITLED` |
| B | B — Mutation-aware (Phase 2 rank 1) | `CMGGMNRRPILTIIT` → `CMGGMNRSPILTIIT` | DRB1_0101 | 32.79 → 31.87 | 89.94 → 66.84 | `MNRRPILTI` → `MNRSPILTI` |
| C | C — Experiment-informative (Phase 2 rank 1) | `MGGMNRRPILTIITL` → `MGGMNRSPILTIITL` | DRB1_0401 | 30.24 → 35.01 | 87.59 → 64.24 | `RRPILTIIT` → `MNRSPILTI` |
| D | Exploratory CTNNB1 coverage (not Phase 2 Top-K) | `LDSGIHSGATTTAPS` → `LDSGIHSGATTTAPP` | DRB1_0401 | 41.23 → 46.85 | 3.53 → 3.14 | `IHSGATTTA` → `IHSGATTTA` |

Lower `%Rank_BA` and `%Rank_EL` indicate more favourable **model-specific** percentile scores; BA and EL optimise different endpoints. `ΔBA = WT %Rank_BA − MUT %Rank_BA` (positive favours MUT within that endpoint). It is **not** a physical binding-affinity difference or evidence of mutant-selective CD4 recognition. A, B and C are locked Phase 2 top-ranked representatives. D was added post hoc for CTNNB1 coverage and is **not** part of the locked Top-K comparison.

## 2. Candidate-specific hypotheses and differential experiments

### Case A — BA-first / TP53 R249S, start 248, DRB1*04:01

**Case-specific hypothesis (proposed):** Primary hypothesis: measurable binding to both WT and MUT can coexist with an absence of mutant binding advantage; binding cannot substitute for recognition.

**Frozen comparison:** WT RRPILTIITLEDSSG vs MUT RSPILTIITLEDSSG; residue 249 is at position 2, outside identical predicted core ILTIITLED.

**Initial experiment and branching logic:** First run paired DRB1*04:01 WT/MUT binding. If technically valid binding exists, assess CD4 mutant discrimination with WT counter-stimulation. Escalate to allele-dependence only for a reproducible CD4 response, and to full-length studies only with validated responder and APC systems.

**Evidence that would disfavor the working hypothesis / stop rule:** If both peptides fail binding after valid assay QC, deprioritise A. Shared WT/MUT CD4 recognition would contradict a mutant-specific recognition claim in that tested responder, regardless of binding.

**Feasibility prerequisite:** DRB1*04:01 binding reagents or validated DR4 monoallelic binding system, WT/MUT synthetic peptides, allele-specific reference ligand, quantitative SOP. Replication count, positivity cut-off, minimum detectable difference and uncertainty conventions must be fixed from the validated assay SOP and controls before any outcome is reviewed.

**Interpretation limit:** Mutation is outside the predicted 9-mer core; shared core is a prediction, not proof of identical presentation or T-cell recognition.

### Case B — Mutation-aware / TP53 R249S, start 242, DRB1*01:01

**Case-specific hypothesis (proposed):** Primary hypothesis: a small percentile rank improvement need not yield reproducible biochemical binding differences, particularly with weak absolute predicted rank.

**Frozen comparison:** WT CMGGMNRRPILTIIT vs MUT CMGGMNRSPILTIIT; mutation is within the predicted binding core, but relative percentile changes are not physical affinity differences.

**Initial experiment and branching logic:** Begin with quantitative paired DRB1*01:01 binding; compare WT/MUT differences against measured assay variability, rather than assigning physical significance to +0.92 BA-rank units. Functional work is conditional on binding and reagents.

**Evidence that would disfavor the working hypothesis / stop rule:** No reliable binding or no distinguishable WT/MUT difference after QC stops expensive B follow-up; that outcome cannot prove all R249S windows are irrelevant.

**Feasibility prerequisite:** DRB1*01:01 binding system, paired intact peptides, validated binding range and reference controls. Replication count, positivity cut-off, minimum detectable difference and uncertainty conventions must be fixed from the validated assay SOP and controls before any outcome is reviewed.

**Interpretation limit:** Low absolute BA/EL predictions; positive WT–MUT rank delta is not evidence of mutant specificity.

### Case C — Experiment-informative / TP53 R249S, start 243, DRB1*04:01

**Case-specific hypothesis (proposed):** Primary hypothesis: a model-predicted register switch is investigable only if binding is experimentally supported; BA and EL directions do not by themselves establish contradictory biology.

**Frozen comparison:** WT MGGMNRRPILTIITL vs MUT MGGMNRSPILTIITL; predicted cores RRPILTIIT vs MNRSPILTI; MUT BA %Rank 35.01, EL %Rank 64.24.

**Initial experiment and branching logic:** First establish paired binding. If repeatable binding supports the case and a register-resolving platform is available, prospectively define a separate exploratory register-mapping comparison for predicted WT RRPILTIIT versus MUT MNRSPILTI, with controls addressing artefacts caused by truncation. The four main evidence gates do not themselves prove register usage.

**Evidence that would disfavor the working hypothesis / stop rule:** If reliable binding is absent, stop register mapping and do not claim the predicted switch is biological. If binding is present but register mapping ambiguous, leave register attribution unresolved.

**Feasibility prerequisite:** DRB1*04:01 binding assay, valid comparator ligands, access to register-resolving methodology before planning specialised follow-up. Replication count, positivity cut-off, minimum detectable difference and uncertainty conventions must be fixed from the validated assay SOP and controls before any outcome is reviewed.

**Interpretation limit:** BA and EL optimise different endpoints; opposing shifts do not prove biological model conflict. Predicted core shift is not experimentally established.

### Case D — Post hoc exploratory / CTNNB1 S45P, 15-mer terminal residue, DRB1*04:01

**Case-specific hypothesis (proposed):** Primary hypothesis: the selected 15-mer may have window-dependent behaviour because the mutation is the terminal residue; low EL rank alone does not establish processing.

**Frozen comparison:** WT LDSGIHSGATTTAPS vs MUT LDSGIHSGATTTAPP; substitution is residue 15/15, outside identical predicted core IHSGATTTA.

**Initial experiment and branching logic:** First assay the exact terminal-substitution 15-mer, acknowledging peptide identity and C-terminal integrity. If exploring a mutation-centred window, predefine WT/MUT sequence pairs and HLA input separately and label them new exploratory candidates.

**Evidence that would disfavor the working hypothesis / stop rule:** Failure to bind for the tested terminal 15-mer cannot be extrapolated to all S45P-containing windows; no post hoc window can be inserted into frozen Phase 2 rankings.

**Feasibility prerequisite:** DRB1*04:01 binding reagents; confirmed peptide C-terminal integrity; ability to predefine alternative windows only if justified. Replication count, positivity cut-off, minimum detectable difference and uncertainty conventions must be fixed from the validated assay SOP and controls before any outcome is reviewed.

**Interpretation limit:** Mutation at peptide C terminus may make conclusions dependent on selected window. Any centred or extended window is a new prelabelled exploratory analysis, not a Phase 2 rerank.

## 3. Sequential evidence gates — not a fixed four-experiment commitment

| Gate | Primary experimental endpoint | Essential control | Prerequisite to progress | What remains unresolved |
|---|---|---|---|---|
| 1. HLA-II binding | Direct or competitive binding of the exact WT/MUT peptide to the chosen allele | Known allele-specific binder, non-binder, background, peptide integrity, assay performance | Reproducible technically valid result; measure uncertainty | CD4 recognition, restriction, natural presentation |
| 2. Functional CD4 recognition | WT versus MUT stimulation using **CD4-defined** cells and antigen-presenting cells | No peptide, irrelevant peptide, matched WT, functional positive control, donor/clonal context, APC viability/HLA expression | Reproducible positive CD4 signal before restriction work | Allele-specific restriction; failure in one repertoire is not a general negative |
| 3. Allele-specific restriction | Demonstrate dependence on the named DRB1 allele | Allele-defined APC or controlled allele perturbation; expression/viability/isotype controls | Positive functional CD4 response and appropriate allele-specific tools | Broad anti-HLA-DR blockade alone cannot assign DRB1*01:01 or *04:01 |
| 4. Processing-dependent recognition | Compare full-length antigen to peptide-pulsed rescue in defined processing-capable system | Full-length WT/MUT expression, antigen expression, HLA-II competence, pulsed positive control, cell viability | Verified responding CD4 system and upstream evidence | Expression in APC model does not establish presentation by patient HCC tumours |

**Order is conditional:** Gate 1 precedes costly downstream work for most candidates. Functional (gate 2) and restriction (gate 3) require a suitable, competent CD4 system, which is **not currently established**. A failed early gate can terminate a branch without implying that the mutation is universally non-immunogenic. The matrix lists all 16 *possible* candidate-gate combinations; it does not prescribe executing all 16 experiments.

## 4. Additional planned mechanistic comparison for Case C (conditional, outside the 16-gate matrix)

- **When:** Only after reliable allele-matched binding is demonstrated for at least one C peptide and specialised experimental resources are available.
- **Question:** Does the experiment support different binding registers for the exact C WT and MUT peptides?
- **Method consideration:** A prospectively defined register-sensitive analysis may use a suitably controlled series of truncated/overlapping peptides or a more direct register-resolving approach, recognising that truncation itself changes HLA-II binding.
- **Interpretation:** Distinct predicted 9-mer cores are not proof of a register switch. General biochemical binding and CD4 assay outcomes alone cannot locate the register. A non-resolving result should remain **unresolved**, not be converted into a mechanism claim.

## 5. Explicit stop/go framework and required pre-specification

**Before any wet-lab work**, record: exact WT/MUT sequences and chemical identity; target HLA allele and allele-defined reagents; validated assay SOP; positive/negative control performance; concentration range; positivity/background definition; repeatability and uncertainty; independent biological versus technical replicates; donor/clone description (if applicable); predefined handling of an inconclusive result; and resource availability. **Do not invent numeric thresholds** in this desk-based plan.

| Case | Primary decision point | GO to next stage only if | STOP / defer scenario |
|---|---|---|---|
| A | Pairwise binding and subsequent CD4 selectivity | Valid binding then reproducible, interpretable functional readout with feasible APC/CD4 resources | No reproducible binding after QC; WT cross-reactivity precludes mutant-selectivity claim |
| B | Is +0.92 rank change experimentally measurable? | Binding and WT–MUT difference are reproducible relative to assay variability | Weak/absent binding; no convincing difference; inadequate functional resources |
| C | Is any binding demonstrable before register follow-up? | Binding with appropriate quality controls and specific register-testing capability | No binding, unresolved register, or impractical specialised resources |
| D | Is terminal-substitution 15-mer interpretable? | The tested exact 15-mer binds reproducibly and suitable WT/MUT functional model exists | No binding for this window; additional windows would require new, separately labelled exploratory design |

**A negative assay is endpoint-limited.** Lack of peptide binding in a defined in vitro system cannot prove absence of CD4 recognition in every HLA or repertoire; absence of peptide-pulsed CD4 response cannot exclude an alternative T-cell repertoire; absence of full-length recognition does not establish failed presentation in all tumours.

## 6. Evidence provenance and limitations

- **Phase 1–2:** 216 overlapping peptide-window–HLA comparisons from just **two mutations**; non-independent windows must not be counted as distinct tumour mutations. Candidate selection was rule-based after a pilot. The three prioritisation schemes optimise different goals; Top-K divergence is not experimental superiority.
- **Phase 3 IEDB/CEDAR:** Six T-cell assay IDs all trace to **WT control peptides** used in TP53 R248Q/R248W experiments, not the R249S mutant peptide. IEDB/CD4 annotation, allele-specific restriction and direct binding for the shortlisted mutant candidates remain unsupported by these search results.
- **Historical Repáraz comparison:** A separate selected HCC cohort of 13 peptide cases, including 10 experimental HOM2 binding-score-positive and 3 non-binders. Selection using an earlier NetMHCIIpan and incomplete model training-data overlap provenance prevent claims of independent validation or general predictive accuracy; this cohort does not directly validate R249S or S45P.
- **Case independence and proposal status:** A, B and C are overlapping windows of the same R249S mutation, **not independent biological discoveries**; D is explicitly post hoc. No proposed experiment has been performed. Donor HLA, clinical tumour processing, CD4 repertoire and resource access are unverified.
- **Cancer specificity:** Responses in a peptide-pulsed or engineered APC model do not establish presentation in a patient's HCC tumour, and WT cross-recognition is a separate safety-relevance uncertainty.

## 7. Decision on project continuation

**GO — finish as an independent, reproducible computational research portfolio**: frozen inference, transparent ranking, source-level evidence audit, bounded experimental literature comparison, and falsifiable *conditional* experiment plans. **NO-GO — infer a novel validated MHC-II neoantigen, superior A/B/C strategy, or out-of-sample model performance.** Further model/peptide expansion is justified only if independently attributable experimental peptide–HLA/CD4 outcomes and a leakage-aware evaluation design become available.

## 8. Audit trail

- Original `phase4_experiment_plan.md` and `phase4_experiment_decision_matrix.csv` should be retained unchanged, alongside original Phase 2 and Phase 3 files.
- Revised outputs: `phase4_experiment_plan_v2.md`, `phase4_experiment_decision_matrix_v2.csv`.
- Source Phase 4 document SHA-256: `33561a63a16c74560d50265c0a55070785ecf5be65815ec24478cb8c1cfb76dc`.
- Source decision matrix SHA-256: `958599fe420e02fd35cbc905d8a23199ad7285014ee0b742213da0a9765f292f`.
- Source material: user-supplied Phase 1–3 analysis outputs, locked Phase 2 ranking, Phase 4 original drafts, and previously reviewed Repáraz (2022) and IEDB/CEDAR evidence. No new empirical results or external citations have been added to this revision.
