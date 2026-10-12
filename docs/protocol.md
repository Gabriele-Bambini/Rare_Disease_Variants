# Evaluation protocol for RQ1 (frozen before any predictor is compared with any label)

Status when written: the label history, the submitter attribution and the AlphaMissense
coordinate join exist. **No code has yet computed any metric that relates a predictor score to a
ClinVar label.** The commit that adds this file precedes the commits that do. Changes after that
point go in the Deviations log at the end, with date and reason, never as silent edits.

Question: after predictors were endorsed for ACMG PP3/BP4, did newly assigned ClinVar labels
become more predictable by the endorsed predictors than by predictors that were not endorsed?

## 1. Data

- ClinVar release of 2026-10-12 and the 40 quarterly `variant_summary` snapshots 2017-01 to 2026-10.
- Universe: GRCh38 missense SNVs from `rdv.clinvar.missense_universe`.
- Scores, all joined on GRCh38 coordinates: AlphaMissense (hg38 file, 2023), REVEL v1.3, ESM1b
  log-likelihood ratio (Brandes et al. 2023), Grantham distance, BLOSUM62.
- Common support: in any comparison, only variants that have a score from every predictor in that
  comparison are used, for all of them. No imputation.

## 2. Cohorts, labels, weights

- Entry: the first snapshot at which a variant is P/LP or B/LB (`first_confident.parquet`).
  Variants already confident in 2017-01 are left out (true date unknown). Label y = 1 for P/LP, 0 for B/LB.
  This is the label the submitter asserted when it entered the confident set, not today's label.
- Submitter: of the SCVs on the variant whose polarity equals the entry label, the one with the
  lowest SCV accession (`first_confident_submitter.parquet`). If none, the stratum is "unknown".
- Why submitters matter: 33 of 5,969 (snapshot, submitter, label) cells hold 68% of all entries,
  and two labs (Labcorp/Invitae, Ambry) supply most of the large ones.
- Weights: each cell (snapshot, submitter, label) gets weight min(1, K / n_cell). Primary K = 500;
  sensitivity K = 250, 1000 and no cap. Weights are deterministic fractions, not random samples.

## 3. Metric: within-gene weighted AUC (WG-AUC)

Strata are genes; only genes with both classes in a window contribute. With P_g, B_g the P/LP and
B/LB variants of gene g in the window and w the weights,

    WG-AUC(s) = sum_g sum_{i in P_g, j in B_g} w_i w_j [1(s_i > s_j) + 0.5 * 1(s_i = s_j)]
                / sum_g (sum_{i in P_g} w_i)(sum_{j in B_g} w_j)

Score orientation (higher means more pathogenic), fixed now: AlphaMissense score, REVEL score,
Grantham distance, minus BLOSUM62, minus ESM1b LLR.

Secondary metric: the unweighted mean over genes with at least 3 variants of each class of the
per-gene AUC.

## 4. Estimand

For a treated predictor T, a control C and a pair of windows, with windows given as ranges of entry snapshots:

    DiD = [WG-AUC(T, post) - WG-AUC(T, pre)] - [WG-AUC(C, post) - WG-AUC(C, pre)]

H1a predicts DiD > 0.

| Test | Milestone | Treated | Control | Pre window | Post window |
|---|---|---|---|---|---|
| **P1** | ClinGen PP3/BP4 calibration (Pejaver et al., Dec 2022) | REVEL | ESM1b | 2021-01 to 2022-10 | 2023-01 to 2024-10 |
| **P2** | AlphaMissense release (Sept 2023) | AlphaMissense | ESM1b | 2022-01 to 2023-07 | 2023-10 to 2025-01 |

ESM1b is the primary control because it is not calibrated for PP3/BP4 before the 2025 ClinGen
extension, and every window above ends before 2025-04. Its skill is close to that of the treated
predictors, which matters: if the label set becomes "easier", strong predictors gain more than weak
ones, and a weak control would then produce a spurious positive DiD.

Secondary (exploratory, no multiplicity claim): controls Grantham and BLOSUM62; REVEL as a stable
comparator in P2; a 2025 test (treated AlphaMissense and ESM1b, control Grantham, pre 2023-10 to
2025-01, post 2025-04 to 2026-10); per-lab DiD for labs with at least 300 within-gene pairs in each
window; only entries with two or more review stars; only entries that came from VUS.

## 5. Event study and placebos

- Event study: WG-AUC per predictor in half-year bins of entry snapshot (Jan+Apr, Jul+Oct), 2019
  onward, with T minus C plotted relative to the milestone bin.
- Placebo: each primary test is repeated with every window shifted 24 and 36 months earlier.
- Pre-trend: slope of (T minus C) over the pre-milestone half-years, with its CI.

## 6. Inference

Cluster bootstrap over genes (resample genes with replacement, recompute every statistic), 2000
replicates, percentile 95% CI, seed 20261012. Holm correction across P1 and P2 only. Report the
bootstrap SE and the minimum detectable effect, 2.8 x SE.

## 7. Decision rule for making RQ1 the headline

- **Go.** For P1 or P2, the 95% CI of the DiD is above 0; both placebo DiDs are below half of the
  primary estimate in absolute value; the sign is unchanged for K = 250, 1000 and no cap, and when
  the two largest labs are excluded.
- **Informative null.** The CI contains 0 and the minimum detectable effect is at most 0.03 AUC.
  Reported as "no detectable feedback, with this power".
- **Underpowered.** Otherwise. RQ1 becomes a descriptive section and the headline moves to the
  circularity decomposition (RQ2 to RQ4).

## 8. Known threats, stated in advance

- Predictor skill and case mix interact (above). Mitigated by the skill-matched control and by placebos.
- Submitter attribution is approximate: a re-evaluated SCV keeps its accession, so it can name the
  lab that later revised the label rather than the one that first made it confident.
- Labs were already using REVEL before 2022 (Color's ClinVar texts from 2020 cite REVEL with
  lab-defined thresholds), so the first milestone may be a formalisation, not an adoption. P2 is
  the cleaner natural experiment, because AlphaMissense did not exist before Sept 2023.
- Predictor training sets overlap ClinVar to different degrees (REVEL's HGMD-era training, no
  clinical labels in ESM1b); the overlap is fixed across windows, so it should difference out.
- DMS reliability bounds any DMS-based check (H1b); see the 2026 MAVE reliability preprint.

## Deviations log

Nothing in sections 1 to 8 was changed after the primary results were seen. What follows is
everything done after the freeze, in order.

**2026-10-12, before any predictor was compared with any label.** `join_scores.py` joined AlphaMissense,
REVEL, ESM1b, Grantham and BLOSUM62. As a correctness check on those joins, whole-universe AUCs
against current labels were computed, outside every protocol window: REVEL 0.968, AlphaMissense
0.961, ESM1b 0.931, Grantham 0.682, BLOSUM62 0.690 (pooled across genes, 208,986 variants with all
five scores). They confirm orientation and identifier mapping and changed no protocol value.

**2026-10-12, after the primary results.** Primary and secondary analyses were run exactly as
specified (`rq1_primary.csv`, `rq1_secondary.csv`, `rq1_event_study.csv`). Then, post hoc:

1. *Window boundary.* An entry at snapshot S was labelled during the quarter before S, so the first
   post-window snapshot straddles each milestone (for AlphaMissense, most of the labels
   in the 2023-10 snapshot were assigned before the 19 September release). This dilutes the post window; it cannot create a null on its
   own. Added a clean-boundary run that starts the post window one snapshot later
   (`rq1_posthoc.csv`). The primary results stay as specified.
2. *Positive control.* Prompted by the null, tested whether the design can see predictor-shaped
   labels at all: labels whose own submission text cites PP3/BP4, or says "in silico", against labels
   whose text does not (`rq1_positive_control.csv`). Exploratory.

Not done yet from this protocol: the secondary macro-average per-gene AUC (section 3), and the
per-lab analyses use the 300-pair threshold only as a reading rule (the table reports the pair counts).
