# Protocol 2 for RQ1: does the score of a VUS predict its promotion, and does that change after a milestone?

**Written after the first RQ1 analysis came back null** ([protocol.md](protocol.md): AUC of newly
labelled variants). That analysis could not tell a missing loop from an insensitive design. This one is
aimed where PP3/BP4 acts: on borderline variants. It was designed after seeing that null, which the paper
must say. Status when frozen: label counts, risk-set covariates and a covariate-only information
calculation exist (`promotion_prefit.py`); **no model relating any score to any promotion has been
fitted.** The commit that adds this file precedes the commit that does.

## 1. Design

Between consecutive ClinVar snapshots p and s, a *case* is a variant that was VUS at p and gets its
first P/LP label at s (direction **up**) or its first B/LB label at s (**down**). Its *controls* are VUS
in the same gene that were still VUS at s and never confident before. Up to 20 controls are drawn at
random from the gene's risk set (seed 20261014). One stratum is one case with its controls. Variants
that went VUS to conflicting or other, were confident before p, or left ClinVar are neither case nor
control: **the test covers promotions of the aggregate label only**, and VUS to conflicting transitions
(3 to 16 times as many as promotions in the big quarters) are outside it.

Covariates, z-scored over the universe and oriented so larger means more pathogenic: REVEL,
AlphaMissense, ESM1b (negated LLR), Grantham, log10 of a frequency proxy (largest of AF_EXAC, AF_TGP,
AF_ESP in the ClinVar VCF; 0 if absent) and an indicator that the proxy exists. For direction down the
first four are negated, so a positive coefficient always means "aligned with the direction of the
promotion". Variants without all four scores are dropped.

Model: weighted conditional logistic regression within strata, covariates plus their interaction with a
post-period indicator, errors clustered by gene (`rdv.clogit`, tested against brute force, finite
differences and a coverage simulation). The estimand is the **change in a covariate's coefficient
between windows, holding the other scores fixed** (a partial effect). This is not diluted by the
correlation between REVEL and ESM1b, which sank the AUC-difference design.

Weights: a case in a (snapshot, submitter, direction) cell of n cases has weight min(1, 500 / n);
submitter is the lowest-accession polarity-matching SCV, as in protocol 1.

## 2. Primary tests (two, both direction up)

| Test | Estimand | Pre events s | Dropped | Post events s |
|---|---|---|---|---|
| **Q1**, Pejaver PP3/BP4 calibration | change in the REVEL coefficient | 2020-04 to 2022-01 | 2022-04 to 2023-01 | 2023-04 to 2025-01 |
| **Q2**, AlphaMissense release | change in the AlphaMissense coefficient | 2022-01 to 2023-07 | 2023-10 | 2024-01 to 2025-01 |

An event at s was labelled during the quarter before s, so intervals containing a milestone are
dropped: the Pejaver preprint (17 Mar 2022) lies in the interval ending 2022-04, the journal version
(20 Oct and 1 Dec 2022) in those ending 2023-01, AlphaMissense (19 Sep 2023) in the one ending 2023-10.
Both post windows end at 2025-01, before the 2025 ClinGen extension. Holm correction over Q1 and Q2.

Why direction up is primary: promotions to benign are mostly decided by population frequency, and
AlphaMissense and REVEL were trained with frequency information that ESM1b lacks. If labs moved to gnomAD
v4 (late 2023) for BA1/BS1/PM2, that would load on AlphaMissense and mimic adoption. A 2015-era
frequency proxy cannot absorb that shift, so the benign direction is secondary.

## 3. Secondary analyses (reported, no multiplicity claim)

Direction down; waves removed (the 10% of intervals with most promotions, defined from counts alone:
2022-07, 2023-04, 2024-04, 2025-04); weight cap 250, 1000 and none; the two largest labs removed;
placebo windows shifted 24 and 36 months earlier (early risk sets are small, so these are weak);
lab-matched controls (the case's own lab must also hold a submission on the control, which isolates one
lab's decision); the change in the other coefficients and in the sum of the three machine-learning
coefficients, as a check on specificity. If the gnomAD slices finish, a green-gene analysis with gnomAD
allele frequency in place of the proxy is added as a further secondary analysis; the primary analysis
does not depend on it.

## 4. Power, computed before any fit

Model-based minimum detectable change (2.8 SE, log-odds per SD, cap 500), from the covariance of the
unpromoted VUS and the case counts only. A lower bound: it ignores gene clustering.

| Test | Direction | Events pre / post | REVEL | AlphaMissense | ESM1b |
|---|---|---|---|---|---|
| Q1 | up | 758 / 1,987 | **0.230** | 0.201 | 0.199 |
| Q1 | down | 1,601 / 8,092 | 0.185 | 0.163 | 0.162 |
| Q2 | up | 1,772 / 915 | 0.226 | **0.197** | 0.195 |
| Q2 | down | 5,647 / 3,344 | 0.225 | 0.194 | 0.191 |

Scale for reading these: PP3 at supporting strength multiplies the odds of pathogenicity by about 2.1
(350^(1/8)), so applied to every promotion across a score threshold it would add roughly 0.4 to 0.7 to a
per-SD slope, and applied to a fraction of promotions proportionally less. A test with a minimum
detectable change near 0.2 can therefore rule out PP3 steering most promotions and nothing subtler.

## 5. Decision rule

- **Go.** For Q1 or Q2, Holm p below 0.05 with a positive change; both placebo changes smaller than
  half of the primary in absolute value; same sign with waves removed, caps 250 and 1000, and the two
  largest labs removed.
- **Bounded null.** Holm p of 0.05 or more and a realized minimum detectable change (2.8 times the
  clustered SE) of at most 0.25. Reported as "no PP3-scale steering of promotions, with this power".
- **Inconclusive.** Anything else. Q1 up may land here a priori: its lower-bound figure is already 0.23.

## 6. Threats stated in advance

- Which VUS get reviewed depends on lab programs, not only on scores; within-gene, within-interval
  strata and the cap remove gene, time and bulk-wave effects but not a lab's choice of which variants to
  look at (the lab-matched analysis is the partial answer).
- Promotion rests on evidence we do not observe (segregation, de novo status, functional assays); the
  test measures score dependence, not mechanism.
- The frequency proxy is coarse (23.6% of the universe has a value).
- Submitter attribution is approximate, as in protocol 1.

## Deviations log

**2026-10-14, performance only.** The lab-matched analysis rebuilt a large set per case and was
killed after several minutes; it now keeps one sorted array per lab. Primary results are byte-identical
before and after (same checksum). No analysis choice changed.

**2026-10-14, after the primary results (post hoc, `promotion_posthoc.py`).**
1. The pre-specified waves-out list (four intervals) removes waves from the Q1 post window but not
   from its pre window, so it compares unlike things (and the benign direction gave +0.51, p = 0.0003).
   Added a symmetric rule from counts alone: remove every interval with more than 5 times the median
   number of promotions (8 intervals) from both windows.
2. Because that result was large, added: all coefficients, the largest case suppliers per period,
   leave-one-lab-out, one lab at a time, and a rolling placebo that slides the same pre / gap / post
   geometry across the whole history.
3. Outcome: no change to any primary or secondary conclusion; the benign-direction lead is recorded in
   `results/README.md` as exploratory. A confirmation on later snapshots must be pre-specified before it
   is run.
