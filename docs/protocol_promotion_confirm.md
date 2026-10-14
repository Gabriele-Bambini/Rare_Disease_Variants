# Protocol 3: confirming the Ambry benign-direction lead on intervals not used to find it

Context: [protocol_promotion.md](protocol_promotion.md) found, after the fact, that benign promotions
by one lab (Ambry), with bulk quarters removed, align much more with REVEL after 2022 than before
(+1.88 per SD, SE 0.39), while all other labs together are flat. Two checks on the same data held
(AlphaMissense does not move with it; the effect survives controls restricted to Ambry's own catalogue).
Status when frozen: nothing here has been fitted. Counts and a power figure were computed from risk sets
alone (`promotion_confirm_prefit.py`).

## What is tested, and what is not

**Hypothesis (one-sided):** in benign promotions from non-bulk intervals 2025-07, 2025-10, 2026-01,
2026-07 and 2026-10, the REVEL coefficient, with AlphaMissense, ESM1b, Grantham and the frequency proxy
held fixed, is larger for Ambry's cases than for other labs' cases.

This compares Ambry with other labs *in the same intervals*. It does not compare these intervals with
Ambry's earlier, 85-case baseline, which would carry regression to the mean. So a confirmation says
"Ambry's benign calls track REVEL more than other labs' do", **not** "Ambry changed in 2022": the rolling
placebo already weakened the timing claim and nothing here rescues it. It also does not separate a lab
that steers by REVEL from one whose decisions happen to follow what REVEL captures beyond AlphaMissense
and ESM1b; only an external anchor (deep mutational scans) can.

Disclosure: these intervals were inside the pooled (all labs together) rolling-placebo summary, but no
Ambry-specific estimate was computed for them before this protocol.

## Design (identical to the discovery, nothing re-tuned)

Benign direction; bulk intervals removed by the rule "more than 5 times the median number of promotions"
(this removes 2025-04 and 2026-04); weight cap 500; covariates, 20 controls per case, seed and submitter
attribution as in protocol 2. One model with an Ambry indicator: coefficients for the other labs, and
Ambry's difference as the interaction; errors clustered by gene.

- **Primary:** lab-matched controls (the control must be a variant on which the case's own lab holds a
  submission). Cases with no attributed lab are dropped from this version.
- **Secondary:** controls from all VUS in the gene.
- Specificity, reported for both: the same difference for AlphaMissense, ESM1b and Grantham.

## Power, computed before any fit

Benign promotions in these intervals: **Ambry 775, all other labs 191.** Model-based minimum detectable
lab difference (2.8 SE, per SD, lower bound): **REVEL 0.43**, AlphaMissense 0.37, ESM1b 0.36. The discovered
gap was about 1.6 (Ambry +1.67 after, others +0.09). The rule for running was a figure no larger than half
the discovered gap (0.8): met.

## Decision rule

- **Confirmed:** one-sided p below 0.05 for the REVEL difference in the primary analysis, and a positive
  difference in the secondary.
- **Not confirmed:** otherwise. The lead is dropped and reported as an exploratory result that did not replicate.
- Judged on sign and significance, not on reproducing +1.88. If the AlphaMissense difference is as large as
  REVEL's, the reading is general alignment with machine-learning scores, not REVEL.
- **Stop rule:** the VUS-promotion line of work (RQ1, option B) ends after this analysis, whatever the result.

## Deviations log

None to the analysis. **Result (2026-10-14): not confirmed** (REVEL difference +0.14, SE 0.61, one-sided
p = 0.41; secondary +0.32, SE 0.58, p = 0.29). The outcome-blind power figure above (0.43) was a lower
bound that ignored gene clustering; the realized SE was about four times the one it implied, so the run
was less able to detect a difference than the "runnable" criterion suggested. Reported in `results/README.md`.
