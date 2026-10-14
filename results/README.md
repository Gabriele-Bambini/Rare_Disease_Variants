# Week 2: ClinVar label history and what labs write

Data snapshot: ClinVar current release of 2026-10-12 and the quarterly `variant_summary`
archives 2017-01 … 2026-10. Missense universe: GRCh38 single-nucleotide variants whose HGVS
protein change is one amino-acid substitution (start-loss and stop changes excluded).

Everything here is produced by code in this repository:

```bash
pip install -e '.[dev]'
export DATA_DIR=...                       # see scripts/hpc/README.md
bash scripts/hpc/download_all.sh clinvar_current
CLINVAR_ARCHIVE_KINDS=variant_summary bash scripts/hpc/download_all.sh clinvar_archives
python scripts/clinvar/build_history.py --data-dir $DATA_DIR   # universe + 40 snapshots, ~6 min on 4 cores
python scripts/clinvar/scv_mentions.py  --data-dir $DATA_DIR   # SCV text analysis, ~1 min
python scripts/clinvar/label_flow.py    --data-dir $DATA_DIR   # label flow + summary.csv
pytest
```

| File | Content |
|---|---|
| `summary.csv` | headline counts quoted below |
| `label_flow_by_quarter.csv`, `fig_label_flow.png` | variants first seen with a confident label, per snapshot, by route |
| `scv_mentions_{labs,all}.csv` | share of text-bearing SCVs naming each term, per quarter, with Wilson 95% CI, number of submitters and top-submitter share |
| `scv_standardized.csv` | the same share with every lab held at its 2021-22 volume weight |
| `scv_decomposition.csv` | pooled change split into within-lab, between-lab and entry/exit parts, for two milestones |
| `scv_by_lab.csv`, `scv_composition.csv` | the 15 largest labs per year; lab shares of volume per quarter |
| `fig_scv_mentions.png` | the H1c figure |

## 1. The label stream (`fig_label_flow.png`)

- 2,552,260 missense variants. Today: 70,172 P/LP, 145,731 B/LB, 2,155,628 VUS, 120,436
  conflicting, 60,293 other.
- 262,711 were P/LP or B/LB in at least one snapshot. 38,420 of them already were in 2017-01
  (left-censored: true date unknown). The other 224,291 became confident later: 195,216 entered
  ClinVar already labelled, 23,812 were resolved from VUS, 5,263 from a conflict.
- Only 7,395 (3.3%) of those 224,291 had two or more review stars when first confident.
  Robustness checks restricted to 2-star or expert-panel labels will have little power.
- Direct P/LP to B/LB flips are rare: 297 in ten years. Most transitions go through
  "conflicting" (B/LB to conflicting 33,133; P/LP to conflicting 13,643).
- **The stream is lumpy.** Single quarters carry bulk batches: 24,249 new B/LB in 2023-04
  (19,172 of them entering already labelled), 20,465 in 2020-01, 16,316 in 2024-04, while most
  quarters add well under 6,000. P/LP has a smaller burst every April (3,699 in 2022-04,
  4,292 in 2023-04, 4,103 in 2024-04). Any quarter-cohort comparison is dominated by whoever
  submitted in bulk that quarter unless the design accounts for the submitter.

## 2. What labs write (`fig_scv_mentions.png`)

3,731,912 SCVs sit on missense variants; 2,807,323 come from clinical-testing labs with
assertion criteria and carry at least 30 characters of text. Time axis: `DateLastEvaluated`.

- **The placebo passes.** AlphaMissense is named in 0 of 963,912 texts evaluated before its
  2023-09 release.
- **AlphaMissense is almost never named afterwards:** 93 of 1,966,591 texts from 2023 to 2025
  (0.005%).
- **The largest labs do not name tools.** Ambry supplies 1,582,112 SCVs (42.4% of all missense
  SCVs, 56.4% of the text-bearing lab SCVs) and names REVEL or PP3/BP4 in none of them. It does
  write "predicted to be deleterious by in silico analysis" (19.6% of its texts), so computational
  evidence is used but unattributed. Labcorp (Invitae) and GeneDx also name neither.
- **REVEL mentions come from a few labs.** 82.4% of the 41,580 post-2023 mentions come from
  three labs (All of Us, Color, 3billion). 165 of 711 labs ever name REVEL.
- **The pooled series mostly tracks lab mix.** Change from 2021-01..2022-09 to 2023-01..2025-12,
  in percentage points, over the 36 labs with at least 50 texts in both windows (98 to 99% of
  volume), plus entry/exit of the remaining labs:

| Term | Pooled before | Pooled after | Pooled change | Within labs | Mix | Entry / exit |
|---|---|---|---|---|---|---|
| PP3 / BP4 cited | 0.75 | 0.43 | -0.33 | -0.30 | -0.13 | +0.10 |
| REVEL named | 1.09 | 2.11 | +1.02 | -0.39 | +0.48 | +0.93 |
| "in silico" (any tool) | 19.19 | 20.75 | +1.56 | -4.44 | +4.97 | +1.03 |
| SIFT / PolyPhen etc. (control) | 38.13 | 6.90 | -31.24 | -15.41 | -15.97 | +0.14 |
| Conservation (control) | 31.13 | 25.29 | -5.84 | +3.54 | -9.06 | -0.33 |
| Population frequency (control) | 58.36 | 34.14 | -24.22 | -1.38 | -23.62 | +0.78 |

The pooled REVEL share doubles after the first milestone, yet labs that kept submitting
mention it *less* (-0.39), and the entry term is largely All of Us joining with REVEL in about
94% of its texts (16,824 mentions after 2023 are 0.86 points of the 0.93). The control
"SIFT/PolyPhen" falls by half within labs (Labcorp: about 71% of its texts in 2022, about 34% from
2023) and by half through the mix (volume moved to labs that never name it).

## What this means for the design

1. **H1c as specified is not usable as evidence of the feedback loop.** Free text measures how
   labs *write*, and the labs that supply most labels do not name their tools. A flat or falling
   mention rate says nothing about use. Keep the figure as a descriptive result with the
   composition panel attached, and put the weight of RQ1 on the label-based tests (H1a, H1b).
2. **The difference-in-differences needs a submitter-aware design.** Bulk benign batches and
   annual P/LP bursts will otherwise drive any cohort contrast. Next step: attach the first
   confident SCV's submitter to each variant, then cluster, cap or exclude by submitter-quarter.
3. **Power is limited for the strict-evidence strata** (2-star and expert-panel), see section 1.
4. **Text measures need a lab-weighted or lab-fixed-effects form** whenever they are used.
   `rdv.scv_text.decompose` and `standardized_share` do this.

Known limits: the current `submission_summary` keeps only each SCV's latest text and date, so a
re-evaluated SCV moves to its newest quarter; variants withdrawn from ClinVar are invisible;
snapshot resolution is one quarter.

# Week 3 to 5: RQ1, did endorsed predictors gain on new labels? (`fig_rq1.png`)

Analysis exactly as frozen in [docs/protocol.md](../docs/protocol.md) (commit `ec2fd38`), with the
post-hoc additions listed in its Deviations log. Cohort: 224,291 variants that entered ClinVar as
P/LP or B/LB after 2017-01. Scores joined on GRCh38 coordinates: AlphaMissense 97.6% of the
universe, REVEL 99.1%, ESM1b 97.5% (1,073 wild-type residue mismatches, left unscored). Metric:
within-gene weighted AUC with lab-quarter weights (cap 500), gene-cluster bootstrap with 2,000 replicates.

## Primary tests

| | Treated vs control | Difference-in-differences | 95% CI | Detectable effect | Holm p |
|---|---|---|---|---|---|
| P1, Dec 2022 | REVEL vs ESM1b | -0.0100 | -0.034 to +0.010 | 0.031 | 0.72 |
| P2, Sep 2023 | AlphaMissense vs ESM1b | -0.0004 | -0.019 to +0.018 | 0.026 | 0.93 |

Within-gene AUC, before to after: REVEL 0.912 to 0.894, ESM1b 0.888 to 0.881 (P1); AlphaMissense
0.918 to 0.921, ESM1b 0.878 to 0.881 (P2). Nothing changes the sign or the conclusion when the cap is
250, 1000 or absent, when the two largest labs are removed, or when the post window starts one
snapshot later (clean boundary: -0.009 and -0.005). **Decision rule of the protocol:** neither test
reaches Go. P2 is an informative null (detectable effect 0.026, below 0.03); P1 is underpowered by the
letter of the rule (0.031).

## How much to trust the null

- **The design's noise floor is above its bootstrap error.** Placebo windows give -0.008 and -0.025
  (P1; the second has a 95% interval that excludes zero) and -0.001 and -0.005 (P2). Variation across
  windows from lab and batch composition is as large as either primary estimate, so effects under
  about 0.02 to 0.03 are out of reach.
- **We have not shown that the design can see predictor use.** Among labels entering from 2019 on,
  those whose own submission text cites PP3/BP4 (8,974, against 100,276 with other text) show no
  larger REVEL-over-ESM1b advantage than the others: -0.003 (-0.032 to +0.022). The same for "in
  silico" wording (72,517 against 31,675): +0.008 (-0.024 to +0.036). This is not a clean positive
  control: PP3 can rest on BayesDel, MutPred2, VEST4 or AlphaMissense rather than REVEL, so a flat
  REVEL gap there does not prove the design is blind. It shows only that this particular check did not
  light up. Case mix moves every predictor together ("in silico" labels are easier for all of them,
  Grantham included, +0.03), which is why a difference between predictors cancels it. A real positive
  control needs labels known to be assigned with REVEL; only 9 of 201,081 labels are first asserted by
  a submission that names REVEL, because the REVEL-citing labs (All of Us, Color, 3billion) are rarely
  the first to label a variant.
- **A correlated control absorbs signal.** Within a gene, REVEL and ESM1b ranks correlate at 0.57
  (REVEL and AlphaMissense 0.66, AlphaMissense and ESM1b 0.67). If labs steered labels along REVEL,
  ESM1b would partly move with them, so a differential recovers well under half of such a shift.
  This is a back-of-envelope argument, not something measured here.
- **The only positive-looking signal is not usable.** AlphaMissense minus ESM1b is 0.125 (0.060 to
  0.173) in the 2023H2 bin, but that bin holds 770 weighted pairs, and it covers labels that entered
  between April and September 2023 (bins are named by snapshot month: H1 is the January and April
  snapshots, H2 the July and October ones), nearly all before the 19 September release. The same
  difference sits at 0.02 to 0.06 in every other bin, including years when AlphaMissense did not exist.
- **Subgroup results are noise-level.** 3 of the 19 secondary intervals exclude zero (VUS-resolved
  entries for P1, -0.059; two labs with under 120 pairs in one window), about what chance gives.

## What this means

No differential gain of an endorsed predictor on new labels is detectable at about 0.02 to 0.03 AUC. That
is not evidence that no loop exists: the design has no validated positive control, and the loop may act
mainly on borderline variants (which VUS get promoted), not on how already-confident labels are ordered.
Options, with the first three combinable:

0. Test promotion directly, which uses the large denominator the AUC design discards: among the ~2.1M
   VUS in each snapshot, after Dec 2022 do high-REVEL VUS leave VUS toward P/LP faster than
   comparable high-ESM1b VUS (and the same for benign)? About 24,000 variants were resolved from VUS
   in ten years. Needs its own pre-specified protocol, written before it is run.

1. Replace the AUC difference with a model of the label that conditions on the control: within gene,
   does the label depend on REVEL beyond ESM1b, and does that partial dependence change after the
   milestone? This uses continuous variation and is not diluted by the control's correlation with the
   treated predictor. It needs its own pre-specified protocol, written before it is run.
2. Run the simulation of section 3.4 of the plan to show what size of loop this design can and cannot
   detect, and report the null as a bound.
3. Make the circularity decomposition (RQ2 to RQ4) the headline and keep RQ1 as a bounded, honest
   negative result, as the plan's section 3.5 allows.

# Protocol 2: does a VUS's score predict its promotion? ([protocol](../docs/protocol_promotion.md))

Planned analysis, run as frozen (commit `2a75dae`). Cases: 758 / 1,987 (Q1) and 1,772 / 915 (Q2) VUS
promoted to P/LP before / after the milestone, each with up to 20 unpromoted VUS from the same gene
and interval. Change in the coefficient per SD, with the other scores held fixed.

| Primary test (promotion to P/LP) | Change | 95% CI | Holm p | Realized detectable change |
|---|---|---|---|---|
| Q1, REVEL around Dec 2022 | -0.030 | -0.232 to +0.171 | 0.77 | 0.288 |
| Q2, AlphaMissense around Sept 2023 | +0.081 | -0.099 to +0.262 | 0.75 | 0.258 |

**By the rule written beforehand, both are inconclusive**, not a bounded null: the test could only have
seen changes above 0.26 to 0.29 (threshold 0.25). It rules out PP3 steering most promotions and nothing
subtler. Caps of 250 and 1,000, no cap, removing the two largest labs, lab-matched controls and both
placebo shifts all agree (all 28 secondary results are nulls except one, below). Level fact: given each
other and ESM1b, REVEL and AlphaMissense each carry a partial effect of about +0.65 per SD on promotion
to P/LP, **before and after** (REVEL 0.65 to 0.62, AlphaMissense 0.62 to 0.70). A loop that was already
there in 2020 cannot show up as a change.

## One result that was not planned, and how far to trust it

Benign direction, with the bulk-submission quarters removed from both windows (post hoc), around Dec
2022: the REVEL coefficient goes from +0.03 to +0.74 per SD (change +0.71, SE 0.19). Nearly all of it
sits in one lab, found after the fact:

| Benign promotions, bulk quarters removed | REVEL, before | after | change (SE) | cases |
|---|---|---|---|---|
| Ambry only | -0.21 | +1.67 | **+1.88 (0.39)** | 85 / 223 |
| Everyone except Ambry | +0.13 | +0.09 | -0.04 (0.21) | 300 / 173 |

AlphaMissense, ESM1b, Grantham and the frequency proxy do not move in the same fit, so it is specific to
REVEL, which is what steering by BP4 would look like. Reasons for caution, all real:

- It was found by looking: direction, bulk-quarter rule and lab were all chosen after seeing the data.
- Slid across every position in the history with the same geometry, the benign-direction change in REVEL
  has median +0.27 and exceeds |z| = 2 in 4 of 20 positions, so it is not a sharp break at the milestone.
- The benign direction is the one population frequency dominates, and REVEL carries frequency information.
- The other lab with enough cases (GeneDx) points the same way but is imprecise (+0.57, SE 0.45, 27 cases after).

Status: a lead, not a result. It can be tested on snapshots not used in choosing it (2025-04 onward).

## Protocol 3: the Ambry lead, tested on intervals not used to find it ([protocol](../docs/protocol_promotion_confirm.md))

Benign promotions in the non-bulk intervals 2025-07, 2025-10, 2026-01, 2026-07 and 2026-10, Ambry (771 cases)
against all other labs (156, lab-matched controls). One-sided test that Ambry's REVEL coefficient is larger.

| | Others | Ambry | Difference (SE) | One-sided p |
|---|---|---|---|---|
| REVEL, primary (lab-matched controls) | +0.58 | +0.71 | +0.14 (0.61) | 0.41 |
| REVEL, secondary (all controls) | +0.38 | +0.70 | +0.32 (0.58) | 0.29 |
| AlphaMissense, primary | +0.16 | +0.73 | +0.57 (0.25) | |
| ESM1b, primary | +0.21 | +0.11 | -0.10 (0.16) | |

**Verdict by the rule written beforehand: not confirmed.** The REVEL-specific claim is dropped. Three things to
keep straight:

- **The test was weaker than I had said.** The outcome-blind figure (a detectable difference of 0.43) was a lower
  bound. The realized standard error of 0.61 is about four times the one it implied: about 1.2 times from
  evaluating the information at the fit rather than at zero, and about 3 times from gene clustering. The
  clustering is almost all one gene: BRCA2 supplies 86% of the clustered variance of the REVEL difference (for
  ESM1b and Grantham clustering changes almost nothing). Only differences above about 1.7 would have been
  detected. The 95% interval for the REVEL difference (-1.06 to +1.34) excludes the originally discovered gap
  of about 1.6 but not smaller ones. "Not confirmed" is not "refuted". Whether the discovery estimate also
  leans on a few genes was not checked.
- **The pattern moved.** In the discovery windows Ambry's benign promotions loaded on REVEL (+1.67) and not on
  AlphaMissense (+0.12). In 2025-26 they load on both (+0.71, +0.73). The two scores are correlated within
  Ambry, so their split is unstable, while their sum is not: about 1.8 then and 1.4 now, against about 0.7 for
  the other labs. This is a post-hoc reading and the history ends in 2026-10, so there is no fresh data to test it.
- **Stop rule applied.** The VUS-promotion line of work (RQ1, option B) ends here.

What survives as a hypothesis, not a finding: some labs' benign promotions follow machine-learning scores more
closely than others', and Ambry looks like one. Separating a lab that steers by a predictor from one whose
decisions merely follow what the predictor captures needs an outside anchor (deep mutational scans, H1b).
