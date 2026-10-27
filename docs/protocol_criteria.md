# Protocol 4 for RQ1: how much of ClinVar's labelling rests on computational evidence, and what that does to a benchmark

**Status when frozen (2026-10-27).** Counts, the share of submissions whose listed criteria reproduce their own
class (E0 below) and a parser debugging pass exist. **No estimate of E1 to E5 has been computed with the final
parser.** The commit that adds this file precedes the commit that does.

## 0. Why this protocol exists, and what it can and cannot claim

Protocols 1 to 3 looked for a *change* in how ClinVar labels line up with predictor scores after predictor
milestones. They returned a bounded null, an inconclusive result and a lead that did not replicate
([results/README.md](../results/README.md)). They could not see explicit use. This protocol measures explicit
use directly. Many submitters list the ACMG/AMP evidence codes behind a classification, and two of those codes
*are* a computational prediction: **PP3** (supports pathogenic) and **BP4** (supports benign). That the codes
appear is certain, so the question is not whether a loop exists. It is **how much of a label hangs on them**
(a prevalence) and **what a benchmark looks like without those labels** (an effect size).

What a positive result would say: *among ClinVar missense labels whose evidence is documented, X% would lose
their class if the computational code were removed, and benchmarking REVEL (AlphaMissense) on ClinVar without
them lowers its AUC by Y.* What it would **not** say: that those labels are wrong, that a predictor's true
accuracy is Y lower, or anything about the 89% of pathogenic and 93% of benign submissions that list no codes
(sections 5 and 9). Removing agreement-by-construction labels also removes cases where predictor and truth agree, so Y
is a measure of dependence, not of bias; only an external anchor (deep mutational scans) separates the two.

### Disclosures, in the order they happened

1. **An exploratory pass with an earlier, buggy parser** (65,279 submissions with two or more codes; it missed
   codes glued to a strength by an underscore and let a negation spill across list items) gave a PP3-decisive
   share of about 15%. The bugs were found by unit tests, not by that figure. The thresholds below were set
   after seeing it; the 5% materiality line is deliberately far from it and constrains little. What this
   protocol fixes is the definitions.
2. **Parser debugging on rows I read** (no aggregate computed). A survey of the words that follow a code in
   115,902 submissions found `_STR`, `_SUPP`, `_VSTR`, `moderated`, `Supporti`, "downgraded to" and "upgraded
   to" notations that the first version missed, and a rule that let a later bare mention overwrite a stated
   strength. Two random samples of parsed rows were read for errors (seeds 20261027, 50 rows, and 8675309, 60
   rows). **The first sample printed, per row, which codes were decisive: in 6 of its 50 rows a computational
   code was decisive (PP3 in 3, BP4 in 3).** That is a row-level glimpse of what E1 and E2 will count, not an
   aggregate; it was used only to read the parse. The second sample printed derived classes and no decisive
   codes, and showed no further parse error. Both samples showed many benign calls that list a single code
   (`GENE: BP4`, `GENE: BS1`), which reproduces the class under the points system only. I therefore knew
   before freezing that benign calls resting on BP4 alone exist, and the 5% line below is not blind to that.
3. **E0 was computed before freezing**, with the exclusions of section 1 applied (it concerns whether a list
   reproduces its own class, not PP3 or BP4) and it fixes the primary reading in section 3.

## 1. Data

- ClinVar `submission_summary` and `variant_summary` as downloaded on 2026-10-12 (same release, MD5 in the
  manifest). Missense universe: GRCh38 single-nucleotide missense variants (2,552,260; `missense_universe`).
  Aggregate label and review stars from `variant_summary`.
- Scores joined earlier (`scores.parquet`): REVEL v1.3, AlphaMissense (hg38), ESM1b LLR.
- Text per submission = `Description` + `ExplanationOfInterpretation`. Submissions with review status
  "flagged submission" or "no classification provided" are excluded (they do not enter the aggregate).

## 2. Unit and parsing

- **Submission (SCV)** on a universe variant. Stated class = `classify_significance` of its classification:
  PLP, BLB, VUS or other. Lab = `Submitter`. Year = year of `DateLastEvaluated`.
- **Listing submission:** its text contains at least one non-negated ACMG/AMP code, read by
  `rdv.criteria.parse_criteria` **as committed with this file** (61 unit tests, hand-computed). Features: codes
  `PVS1, PS1-4, PM1-6, PP1-5, BA1, BS1-4, BP1-7` with a stated strength (`_Strong`, `:Moderate`, `(supporting)`,
  `_STR`, `_SUP`, `_P`, "downgraded to...") or the code's default strength; negation scoped to a list item;
  a stated strength outranks a bare mention. A change to this file after the freeze is a deviation and requires
  rerunning everything.
- `scripts/clinvar/criteria_prepare.py` found 3,731,912 submissions on the universe; after the exclusion above,
  3,712,734, of which 115,803 are listing, on 100,374 variants from 781 labs. Listing submissions: 26,093 PLP,
  36,799 BLB, 52,812 VUS, 99 other.

## 3. Reading a list as a class, and the primary readings

Four readings of a list: Richards 2015 combining rules (`r`) or the points of Tavtigian 2018/2020 (`t`: 8/4/2/1,
P at 10 or more, LP 6-9, VUS 0-5, LB -1 to -6, B -7 or less), each with an unqualified PM2 read as moderate or
(`_pm2s`) as supporting, as the ClinGen SVI recommendation of 2020 allows. A stated strength is never changed.

A listing PLP or BLB submission is **consistent** under a reading if the class family (PLP / VUS / BLB) the
list derives equals the stated family. **E0 (computed before freezing):**

| Stated | Listing | `r` | `r_pm2s` | `t` | `t_pm2s` | any of the four |
|---|---|---|---|---|---|---|
| PLP | 26,093 | **94.3%** | 84.7% | 93.9% | 84.6% | 94.5% |
| BLB | 36,799 | 58.0% | 58.0% | 99.1% | **99.4%** | 99.7% |

The reading that reproduces the most stated classes is the best available model of the rule submitters used. It
is the **primary reading**: **`r` for PLP, `t_pm2s` for BLB**. (The benign gap is the points system: a single
supporting or strong benign code reaches likely benign there but not under the 2015 rules, and "GENE: BP4" or
"GENE: BS1" are common. For BLB, `t` and `t_pm2s` differ by 0.3 points because benign lists rarely contain PM2.)
All four readings are reported for everything below.

## 4. Decisive

For a consistent submission and a code c on its list, c is **decisive** if deleting c changes the derived class
family (`rdv.criteria.is_decisive`). The **computational** code is PP3 for a PLP submission and BP4 for a BLB
submission. A PP3 listed by a BLB submission, or a BP4 by a PLP one, still enters the derivation but is not
scored as the computational code. "Other" decisive: some non-computational code is decisive (BA1 alone is
decisive for its own class).

## 5. Estimands

Everything is a ratio over submissions with a **lab-cluster bootstrap** (labs resampled with replacement, 4,000
replicates, seed 20261027, percentile 95% interval) unless stated.

- **E0 coverage and consistency.** Share of all PLP (BLB) submissions on the universe that list a code (from
  `scv_denominator`); the table in section 3.
- **E1 (headline, pathogenic).** Among consistent PLP listing submissions, primary reading: the share in which
  **PP3 is decisive**. **E1c:** among those that list PP3, the share in which it is decisive.
- **E2 (headline, benign).** Same for BLB and BP4. **E2c** likewise.
- **E1v, E2v (variant level).** Among benchmark variants (section 6) labelled PLP (BLB), the share that have at
  least one consistent submission of the same family in which PP3 (BP4) is decisive; and the same share among
  variants with at least one consistent listing submission.
- **Bracket for all submissions.** With N all submissions of the family, K the consistent listing ones and D the
  computationally decisive ones: the *documented* share D/N is a lower bound and D/K is the value if unseen
  submissions behave like seen ones. Nothing in between is claimed.
- **E3 (does the parse track the scores? a validity gate).** Among listing submissions of any class that list
  exactly one of PP3 and BP4, with a score: the within-gene weighted AUC (`rdv.metrics`, unit weights, gene
  bootstrap with 2,000 replicates, seed 20261028) of REVEL for "PP3 vs BP4"; also AlphaMissense and ESM1b; and
  the same restricted to VUS submissions (so the class label cannot carry it). **Gate:** if REVEL's AUC in the
  VUS-only version is below 0.65, the audit stops: either the parse or the premise fails and E1, E2, E5 are not
  interpreted.
- **E4 (secondary; one pre-specified direction).** Alignment with the calibrated REVEL threshold: among listing
  submissions that list PP3 (BP4) with a REVEL score, the share with REVEL at least 0.644 (at most 0.290), the
  supporting-strength bounds of the 2022 ClinGen calibration (Pejaver et al.). Contrast **2023-01-01 to
  2026-10-12** minus **2019-01-01 to 2022-02-28** by evaluation date; the preprint (2022-03-17) and journal
  (2022-10-20, 2022-12-01) months are dropped. Pooled over submissions, and a within-lab version (labs with at
  least 30 qualifying submissions in both windows, labs weighted equally); lab bootstrap; one-sided p for a
  positive contrast. Pre-trend check: 2019-01 to 2020-06 vs 2020-07 to 2022-02. Reading as an effect of the
  calibration requires p below 0.05 and a pre-trend contrast below half the main one.
- **E5 (consequence for a benchmark), section 6.**

## 6. E5: what a benchmark looks like without the dependent labels

**Benchmark V:** universe variants with aggregate label PLP or BLB, one star or more, REVEL, AlphaMissense and
ESM1b all present. Statistic: the within-gene weighted AUC of `rdv.metrics` (genes with both classes; unit
weights; pairs weighted as in protocols 1 to 3).

- **U** = variants of V with at least one consistent listing submission of the variant's family (primary reading).
- **A** = variants of U with at least one such submission in which the computational code is decisive.
- **A_pre** = variants of A whose computationally decisive submissions are all dated before **2023-09-19**
  (AlphaMissense release), no date missing. **U_pre** = variants of U whose consistent submissions are all dated
  before it. AlphaMissense cannot have been used for any label in A_pre: this is label **inheritance** from
  earlier tools.
- **O** = variants of U not in A with at least one submission in which a non-computational code is decisive
  (placebo).

For a tool t and a removed set X: Δ_t(X) = AUC_t(V \ X) - AUC_t(V). **Null:** remove at random the same number
of variants in each (gene, class) cell as X has, drawn from U (from U_pre for A_pre), 1,000 draws, seed
20261029; the same draws for all tools. **D_t(X) = Δ_t(X) - mean of the null Δ.** The one-sided permutation p is
(1 + number of null Δ at or below the observed) / 1,001. The 95% interval for D comes from a gene bootstrap
(2,000 replicates) applied to the observed and the null together.

- **P1:** t = REVEL, X = A. **P2:** t = AlphaMissense, X = A_pre. Holm over the two, alpha 0.05, and an effect
  floor |D| of at least 0.005 AUC (a judgement: small beside the gaps usually reported between leading
  predictors, large enough not to be a rounding artefact; it is not derived from data).
- Secondary: ESM1b with A and with A_pre (ESM1b is rarely cited, so the expectation is a smaller drop than
  REVEL's); REVEL and AlphaMissense with A within U only (where contamination is observed); X = O (the drop
  **per 1,000 variants removed** for A, O and the null, to compare removals of different sizes); stars of two or
  more only; expert panels excluded.

## 7. Sensitivity analyses (reported together, no multiplicity claim)

All four readings for E1, E2, E1v, E2v and P1; expert panels (`reviewed by expert panel`, `practice guideline`)
excluded; clinical-testing submissions only; evaluation years 2019 to 2026 only; the five labs with the most
consistent submissions of the family excluded; a lab-level view (median and interquartile range of lab-level
shares over labs with at least 30 consistent submissions, and the share of those labs with a lab-level share of
5% or more); E1/E2 by evaluation year with a decomposition into within-lab change and change of which labs
submit; the codes that co-occur when PP3 / BP4 is decisive. **E2d, pre-specified because it will be the first question
asked of E2:** among BP4-decisive benign submissions, the share in which BP4 is the only code listed (under the
points reading such a list is consistent and decisive by construction); and its counterpart **E1d**, the
distribution of the number of listed codes among PP3-decisive pathogenic submissions.

## 8. Decision rule (fixed now)

- **Positive on prevalence** (separately for E1 and E2): the lower 95% bound of the estimate under the primary
  reading is **5% or more**, the point estimate under each of the other three readings is 5% or more, and the
  E3 gate passes. It is **broad** if it still holds with the five largest labs excluded, and **concentrated**
  if it does not. Below that: **not positive**, and the number is reported as it is.
- **Positive on consequence:** P1 or P2 passes (Holm p below 0.05 and D at or below -0.005). It is **specific to
  computational evidence** only if the drop per 1,000 removed variants for A is larger than for O, with a gene-
  bootstrap interval for the difference that excludes zero; otherwise the text says removal of evidence-rich
  variants lowers AUC whatever the code.
- **Language.** A prevalence claim is made for submissions that list their criteria only. The bracket in
  section 5 is the only statement about all submissions.

## 9. Threats stated in advance

- **Listing submissions are 10.9% of PLP and 7.4% of BLB submissions** (E0) and come from labs that are not a
  random sample (larger and commercial labs list more); nothing here covers the rest beyond the bracket.
- **"Decisive" is arithmetic given the list**, under a published rule. A submitter may use a gene-specific
  specification, list a subset of its evidence, or reach its class by judgement. Requiring the list to
  reproduce the class removes the worst cases (E0) but not these. It measures dependence on the code, not the
  submitter's reasoning.
- **Template labs.** One laboratory applying the same list to thousands of variants moves a pooled share; the
  lab bootstrap and the top-five exclusion are the guards, and the lab-level view shows breadth.
- **BP4-alone benign calls** are consistent only under the points system and are decisive by construction there;
  the 2015 reading drops them. Both are reported; the primary reading for BLB is the one that reproduces nearly
  all benign lists (99.4%).
- **E5 measures dependence, not bias** (section 0). Selection works in both directions: removed labels are the
  ones where tools agree, and what remains is the evidence-rich, easier-for-experts set. The matched random
  removal and the placebo O bound the selection from outside; neither removes it.
- **Dates.** `DateLastEvaluated` is the last evaluation; an old date means the text was not rewritten since.
- **Parse errors** remain in unusual notations and in prose that never names a code. E3 tests whether the
  parse tracks the scores; it does not give a recall figure.
- **Multiple readings.** The primary reading was chosen on E0 (validity), not on any decisive share.
- **Release.** One snapshot (2026-10-12). A different release would change the counts.

## 10. Prior work, and what is new (searched 2026-10-27 with PubMed keyword queries and web search; both weak)

Circularity in predictor benchmarking was named by Grimm et al. ([DOI](https://doi.org/10.1002/humu.22768),
Hum Mutat 2015): the same variants, or different variants of the same protein, appear in training and evaluation.
That is circularity between the *predictor's* data and the benchmark. Here the loop runs through the
**labels**: predictors feed the evidence that makes the labels. I found (as summarised by a text fetch of the
full texts, not checked against the originals) **the nearest precedent**, a single-laboratory study of
recalibrated PP3/BP4 thresholds on missense VUS (Dubé et al., Genet Med Open 2026). It reports PP3 in 19-44% and
BP4 in 7-33% of the laboratory's missense variants and that removing the laboratory's own PP3/BP4 evidence
reclassified about 6% of those VUS; the same summary flagged an inconsistency in the paper (it says the analysis
was limited to variants whose initial classification included "BP3 or PP4"). E1 and E2 are the ClinVar-wide,
multi-laboratory analogue of that counterfactual, for P/LP and B/LB labels instead of VUS, and 6% is an outside
magnitude to compare with. I also found a
benchmark of predictors in five cancer genes that says dependence on in silico evidence "could not be determined
from public data" for BRCA1, BRCA2 and ATM and builds an alternative truth set excluding PP3-reliant calls for
TP53 and ClinGen (Nelson et al., Hum Mutat 2026). I found no ClinVar-wide count of labels that depend on PP3/BP4
or measure of its effect on a benchmark. Absence in two weak searches is not proof of absence.

## 11. Outputs

`results/criteria_*.csv` and figures from `scripts/clinvar/criteria_audit.py` (E0 to E4) and
`scripts/clinvar/criteria_benchmark.py` (E5), written after this commit; a section in `results/README.md`.

## Deviations log

None at freezing.
