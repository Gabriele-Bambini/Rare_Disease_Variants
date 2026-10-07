# Performative labels in variant effect prediction

*Research plan. BSML 2026/27, project 15 ("Genetic variation and rare disease").
Supersedes the first brainstorm (see git history for the earlier version).*

## 0. Pitch

Clinical labs now use computational predictors as evidence when they classify genetic
variants. The ACMG criterion PP3/BP4 was calibrated for REVEL and others in late 2022, and
for AlphaMissense in 2025. Those classifications land in ClinVar. ClinVar is then the main
benchmark the same predictors are judged on. **The predictor helps write the answer key it
is graded against.**

This is *performative prediction* (Perdomo et al., 2020) observed in the wild, with
timestamps.

We measure three things:

1. whether this feedback loop is already visible in ClinVar;
2. how much of the apparent accuracy of variant effect predictors it, and older forms of
   circularity, account for;
3. which biological signals survive once the circular ones are removed.

We do this with an **interpretable model**, so the circularity can be *seen* in the learned
feature effects, not only inferred from a drop in AUROC.

**Target outputs:** a workshop paper (ML for biology/health), a preprint, and an
open-source package with leakage-aware benchmark splits.

## 1. Positioning: what is known, what is new

| Already established | Reference |
|---|---|
| Two classic circularity types: variant overlap and gene-level label imbalance | Grimm et al. 2015 |
| Gene identity dominates ClinVar benchmarks; predictor rankings change under within-gene evaluation | bioRxiv preprint, Aug 2026 ("Gene identity, not variant effect, dominates ClinVar benchmarks…"). **Read in week 1.** |
| DMS (deep mutational scanning) assays are a non-circular check on predictors | ProteinGym (Notin et al. 2023); Genome Biology 2025 |
| Calibrated predictor thresholds formally enter clinical classification | Pejaver et al. 2022; ClinGen 2025 extension incl. AlphaMissense |
| Interpretable predictors exist, but they aim at accuracy, not at measuring circularity | e.g. DAVE (medRxiv 2025), SHAP-on-GBM tools |
| Performative prediction: deployed models shift the data distribution they are evaluated on | Perdomo et al. 2020 and follow-ups. Almost all evidence is theoretical or simulated. |

**Our contribution (gaps found in a quick literature search; confirm with a proper
search in week 1):**

- **C1, the main one.** Empirical evidence, with a causal-style design, of a *temporal
  feedback loop* between predictors and ClinVar labels.
- **C2.** A decomposition of reported accuracy into its circular components: gene
  (replicating the 2026 preprint), homology, allele frequency, predictor scores, and
  feedback.
- **C3.** An interpretable model (EBM and a learned integer scorecard) that makes the
  leakage visible, and lets us compare data-learned evidence weights with the ACMG points
  system.
- **C4.** Released splits and code, so other groups can evaluate without these leaks.

## 2. Research questions and hypotheses

**RQ1, feedback loop (the headline).** Have labels become more predictable *specifically* by
the predictors that clinical guidelines endorse, after those endorsements?
- *H1a.* After the calibration milestones (Dec 2022 for REVEL and others; Sept 2023 for the
  AlphaMissense release; 2025 for its ClinGen calibration), agreement between new ClinVar
  labels and endorsed predictors rises **more than** agreement with non-endorsed control
  signals.
- *H1b.* Variants resolved from VUS to P/LP or B/LB after these dates are resolved in the
  direction of the endorsed predictor more often than DMS evidence would justify.
- *H1c.* Mentions of PP3/BP4, REVEL and AlphaMissense in ClinVar submission texts increase
  over time. This is descriptive, and gives direct evidence of usage.

**RQ2, allele-frequency leakage.** Rarity drives benign labels (BA1/BS1) and
pathogenic support (PM2).
- *H2.* The learned effect of AF is strongest on the full test set and collapses on the
  rare-only stratum (AF < 1e-4 in every ancestry group).

**RQ3, gene identity (replication, not novelty).**
- *H3.* Random split ≫ gene-held-out ≳ family-held-out. A gene-prior-only baseline reaches
  a high global AUROC.

**RQ4, what survives.**
- *H4.* After all leaks are removed, conservation, structure and the ESM-2 zero-shot score
  keep most of the signal measured against DMS. A small interpretable model comes close to
  AlphaMissense on the honest benchmark.

## 3. RQ1 in detail: measuring the feedback loop

This is the most original part, and also the one most exposed to confounding. The design
has to be convincing.

### 3.1 Building the label history
- From the monthly ClinVar archives (`tab_delimited/archive/variant_summary_YYYY-MM.txt.gz`
  or `vcf_GRCh38/archive_2.0/`), reconstruct for each missense variant:
  - the date it first received a confident P/LP or B/LB classification;
  - any later reclassification (VUS → P/LP, VUS → B/LB, P ↔ B).
- Group variants into **quarterly cohorts** by date of first confident classification,
  2017 → 2026.
- Get exact submission dates from `submission_summary.txt.gz`, which has
  `DateLastEvaluated` per SCV.

### 3.2 The central comparison: difference-in-differences

Raw AUROC of AlphaMissense per cohort is **not** enough. The mix of variants classified
changes over time (new genes, easier or harder cases), and AlphaMissense might simply be
right. So we compare **endorsed** predictors with **control** signals:

| Endorsed (used for PP3/BP4) | Controls (not used for PP3/BP4, or not used directly) |
|---|---|
| REVEL, AlphaMissense, BayesDel, CADD | ESM-2 / ESM1b zero-shot (precomputed, Brandes et al. 2023), a conservation-only score (phyloP), Grantham, **DMS scores** where available |

Estimand, roughly:

```
Δ = [AUROC(endorsed, after) − AUROC(endorsed, before)]
  − [AUROC(control,  after) − AUROC(control,  before)]
```

All AUROCs are computed **within gene** (macro over genes), so a changing gene mix cannot
drive Δ. Results are shown as an **event-study plot**: Δ per quarter, relative to each
milestone. Confidence intervals come from a cluster bootstrap over genes.

How to read the outcomes:
- If labels just became "easier", every signal improves and Δ ≈ 0.
- If labels were partly *written by* the endorsed predictors, those gain more and Δ > 0.

### 3.3 Additional evidence
- **DMS anchor (H1b).** For genes with DMS data, take the variants resolved after the
  milestones. Do the labels follow AlphaMissense in the cases where AlphaMissense and the
  DMS disagree? This is the sharpest single test.
- **Text evidence (H1c).** Plot the share of SCVs per quarter whose free-text description
  mentions `PP3|BP4|REVEL|AlphaMissense|in silico`. It is simple and very readable. Use it
  as figure 2 of the paper.
- **Robustness.**
  - ≥2-star labels only, expert panels only (ClinGen VCEPs, which explicitly apply calibrated
    PP3);
  - excluding genes that enter ClinVar for the first time after the milestone;
  - alternative milestone dates (placebo dates should show no effect).

### 3.4 Simulation: the loop with a known ground truth
For ML reviewers, and as a sanity check:
- Use DMS genes as "truth".
- Simulate labelling rounds: a labeller combines noisy true evidence with PP3 points from a
  predictor, with weight *w*, using the Tavtigian 2020 Bayesian points.
- Each round, the predictor is re-evaluated, and optionally re-thresholded, on the
  accumulated labels.
- Show how the benchmark score inflates, and how predictor *rankings* distort, as a function
  of *w* and the number of rounds.
- This connects our empirical Δ to the performative-prediction literature, and shows the
  failure mode cleanly.

### 3.5 What if RQ1 comes back null?
That is still a result: "no detectable feedback yet, with this power", which is useful for
guideline bodies. The rest of the project (RQ2–RQ4) stands on its own. Report the minimum
detectable effect.

## 4. The interpretable model (RQ2–RQ4)

### 4.1 Feature groups (one prefix per group, so ablations are column selections)

| Group | Feature | Biology in one line | Circularity risk |
|---|---|---|---|
| A | AF: global, `grpmax` FAF95, per-ancestry, nhomalt, absent flag | Variants common in healthy people rarely cause rare disease | **High** (BA1/BS1/PM2) |
| B | Conservation: phyloP (241 mammals), GERP | Positions unchanged across evolution matter | Low |
| C | Chemistry: Grantham, BLOSUM62, domain | Drastic substitutions in functional regions break proteins | Low |
| D | 3D structure graph from AlphaFold: contact degree, centrality, distance to annotated functional sites, pLDDT | Buried and central residues are fragile | Low |
| E | ESM-2 zero-shot log-odds (one number) | Protein "grammar" learned from evolution | Low (never saw clinical labels) |
| G | Gene: inheritance mode, LOEUF, missense Z | Recessive genes tolerate higher carrier frequency | Medium (gene shortcut) |
| P | Endorsed predictors: REVEL, CADD, AlphaMissense | Used directly in PP3/BP4 | **High** (comparator and ablation only) |
| H | *Experiment only:* count of known pathogenic variants within 8 Å | = the PM1 "hotspot" criterion | **Label leakage by construction** |

### 4.2 Models
- **EBM (Explainable Boosting Machine, `interpret`)** is the main model. It gives one shape
  function per feature plus a few pairwise interactions, notably AF × inheritance mode.
  Every prediction decomposes into readable contributions.
- **Learned integer scorecard.** Sparse logistic regression with rounded coefficients
  (or RiskSLIM if feasible), with 5–8 features. We compare the learned points with the ACMG
  points system, before and after removing circular labels.
- **Reference models**, used to measure the cost of interpretability: XGBoost on the same
  features, AlphaMissense, REVEL.
- **Calibration**: isotonic or Platt inside nested CV, stated as calibrated to ClinVar's
  P:B ratio and not to a clinical prior.

### 4.3 How interpretability shows the leakage
- AF shape function on the full set vs the rare-only stratum: it should flatten (H2).
- The scorecard's points for AF and for predictors, before and after removing circular
  labels, compared with the ACMG points.
- Adding feature group H inflates performance under random splits and vanishes under
  gene-held-out splits. This makes the circularity directly visible in the model.

## 5. Evaluation protocol (frozen before any fitting, in week 3)

- **Splits:**
  - random (only as the inflated reference);
  - gene-held-out (primary);
  - family-held-out, using connected components of a gene-similarity graph (MMseqs2 at ≥30%
    identity, or HGNC gene groups);
  - temporal.
- **Strata:**
  - rare-only;
  - predictor-indeterminate band (REVEL ≈ 0.29–0.64, check against Pejaver 2022);
  - SCVs whose text does or does not mention PP3/BP4;
  - ≥2 stars;
  - DMS genes.
- **Metrics:**
  - AUPRC (primary), AUROC, per-gene macro AUROC;
  - Brier score, reliability diagram;
  - Spearman correlation with DMS;
  - cluster bootstrap over genes for every CI.
- **Sanity baselines:**
  - gene prior only;
  - AF only;
  - Grantham only;
  - labels shuffled within gene.

The final artefact is a **circularity waterfall**: reported AUPRC, minus gene, minus
homology, minus AF, minus predictors, minus feedback, giving the honest AUPRC, checked
against DMS.

## 6. Data

All open. **Download everything on the cluster with `scripts/hpc/download_all.sh`**
(see `scripts/hpc/README.md`).

| Data | Use | Access |
|---|---|---|
| ClinVar `variant_summary` monthly archives and `submission_summary` | Labels, label history, SCV texts | NCBI FTP. Stream with polars/duckdb. |
| PanelApp green genes | Gene set, inheritance mode | REST API |
| gnomAD v4.1 exomes+genomes | Group A | Remote tabix on gene regions (`pysam`) |
| phyloP / GERP bigWig | Group B | UCSC, `pyBigWig` |
| AlphaFold DB structures, UniProt features | Groups C, D | REST; contact graph with `biopython` + `networkx` |
| ESM-2 650M / precomputed ESM1b (Brandes 2023) | Group E and the RQ1 control | Hugging Face / published downloads |
| gnomAD constraint | Group G | gnomAD downloads |
| REVEL, CADD, BayesDel (dbNSFP fields) | Group P, RQ1 endorsed | myvariant.info batch API |
| AlphaMissense | Group P, RQ1 endorsed | Zenodo |
| ProteinGym DMS, MaveDB | External truth, RQ1 anchor, simulation | Downloads / API |

**Identifier hygiene.** Use MANE Select transcript ↔ UniProt canonical isoform throughout.
Check that the reference amino acid matches, and log the dropped variants.

## 7. Timeline (~13 weeks)

The fastest path to the headline result comes first. RQ1 needs **no gnomAD and no GPU**:
only ClinVar archives, AlphaMissense, REVEL and precomputed ESM1b.

| Week | Milestone |
|---|---|
| 1 | Read the Aug 2026 gene-identity preprint, Perdomo 2020, Pejaver 2022, ACMG 2015 and Grimm 2015. Do a proper literature search on "ClinVar circularity temporal / PP3 feedback". Write the genetics vocabulary sheet. |
| 2 | Gene set. ClinVar archives parsed into a label-history table. **First figure: SCV text mentions over time (H1c).** |
| 3 | Freeze the evaluation protocol. Join AlphaMissense, REVEL and ESM1b. |
| 4–5 | **RQ1 MVP: event-study / DiD plot (H1a).** DMS anchor (H1b). Decide go/no-go on RQ1 as the headline. |
| 6–7 | Feature groups A–D and G (gnomAD tabix, bigWig, AlphaFold graph). Baselines. Gene-held-out vs random (H3). |
| 8–9 | EBM and scorecard. AF strata (H2). Group-H leakage experiment. XGBoost reference. |
| 10 | DMS external evaluation (H4). Circularity waterfall. |
| 11 | Simulation of the loop (§3.4). Robustness checks for RQ1. |
| 12–13 | Workshop paper draft (4–6 pages), preprint, package release, cleanup. |

**Cut order if short on time:** simulation, family-held-out split, group D, scorecard
(keep EBM).

**Roles for 2 people:**
1. Label history, RQ1 and DMS.
2. Features, interpretable models and evaluation.

Code review is weekly and swapped between the two.

## 8. Deliverables

- **Paper figures (plan them now):**
  1. Schematic of the loop: predictor → PP3 → ClinVar → benchmark → predictor.
  2. SCV mentions of PP3/REVEL/AlphaMissense over time.
  3. Event-study plot of Δ around the milestones.
  4. Circularity waterfall.
  5. EBM shape function of AF: full set vs rare-only.
  6. Learned scorecard vs ACMG points.
- **Venues:** ML for computational biology or health workshops (NeurIPS / ICLR workshops,
  MLCB, ML4H). Check the deadlines in week 1 and pick one as the hard internal deadline.
- **Code:** `pip install`-able package with the splits, strata and evaluation functions,
  compatible with the ProteinGym format.

## 9. Risks

| Risk | Mitigation |
|---|---|
| RQ1 confounded by changing case mix | Within-gene AUROC, DiD against control signals, placebo dates, excluding newly added genes |
| The control signals are themselves partly used by labs (e.g. conservation under PP3) | Use DMS as the cleanest control. Report each control separately. |
| Too few post-2023 confident labels in DMS genes | Pool milestones, use 1-star labels with a sensitivity analysis, report the minimum detectable effect |
| Overlap with the Aug 2026 preprint | Treat gene identity as a replication. The headline is RQ1, which it does not cover (verify). |
| Coordinate and isoform mismatches | MANE only, ref-aa check, dropped-variant log |
| Ancestry bias in AF | Per-group AF and `grpmax`. Error analysis on ancestry-differential variants. |

## 10. Key references
- Perdomo, Zrnic, Mendler-Dünner, Hardt 2020, *Performative Prediction* (ICML)
- Grimm et al. 2015, *Hum Mutat*: circularity in predictor evaluation
- *Gene identity, not variant effect, dominates ClinVar benchmarks of missense pathogenicity
  predictors*, bioRxiv 2026
- Richards et al. 2015: ACMG/AMP guidelines. Tavtigian et al. 2020: Bayesian points.
- Pejaver et al. 2022, *AJHG*: PP3/BP4 calibration. ClinGen 2025: calibration of additional
  tools.
- Cheng et al. 2023, *Science*: AlphaMissense. Brandes et al. 2023: ESM1b genome-wide.
  Lin et al. 2023: ESM-2.
- Notin et al. 2023: ProteinGym. Esposito et al. 2019: MaveDB.
- Lou et al. 2013 / Nori et al. 2019: GA²M / InterpretML (EBM). Ustun & Rudin 2019: RiskSLIM.
- Whiffin et al. 2017: maximum credible allele frequency. Chen et al. 2024: gnomAD v4.
