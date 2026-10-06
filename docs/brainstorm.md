# Rare disease genetics: implementation brainstorm

Source: BSML Project Proposals 2026/27, project 15, "Genetic variation and rare disease".

## 1. The question, restated

The proposal's stated contribution is to **measure circularity**: how much of the reported
accuracy of missense variant effect predictors comes from labels that were partly built
from the same evidence the models use as input.

Agree on the reframing in week one. We predict pathogenic vs benign for **missense variants
already in ClinVar**, in a fixed set of rare disease genes, using open annotations. We do not
run an association study on patient cohorts.

Working title: *How much of variant effect prediction is the label talking back?*

### Where the circularity comes from

| Type | Mechanism | How we expose it |
|---|---|---|
| **Variant overlap** (Grimm et al. 2015, type 1) | Tools trained on ClinVar or HGMD are then tested on ClinVar | Never use a precomputed score as the reference. Run a temporal split. |
| **Gene prior** (Grimm type 2) | Most genes have mostly-P or mostly-B variants, so a model can learn "which gene" rather than "which variant" | Gene-held-out splits, a gene-prior-only baseline, and per-gene (macro) AUROC |
| **AF → label** | ACMG BA1/BS1 call common variants benign; PM2 favours absent ones. Many ClinVar benign missense are benign *because* they are common | Ablate AF, and evaluate only on rare variants |
| **Predictor → label** | ACMG PP3/BP4 use in silico scores (REVEL, CADD, and since ~2023 AlphaMissense, at the ClinGen thresholds from Pejaver et al. 2022) | Ablate predictor scores. Evaluate inside the "indeterminate" score band. Run a temporal analysis. |
| **Homology** | Paralogs (SCN1A/2A/5A, KCNQ*, COL*) share sequence and pathogenic hotspots | Group genes by family or sequence cluster when splitting |

## 2. Hypotheses we can actually test

- **H1 (AF leakage).** AF is the top SHAP feature on the full test set. On the subset with
  AF < 1e-4 in every ancestry group, where BA1/BS1 could not have been applied, AF's
  contribution collapses and AUPRC drops by a large margin.
- **H2 (gene prior).** A gene-prior-only model, which predicts the training-set pathogenic
  fraction of the gene, gets a high global AUROC under random splits. Most of the gap
  between random and gene-held-out splits is this effect.
- **H3 (predictor leakage).** Models that use REVEL/CADD/SIFT/PolyPhen lose more performance
  between ClinVar and an external lab-measured benchmark (DMS) than ESM-2 zero-shot does.
- **H4 (circularity is rising).** Variants first classified after the ClinGen PP3/BP4
  calibration (late 2022) and the AlphaMissense release (Sept 2023) agree more with
  REVEL/AlphaMissense than variants classified before. This inflates AlphaMissense's measured
  accuracy on recent labels. *This is the most original angle, and cheap: ClinVar keeps
  dated archives.*
- **H5 (ancestry).** A model given only global AF makes more false-pathogenic calls on
  variants that are rare in NFE but common in AFR, SAS or EAS. Using `grpmax` FAF95 reduces
  these errors.

The deliverable can be a **circularity decomposition**, a single waterfall chart:

```
reported AUPRC (random split, all features, all variants)
  − gene prior            (random → gene-held-out)
  − homology              (gene-held-out → family-held-out)
  − AF label leakage      (all → rare-only test set)
  − predictor leakage     (all → indeterminate-band / pre-2022 labels)
  = "honest" AUPRC, checked against DMS
```

## 3. Data pipeline

All data is open and needs no data use agreement.

### 3.1 Genes
- **PanelApp (Genomics England)**: green (high-evidence) genes, through the REST API.
  Take the union of panels, then keep 300–600 genes.
- Selection rule, fixed before looking at results: at least N₁ P/LP **and** at least N₂ B/LB
  missense variants at ≥1 star (start with N₁ = N₂ = 5). Record that this rule biases the
  set towards well-studied genes.
- **Force-include genes that have clinical DMS assays** in ProteinGym/MaveDB (BRCA1, TP53,
  PTEN, MSH2, CBS, GCK, HMBS, LDLR, KCNH2, SCN5A, … check the list). These are always held
  out of training.
- Keep mode of inheritance (AD/AR/XL) from PanelApp. It matters for AF.

### 3.2 Labels: ClinVar
- `variant_summary.txt.gz` (GRCh38), or the VCF. Stream it with **polars/duckdb**, never pandas
  in full.
- Keep: missense on the **MANE Select** transcript; P, LP, B or LB; review status ≥1 star with
  no conflicts. Sensitivity analyses: ≥2 stars only, and P/B only (drop "likely").
- Deduplicate at protein level (gene, position, ref aa, alt aa). Different nucleotide changes
  that give the same amino acid change are one example.
- **Dated archives** (`vcf_GRCh38/archive_2.0/`) give each variant's first-classification date,
  for H4 and the temporal split.
- **Optional, high value:** `submission_summary.txt.gz` has free-text descriptions per
  submission (SCV). Grep them for `PP3|BP4|REVEL|in silico|computational` and
  `BA1|BS1|PM2|gnomAD|frequency` to flag which labels *say* they used predictors or AF.
  Coverage will be partial, but even partial coverage gives a direct circularity
  stratification.

### 3.3 Features
Group the features so the ablations are clean:

| Group | Features | Source / access |
|---|---|---|
| **A. Frequency** | global AF, `grpmax` FAF95, per-ancestry AF, nhomalt, an explicit "absent from gnomAD" flag (absence is informative, not missing at random), AF / max credible AF for the gene's inheritance mode (Whiffin et al. 2017) | gnomAD v4.1 exomes+genomes. Use **remote tabix** on the public per-chromosome VCFs with `pysam`, restricted to the gene regions. Do not download the full release. The GraphQL API is a fallback (rate limited). |
| **B. Conservation** | phyloP (241 mammals and 100 vertebrates), phastCons, GERP++ | UCSC bigWig through `pyBigWig` (remote reads work) |
| **C. Protein / structure** | Grantham, BLOSUM62, aa property deltas, relative position, Pfam/InterPro domain, AlphaFold pLDDT, relative solvent accessibility | UniProt, InterPro API, AlphaFold DB (PDB per protein, then DSSP or a simple neighbour count) |
| **D. Gene / phenotype** | inheritance mode, LOEUF and missense Z, HPO term count and top-level HPO categories | gnomAD constraint table, HPO `genes_to_phenotype.txt`. These are gene-level, so they mostly help across genes, which is exactly what gene-held-out tests. |
| **E. Precomputed predictors** | REVEL, CADD, SIFT, PolyPhen-2, (BayesDel, MetaRNN) | **myvariant.info** batch API serves the dbNSFP fields and avoids the dbNSFP download and licence |
| **F. PLM** | ESM-2 650M zero-shot masked-marginal log-odds (mut vs wt); embedding at the mutated position (wt context), 1280-d, reduce with PCA to ~32 for the tree models | Hugging Face `facebook/esm2_t33_650M_UR50D`, Colab T4 fp16 |
| **Comparator** | AlphaMissense | Zenodo `AlphaMissense_aa_substitutions.tsv.gz`, joined on UniProt accession + protein variant |

**Notes on ESM-2 compute.** True masked marginals need one forward pass per *variant
position*, not every position. Batch the masked copies. For proteins over 1022 aa, use
overlapping windows centred on the position. Budget: a few thousand positions across a few
hundred proteins takes a few GPU-hours on a T4. Fallback: wt-marginals (one pass per protein)
or the precomputed genome-wide ESM1b scores from Brandes et al. (2023).

**Identifier hygiene is the main engineering risk.** Fix one coordinate system early: GRCh38
genomic, plus MANE Select ENST, plus the matching UniProt canonical isoform. Check that the
reference amino acid in ClinVar's HGVS.p matches the protein sequence at that position. Drop
and count mismatches. AlphaMissense and ESM must use the same isoform.

### 3.4 External, non-circular test sets
- **ProteinGym DMS substitutions** (human clinical genes) and **MaveDB**. Labels come from lab
  assays. Metrics: Spearman between score and assay, and AUROC against the assay's own
  functional classes where they exist (e.g. BRCA1 SGE, Findlay et al. 2018).
- Do **not** use the ProteinGym *clinical* benchmark as external. It is ClinVar.
- Caveat: DMS measures one molecular function, not disease. It is a check on direction, not a
  gold standard.

## 4. Evaluation, frozen before any fitting

- **Splits**
  1. random (shown only as the inflated reference)
  2. **gene-held-out** `GroupKFold` (primary)
  3. **family-held-out**, grouping genes by HGNC gene group or MMseqs2 clusters at ≥30% identity
  4. **temporal**: train on labels first classified before 2021, test on labels from 2023 onwards
- **Test-set strata** (the circularity lenses)
  - rare-only: AF < 1e-4 in every ancestry group, or absent
  - predictor-indeterminate: REVEL in the band where ClinGen gives no PP3/BP4 (≈0.29–0.64,
    check Pejaver et al. 2022)
  - text-flagged: SCVs that mention PP3/BP4 vs those that don't (§3.2)
  - ≥2-star only
- **Metrics:** AUPRC (primary, under imbalance), AUROC, **per-gene macro AUROC** over genes
  with both classes, Brier score, reliability diagram and ECE. **Gene-level cluster bootstrap**
  for every CI, since variants within a gene are not independent.
- **Calibration note:** probabilities are calibrated to ClinVar's P:B ratio, not to a clinical
  prior. Report this, and optionally show prior-shift recalibration.
- **Sanity baselines**
  - gene-prior-only
  - "AF only" logistic regression
  - "Grantham only"
  - labels shuffled within gene (should give per-gene AUROC ≈ 0.5)

## 5. Models and ablation grid

Models: logistic regression (standardised, L2), XGBoost (`scale_pos_weight`, tuned with
nested GroupKFold), ESM-2 zero-shot (no training), ESM-2 embeddings + logistic regression,
AlphaMissense (no training).

| Run | A freq | B cons | C prot | D gene | E predictors | F PLM |
|---|:-:|:-:|:-:|:-:|:-:|:-:|
| full | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| −predictors | ✓ | ✓ | ✓ | ✓ |  | ✓ |
| −AF | | ✓ | ✓ | ✓ | ✓ | ✓ |
| −AF −predictors ("clean") | | ✓ | ✓ | ✓ | | ✓ |
| PLM only | | | | | | ✓ |
| classic (no PLM) | ✓ | ✓ | ✓ | ✓ | ✓ | |

Cross every run with every split and every stratum in §4 into one long results table
(`run × split × stratum × metric`), then plot. Use SHAP for the XGBoost runs, especially to
show AF dominance moving between the full and rare-only strata.

**Optional controlled experiment.** A simulated circular labeller. Take the DMS genes, where
the true functional effect is known. Simulate ACMG point-based labelling (Tavtigian et al.
2020) in which a fraction *p* of variants receive PP3/BP4 points from REVEL. Train and test
on the simulated labels and plot measured AUROC against *p*. This shows the inflation
mechanism with a known ground truth.

## 6. Suggested repo layout

```
config/            genes.yaml, thresholds.yaml, splits.yaml (frozen in week 3)
src/data/          panelapp.py, clinvar.py, gnomad_tabix.py, conservation.py,
                   dbnsfp_myvariant.py, alphamissense.py, dms.py, idmap.py
src/features/      protein.py, esm2.py (zero-shot + embeddings)
src/eval/          splits.py, strata.py, metrics.py (cluster bootstrap), calibration.py
src/models/        baselines.py, logreg.py, xgb.py
notebooks/         exploration only, no logic that results depend on
data/              gitignored; raw/ interim/ processed/ (parquet)
results/           results table + figures
```
Tooling: Python 3.11, polars/duckdb, pysam, pyBigWig, requests, scikit-learn, xgboost, shap,
torch + transformers. Use a `Makefile` or Snakemake so `make data features eval` rebuilds
everything.

## 7. Timeline (~13 weeks, team of 2–5)

| Week | Milestone |
|---|---|
| 1 | Agree the reframing. Read ACMG/AMP 2015, Grimm 2015, Pejaver 2022, the AlphaMissense paper and Meier 2021 (ESM-1v). The genetics glossary owner writes a one-page vocab sheet. |
| 2 | Gene list (PanelApp + DMS genes). ClinVar filtered to MANE missense. Identifier mapping with mismatch report. |
| 3 | **Freeze the evaluation**: splits, strata, metrics, baselines, committed before any model is trained. |
| 4–5 | Feature groups A–E (gnomAD tabix, bigWig, myvariant, structure). Gene-prior and single-feature baselines. |
| 6 | Logistic regression + XGBoost full run, SHAP. First random-vs-gene-held-out gap (H2). |
| 7–8 | Ablation grid and strata (H1, H3). ClinVar archive dating for H4. |
| 8–9 | ESM-2 zero-shot + embeddings on Colab. AlphaMissense on identical splits. |
| 10 | DMS external evaluation. Ancestry analysis (H5). |
| 11 | Optional: simulated circular labeller, SCV text flags. |
| 12–13 | Circularity waterfall, write-up, cleanup. |

**Cut order if time runs short:** simulated labeller, then SCV text mining, then
family-held-out split, then ESM embeddings (keep zero-shot).

**Roles (for 2 people):** (1) data and genetics: ClinVar, gnomAD, ID mapping, DMS;
(2) modelling and evaluation: splits, metrics, models, ESM. Swap code review weekly.

## 8. Risks and mitigations

| Risk | Mitigation |
|---|---|
| No truly independent test set | DMS as external check. Strata and temporal splits as internal lenses. State the limitation explicitly. |
| Isoform and coordinate mismatches silently corrupt the joins | MANE-only, ref-aa check, report the dropped counts |
| gnomAD scale | Remote tabix on gene regions only, cached as parquet |
| Too few benign variants in rare genes | ≥1-star B/LB, per-gene macro AUROC only over eligible genes. Report the class balance per gene. |
| Ancestry confounding | Use per-group AF and grpmax FAF95. Error analysis by ancestry-differential AF (H5). |
| DMS ≠ disease | Use it for direction and ranking only. Report per-assay results, never a pooled "accuracy". |
| Leakage through tuning | Nested GroupKFold. The test strata are never used for model selection. |

## 9. Key references
- Richards et al. 2015, ACMG/AMP variant interpretation guidelines
- Grimm et al. 2015, *Hum Mutat*: two types of circularity in predictor evaluation
- Pejaver et al. 2022, *AJHG*: calibration of computational tools for PP3/BP4
- Tavtigian et al. 2020: Bayesian points system for ACMG
- Whiffin et al. 2017: maximum credible allele frequency
- Cheng et al. 2023, *Science*: AlphaMissense
- Meier et al. 2021 (ESM-1v zero-shot); Lin et al. 2023 (ESM-2); Brandes et al. 2023 (ESM1b genome-wide)
- Notin et al. 2023: ProteinGym; Esposito et al. 2019: MaveDB
- Chen et al. 2024, gnomAD v4
