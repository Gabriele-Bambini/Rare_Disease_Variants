# Performative Labels in Variant Effect Prediction

> **Do computational predictors shape the ClinVar labels they are later benchmarked on?**
> We measure this feedback loop, and other forms of circularity, with pre-registered,
> interpretable statistical designs on ten years of ClinVar history.

BSML 2026/27 · project 15 · Python ≥ 3.10 · fully reproducible from public data

---

## Why this matters

Variant effect predictors (REVEL, AlphaMissense, ESM1b, …) are evaluated against ClinVar
labels. But clinical labs can use those same predictors when they assign the labels
(ACMG criteria PP3/BP4). If predictors influence labels, benchmark scores are inflated by
circularity. This project tests whether that loop is visible in the public record.

## Headline results

| Question | Design | Outcome |
|---|---|---|
| **Is the loop visible in what labs write?** | 3.7M submissions (SCVs) on missense variants, text analysis | **No.** The labs supplying most labels never name their tools (Ambry: 0 REVEL/PP3/BP4 mentions in 1.58M SCVs), and pooled trends mostly reflect lab mix. A placebo check passes (AlphaMissense named in 0 of 963,912 texts before its release). |
| **RQ1: did endorsed predictors gain accuracy on new labels?** | Difference-in-differences of within-gene AUC around REVEL (Dec 2022) and AlphaMissense (Sep 2023) milestones | **No detectable gain** (DiD −0.010 and −0.0004). Detectable effect ≈ 0.026–0.031 AUC; the placebo noise floor is of similar size. |
| **Protocol 2: does a VUS's score predict its promotion to P/LP?** | Conditional logistic model, pre-frozen | **Inconclusive** by the rule written beforehand (detectable change 0.26–0.29 vs. threshold 0.25). |
| **Protocol 3: is the post-hoc Ambry/REVEL lead real?** | Confirmation on unseen intervals | **Not confirmed**, and not refuted. Line of work stopped per the pre-written stop rule. |

Negative and inconclusive results are reported as such, with the reasons the design may
have been blind. Full discussion: [`results/README.md`](results/README.md).

## Repository map

```
.
├── docs/                      Research plan and frozen analysis protocols
│   ├── brainstorm.md            full research plan (hypotheses H1a–H1c, RQ1–RQ4)
│   ├── protocol.md              RQ1 protocol, frozen before any metric was computed
│   ├── protocol_promotion.md    Protocol 2: VUS promotion
│   └── protocol_promotion_confirm.md   Protocol 3: confirmation of the Ambry lead
├── src/rdv/                   Library code (installable package)
│   ├── clinvar.py               streaming ClinVar readers, label and review-status normalisation
│   ├── history.py               per-snapshot variant states, first confident label, transitions
│   ├── scv_text.py              predictor mentions in free text, lab-mix decomposition
│   ├── scores.py, aa.py         AlphaMissense / REVEL / ESM1b joins with identifier checks
│   ├── rq1.py, promotion.py     the RQ1 and promotion analyses
│   ├── clogit.py, metrics.py    conditional logit and within-gene weighted AUC
│   └── plots.py                 figures
├── scripts/
│   ├── clinvar/                 command-line drivers for each analysis step
│   └── hpc/                     resumable download kit for every dataset (SLURM-ready)
├── results/                   Output tables, figures and the written results report
└── tests/                     pytest suite (one test module per library module)
```

## Reproduce

```bash
pip install -e '.[dev]'
pytest

export DATA_DIR=/path/to/data                       # see scripts/hpc/README.md
bash scripts/hpc/download_all.sh clinvar_current
CLINVAR_ARCHIVE_KINDS=variant_summary bash scripts/hpc/download_all.sh clinvar_archives

python scripts/clinvar/build_history.py --data-dir $DATA_DIR   # universe + 40 snapshots, ~6 min on 4 cores
python scripts/clinvar/scv_mentions.py  --data-dir $DATA_DIR   # SCV text analysis, ~1 min
python scripts/clinvar/label_flow.py    --data-dir $DATA_DIR   # label flow + summary.csv
```

The RQ1 and promotion drivers (`rq1_did.py`, `promotion_*.py`, `join_scores.py`) are
documented in the corresponding protocol under [`docs/`](docs/).

## Methodological principles

- **Protocols are frozen before outcomes are computed.** Each analysis cites the commit
  that froze it; deviations are logged, never silently applied.
- **Post-hoc findings are labelled as leads**, then tested on data not used to find them.
- **Every claim names its limits**: detectable effect size, placebo noise floor, and whether
  a positive control exists.

## Status

Analysis complete through Protocol 3 (data to ClinVar 2026-10). The results report lists
the options for what follows, including making the circularity decomposition (RQ2–RQ4) the
headline and anchoring against deep mutational scans (H1b).

## Data and licensing

Raw data (ClinVar, gnomAD, AlphaMissense, REVEL, ESM1b, AlphaFold) is **not** redistributed
here; the download kit fetches it from the original sources under their own terms.
Code license: not yet specified.
