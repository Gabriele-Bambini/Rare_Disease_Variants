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
