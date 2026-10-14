# Performative labels in variant effect prediction

BSML 2026/27, project 15. Do computational predictors shape the ClinVar labels they are
benchmarked on? We measure this feedback loop and other forms of circularity, using
interpretable models.

- Research plan: [docs/brainstorm.md](docs/brainstorm.md)
- Data download kit for the cluster: [scripts/hpc/README.md](scripts/hpc/README.md)
- Week 2 results (ClinVar label history, what labs write): [results/README.md](results/README.md)

## Code

```bash
pip install -e '.[dev]'
pytest
```

| Path | What |
|---|---|
| `src/rdv/clinvar.py` | streaming readers for ClinVar TSVs, label and review-status normalisation, missense universe |
| `src/rdv/history.py` | per-snapshot variant states, first confident label, transitions |
| `src/rdv/scv_text.py` | predictor mentions in submitter free text, lab-mix standardisation and decomposition |
| `src/rdv/plots.py` | figures |
| `scripts/clinvar/` | command-line drivers (`build_history.py`, `scv_mentions.py`, `label_flow.py`) |
| `scripts/hpc/` | resumable download of every dataset in the plan |
