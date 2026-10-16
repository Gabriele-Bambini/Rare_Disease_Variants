"""Predictor mentions in ClinVar SCV text, by quarter and by lab (H1c).

    python scripts/clinvar/scv_mentions.py --data-dir $DATA_DIR [--results results]

Needs the output of build_history.py (the missense universe) and
$DATA_DIR/clinvar/current/submission_summary.txt.gz. Writes CSV tables and
fig_scv_mentions.png into --results.
"""

import argparse
import time
from pathlib import Path

import polars as pl

from rdv.plots import fig_scv_mentions
from rdv.scv_text import (
    FLAGS, aggregate, by_lab, composition, decompose, scv_flags, standardized_share,
)

LABS = pl.col("criteria_provided") & (pl.col("method") == "clinical testing")
REF = ("2021-01-01", "2022-10-01")  # lab-mix reference: the 7 quarters before the first milestone
# (milestone, pre window, post window); post windows end before the incomplete current year.
WINDOWS = [
    ("pejaver_2022-12", ("2021-01-01", "2022-10-01"), ("2023-01-01", "2026-01-01")),
    ("alphamissense_2023-09", ("2022-01-01", "2023-07-01"), ("2023-10-01", "2026-01-01")),
]
FIG_TERMS = ["pp3_bp4", "revel", "alphamissense", "insilico_generic", "legacy_tools"]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", required=True, type=Path)
    ap.add_argument("--results", type=Path, default=Path("results"))
    a = ap.parse_args()

    derived = a.data_dir / "derived"
    flags_path = derived / "scv_flags.parquet"
    if not flags_path.exists():
        t = time.time()
        ids = pl.read_parquet(derived / "missense_universe.parquet", columns=["VariationID"])
        scv = scv_flags(a.data_dir / "clinvar/current/submission_summary.txt.gz", ids)
        scv.write_parquet(flags_path, compression="zstd")
        print(f"scv_flags: {scv.height:,} SCVs on missense variants ({time.time() - t:.0f}s)", flush=True)
    scv = pl.read_parquet(flags_path)

    a.results.mkdir(parents=True, exist_ok=True)
    for name, where in {"all": None, "labs": LABS}.items():
        aggregate(scv, where).write_csv(a.results / f"scv_mentions_{name}.csv")

    agg = aggregate(scv, LABS)
    comp = composition(scv, LABS)
    comp.write_csv(a.results / "scv_composition.csv")
    by_lab(scv, LABS).write_csv(a.results / "scv_by_lab.csv")
    std = {t: standardized_share(scv, t, REF, LABS) for t in FLAGS}
    pl.concat([d.with_columns(pl.lit(t).alias("term")) for t, d in std.items()]).write_csv(
        a.results / "scv_standardized.csv"
    )
    rows = [
        {"milestone": m, **decompose(scv, t, pre, post, LABS)}
        for m, pre, post in WINDOWS for t in FLAGS
    ]
    pl.DataFrame(rows).write_csv(a.results / "scv_decomposition.csv")
    fig_scv_mentions(agg, std, comp, a.results / "fig_scv_mentions.png")
    print("wrote tables and fig_scv_mentions.png to", a.results)


if __name__ == "__main__":
    main()
