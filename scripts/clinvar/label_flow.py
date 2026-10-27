"""How variants reach a confident ClinVar label, per quarter (descriptive).

    python scripts/clinvar/label_flow.py --data-dir $DATA_DIR [--results results]
"""

import argparse
from pathlib import Path

import polars as pl

from rdv.plots import fig_label_flow


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", required=True, type=Path)
    ap.add_argument("--results", type=Path, default=Path("results"))
    a = ap.parse_args()
    first = pl.read_parquet(a.data_dir / "derived/first_confident.parquet")
    a.results.mkdir(parents=True, exist_ok=True)
    (
        first.filter(~pl.col("left_censored"))
        .group_by("snapshot", "label", "route")
        .agg(pl.len().alias("n_variants"), (pl.col("stars") >= 2).sum().alias("n_two_star_plus"))
        .sort("snapshot", "label", "route")
        .write_csv(a.results / "label_flow_by_quarter.csv")
    )
    universe = pl.read_parquet(a.data_dir / "derived/missense_universe.parquet", columns=["label", "stars"])
    trans = pl.read_parquet(a.data_dir / "derived/transitions.parquet")
    new = first.filter(~pl.col("left_censored"))
    metrics = {"universe_missense_variants": universe.height}
    metrics |= {f"universe_{k}": v for k, v in universe["label"].value_counts().iter_rows()}
    metrics["ever_confident_variants"] = first.height
    metrics["left_censored_at_first_snapshot"] = first.height - new.height
    metrics["first_confident_two_star_plus"] = int((new["stars"] >= 2).sum())
    for label, route, n in new.group_by("label", "route").len().sort("label", "route").iter_rows():
        metrics[f"first_confident_{label}_{route}"] = n
    flips = trans.filter(
        pl.col("prior_label").is_in(["PLP", "BLB"]) & pl.col("label").is_in(["PLP", "BLB"])
    )
    metrics["direct_plp_blb_flips"] = flips.height
    pl.DataFrame({"metric": list(metrics), "value": [int(v) for v in metrics.values()]}).write_csv(
        a.results / "summary.csv"
    )
    fig_label_flow(first, a.results / "fig_label_flow.png")
    print("wrote label_flow_by_quarter.csv, summary.csv and fig_label_flow.png to", a.results)


if __name__ == "__main__":
    main()
