"""Attach the likely labelling submitter to every first-confident variant.

    python scripts/clinvar/attach_submitters.py --data-dir $DATA_DIR

Needs build_history.py and scv_mentions.py outputs. Writes
$DATA_DIR/derived/first_confident_submitter.parquet and results/label_flow_by_submitter.csv
(entries per snapshot, submitter and label, the cells that the RQ1 protocol caps).
"""

import argparse
from pathlib import Path

import polars as pl

from rdv.history import attach_submitters


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", required=True, type=Path)
    ap.add_argument("--results", type=Path, default=Path("results"))
    a = ap.parse_args()
    d = a.data_dir / "derived"
    first = pl.read_parquet(d / "first_confident.parquet").filter(~pl.col("left_censored"))
    universe = pl.read_parquet(d / "missense_universe.parquet", columns=["AlleleID", "VariationID", "gene"])
    scv = pl.read_parquet(d / "scv_flags.parquet", columns=["variation_id", "scv", "label", "submitter"])
    out = attach_submitters(first, universe, scv)
    out.write_parquet(d / "first_confident_submitter.parquet", compression="zstd")
    cells = (
        out.with_columns(pl.col("submitter").fill_null("unknown"))
        .group_by("snapshot", "submitter", "label")
        .agg(pl.len().alias("n_variants"))
        .sort("n_variants", descending=True)
    )
    a.results.mkdir(parents=True, exist_ok=True)
    cells.write_csv(a.results / "label_flow_by_submitter.csv")
    print(f"{out.height:,} variants, {out['submitter'].null_count():,} without a matching SCV; "
          f"{cells.height:,} cells, {cells.filter(pl.col('n_variants') > 500).height} above 500")


if __name__ == "__main__":
    main()
