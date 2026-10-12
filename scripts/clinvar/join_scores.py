"""Attach predictor scores to the missense universe and report coverage.

    python scripts/clinvar/join_scores.py --data-dir $DATA_DIR [--results results]

Writes $DATA_DIR/derived/scores.parquet (one row per universe variant, null where a predictor has
no score) and results/score_coverage.csv. Each score family is skipped, with a message, if its
download is missing; rerun after fetching it.
"""

import argparse
import time
from pathlib import Path

import polars as pl

from rdv.aa import THREE_TO_ONE, substitution_scores
from rdv.scores import alphamissense, mane_table


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", required=True, type=Path)
    ap.add_argument("--results", type=Path, default=Path("results"))
    a = ap.parse_args()
    d = a.data_dir
    universe = pl.read_parquet(d / "derived/missense_universe.parquet")
    out = universe.select(pl.col("AlleleID").alias("allele_id"), "VariationID", "gene", "hgnc_id")

    # chemistry-only baselines: need nothing but the substitution
    sub = universe.select(
        pl.col("ref_aa").replace_strict(THREE_TO_ONE, default=None).alias("r"),
        pl.col("alt_aa").replace_strict(THREE_TO_ONE, default=None).alias("a"),
    )
    chem = substitution_scores(sub["r"], sub["a"])
    out = out.with_columns(chem["grantham"], chem["blosum62"])

    mane_files = sorted((d / "mane").glob("MANE.GRCh38.v*.summary.txt.gz"))
    am_path = d / "alphamissense/AlphaMissense_hg38.tsv.gz"
    if mane_files and am_path.exists():
        t = time.time()
        am = alphamissense(universe, am_path, mane_table(mane_files[-1]))
        am.write_parquet(d / "derived/am_scores.parquet", compression="zstd")
        print(f"alphamissense: {am.height:,} of {universe.height:,} variants ({time.time() - t:.0f}s)", flush=True)
        out = out.join(am, on="allele_id", how="left")
    else:
        print("skipping AlphaMissense: run download_all.sh mane alphamissense")

    out.write_parquet(d / "derived/scores.parquet", compression="zstd")
    cov = pl.DataFrame({
        "column": [c for c in out.columns if c not in ("allele_id", "VariationID", "gene", "hgnc_id")],
    }).with_columns(
        pl.col("column").map_elements(lambda c: out.height - out[c].null_count(), return_dtype=pl.Int64).alias("n_scored"),
    ).with_columns((pl.col("n_scored") / out.height).alias("fraction"))
    a.results.mkdir(parents=True, exist_ok=True)
    cov.write_csv(a.results / "score_coverage.csv")
    print(cov)


if __name__ == "__main__":
    main()
