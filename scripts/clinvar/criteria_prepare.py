"""Prepare the per-submission criteria tables for the computational-evidence audit (docs/protocol_criteria.md).

    python scripts/clinvar/criteria_prepare.py --data-dir $DATA_DIR

One pass over submission_summary, restricted to submissions on the missense universe. Writes, under
$DATA_DIR/derived:

  criteria_scv.parquet     one row per submission whose text lists at least one non-negated ACMG/AMP code,
                           with the parsed codes, the derived class under four readings of the rules and
                           the decisiveness flags (rdv.criteria.describe); the text itself is kept for audit.
  scv_denominator.parquet  counts of ALL submissions on the universe by stated class, lab, evaluation year,
                           review status and collection method, so the audit can say what share of
                           submissions it can see at all.

Computes no estimate.
"""

import argparse
import time
from pathlib import Path

import polars as pl

from rdv.clinvar import classify_significance, iter_tsv_chunks, map_unique, parse_date
from rdv.criteria import describe

COLS = ["VariationID", "ClinicalSignificance", "DateLastEvaluated", "Description", "ExplanationOfInterpretation",
        "ReviewStatus", "CollectionMethod", "Submitter", "SCV"]
# A superset filter with no word edges (the polars regex engine has no lookaround); the Python parser decides.
ANY_CODE = r"(?i)(?:pvs1|ps[1-4]|pm[1-6]|pp[1-5]|ba1|bs[1-4]|bp[1-7])"
TEXT_KEEP = 2000  # characters of submitter text stored per row


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", required=True, type=Path)
    a = ap.parse_args()
    d = a.data_dir / "derived"
    ids = pl.read_parquet(d / "missense_universe.parquet", columns=["VariationID"])

    t = time.time()
    candidates, counts = [], []
    for df in iter_tsv_chunks(a.data_dir / "clinvar/current/submission_summary.txt.gz", COLS, chunk_lines=200_000):
        df = df.with_columns(pl.col("VariationID").cast(pl.Int64, strict=False)).join(ids, on="VariationID", how="semi")
        text = pl.concat_str(pl.col("Description").fill_null(""), pl.lit(" | "), pl.col("ExplanationOfInterpretation").fill_null(""))
        df = df.with_columns(
            text.alias("text"),
            map_unique("ClinicalSignificance", classify_significance, pl.String).alias("class"),
            parse_date("DateLastEvaluated").dt.year().alias("year"),
        )
        counts.append(df.group_by("class", "Submitter", "year", "ReviewStatus", "CollectionMethod").len())
        candidates.append(df.filter(pl.col("text").str.count_matches(ANY_CODE) >= 1).select(
            "VariationID", "ClinicalSignificance", "DateLastEvaluated", "ReviewStatus", "CollectionMethod", "Submitter", "SCV", "text"))
    denominator = (pl.concat(counts).group_by("class", "Submitter", "year", "ReviewStatus", "CollectionMethod")
                   .agg(pl.col("len").sum().alias("n")))
    denominator.write_parquet(d / "scv_denominator.parquet", compression="zstd")
    raw = pl.concat(candidates)
    print(f"{denominator['n'].sum():,} submissions on the universe; {raw.height:,} mention a code ({time.time() - t:.0f}s)", flush=True)

    rows = []
    for r in raw.iter_rows(named=True):
        info = describe(r["text"], r["ClinicalSignificance"])
        if info is None:
            continue
        rows.append({"variation_id": r["VariationID"], "scv": r["SCV"].split(".")[0], "lab": r["Submitter"],
                     "review": r["ReviewStatus"], "method": r["CollectionMethod"], "date": r["DateLastEvaluated"],
                     **info, "text": r["text"][:TEXT_KEEP]})
    out = pl.DataFrame(rows)
    out.write_parquet(d / "criteria_scv.parquet", compression="zstd")
    print(f"{out.height:,} submissions list at least one non-negated code, on {out['variation_id'].n_unique():,} variants "
          f"from {out['lab'].n_unique()} labs ({time.time() - t:.0f}s)")


if __name__ == "__main__":
    main()
