"""Label history of ClinVar missense variants across quarterly snapshots (RQ1).

Two steps, both resumable:

1. `extract_snapshot`: one archived `variant_summary_YYYY-MM.txt.gz` -> one small
   parquet with the state (label, stars, last-evaluated date) of every variant in
   the current missense universe.
2. `derive`: all snapshot parquets -> per-variant first-confident-label table and
   a table of label transitions.

Variants are keyed by AlleleID, not VariationID: the early archives (2017) have no
VariationID column, and for single-nucleotide variants the two are 1:1.

Resolution is the snapshot interval. A variant "first confident in 2023-01" got its
P/LP or B/LB label between the 2022-10 and 2023-01 releases.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from .clinvar import (
    BLB, PLP, VUS, CONFLICT, OTHER,
    classify_significance, map_unique, parse_date, read_tsv, review_stars, snapshot_month,
)

SNAP_COLUMNS = ["AlleleID", "ClinicalSignificance", "ReviewStatus", "LastEvaluated"]
LABELS = [VUS, PLP, BLB, CONFLICT, OTHER]
LABEL_DTYPE = pl.Enum(LABELS)
CONFIDENT = [PLP, BLB]


def snapshot_path(out_dir: str | Path, month: str) -> Path:
    return Path(out_dir) / f"snapshot={month}.parquet"


def extract_snapshot(variant_summary: str | Path, allele_ids: pl.DataFrame) -> pl.DataFrame:
    """State of every variant in `allele_ids` (column `allele_id`) at one snapshot."""

    def keep(df: pl.DataFrame) -> pl.DataFrame:
        return (
            df.with_columns(pl.col("AlleleID").cast(pl.Int64, strict=False).alias("allele_id"))
            .join(allele_ids, on="allele_id", how="semi")
            .unique("allele_id", keep="first")  # GRCh37 + GRCh38 rows carry the same label
        )

    raw = read_tsv(variant_summary, SNAP_COLUMNS, keep).unique("allele_id", keep="first")
    return raw.select(
        pl.col("allele_id").cast(pl.Int64),
        pl.lit(snapshot_month(variant_summary)).alias("snapshot"),
        map_unique("ClinicalSignificance", classify_significance, pl.Utf8)
        .cast(LABEL_DTYPE)
        .alias("label"),
        map_unique("ReviewStatus", review_stars, pl.Int8).alias("stars"),
        parse_date("LastEvaluated").alias("last_evaluated"),
    ).sort("allele_id")


def derive(snapshot_dir: str | Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Return (first_confident, transitions) from a directory of snapshot parquets.

    first_confident: one row per variant that is P/LP or B/LB in any snapshot:
        first_seen         first snapshot in which the variant exists at all
        snapshot, label    first snapshot with a confident label, and that label
        stars, last_evaluated   review stars / evaluation date at that snapshot
        prior_label        label in the previous snapshot where the variant was
                           present; null if it entered ClinVar already confident
        left_censored      confident already in the earliest snapshot, so the true
                           date is earlier than we can see
        current_label      label in the latest snapshot (null if since removed)
    transitions: every change of label between consecutive snapshots in which the
        variant is present, for variants that are ever confident.
    """
    files = sorted(Path(snapshot_dir).glob("snapshot=*.parquet"))
    if not files:
        raise FileNotFoundError(f"no snapshot parquets in {snapshot_dir}")
    lf = pl.concat([pl.scan_parquet(f) for f in files])

    ever = lf.filter(pl.col("label").is_in(CONFIDENT)).select("allele_id").unique()
    h = (
        lf.join(ever, on="allele_id", how="semi")
        .sort("allele_id", "snapshot")
        .with_columns(
            pl.col("label").shift(1).over("allele_id").alias("prior_label"),
            pl.col("snapshot").shift(1).over("allele_id").alias("prior_snapshot"),
        )
        .collect()
    )
    earliest, latest = files[0].stem.split("=")[1], files[-1].stem.split("=")[1]

    confident = pl.col("label").is_in(CONFIDENT)
    first = (
        h.filter(confident)
        .group_by("allele_id")
        .agg(pl.all().sort_by("snapshot").first())
        .drop("prior_snapshot")
    )
    seen = h.group_by("allele_id").agg(pl.col("snapshot").min().alias("first_seen"))
    current = h.filter(pl.col("snapshot") == latest).select(
        "allele_id", pl.col("label").alias("current_label")
    )
    first_confident = (
        first.join(seen, on="allele_id")
        .join(current, on="allele_id", how="left")
        .with_columns((pl.col("snapshot") == earliest).alias("left_censored"))
        .select(
            "allele_id", "first_seen", "snapshot", "label", "stars", "last_evaluated",
            "prior_label", "left_censored", "current_label",
        )
        .sort("allele_id")
    )
    transitions = (
        h.filter(pl.col("prior_label").is_not_null() & (pl.col("label") != pl.col("prior_label")))
        .select("allele_id", "prior_snapshot", "snapshot", "prior_label", "label", "stars")
        .sort("allele_id", "snapshot")
    )
    return first_confident, transitions


def entry_route(first_confident: pl.DataFrame) -> pl.DataFrame:
    """Add `route`: how the variant reached its first confident label.

    new          entered ClinVar already P/LP or B/LB
    from_vus     was a VUS (or uncertain/other) in the previous snapshot
    from_conflict  was in conflict
    left_censored  already confident in the earliest snapshot (route unknown)
    """
    return first_confident.with_columns(
        pl.when(pl.col("left_censored"))
        .then(pl.lit("left_censored"))
        .when(pl.col("prior_label").is_null())
        .then(pl.lit("new"))
        .when(pl.col("prior_label") == CONFLICT)
        .then(pl.lit("from_conflict"))
        .otherwise(pl.lit("from_vus"))
        .alias("route")
    )


def quarter_of(snapshot: str) -> str:
    """'2023-01' -> '2023Q1' (ClinVar archive months are quarter starts here)."""
    y, m = snapshot.split("-")
    return f"{y}Q{(int(m) - 1) // 3 + 1}"


def attach_submitters(first_confident: pl.DataFrame, universe: pl.DataFrame, scv: pl.DataFrame) -> pl.DataFrame:
    """Add the likely labelling submitter to each first-confident variant.

    Among the variant's SCVs whose polarity (P/LP or B/LB) equals the entry label, take the one with
    the lowest accession number: accessions grow with time of first submission, so this is the first
    lab to assert the label the variant now carries. Approximate: a re-evaluated SCV keeps its
    accession. Variants with no polarity-matching SCV get a null submitter.

    `universe` needs AlleleID, VariationID, gene; `scv` needs variation_id, scv, label, submitter.
    """
    f = first_confident.with_columns(pl.col("label").cast(pl.Utf8)).join(
        universe.select(pl.col("AlleleID").alias("allele_id"), "VariationID", "gene"), on="allele_id"
    )
    s = scv.select(
        "variation_id", "scv", pl.col("label").cast(pl.Utf8), "submitter",
        pl.col("scv").str.slice(3).cast(pl.Int64).alias("acc"),
    )
    return (
        f.join(s, left_on=["VariationID", "label"], right_on=["variation_id", "label"], how="left")
        .sort("allele_id", "acc", nulls_last=True)
        .group_by("allele_id", maintain_order=True)
        .first()
        .drop("acc")
        .sort("allele_id")
    )
