"""What do submitters write? Predictor mentions in ClinVar SCV free text (H1c).

`submission_summary.txt.gz` holds one row per submission (SCV) with the lab's free-text
`Description` and the date it last evaluated the variant. We flag which predictors and
ACMG codes each text mentions, then aggregate by quarter of `DateLastEvaluated`.

Caveats that matter when reading the series:
* The current file only keeps each SCV's latest text. An SCV re-evaluated in 2025 lives
  in 2025 even if its first text was written in 2019.
* Labs write from templates, so one lab changing its template moves the whole series.
  `aggregate` therefore also reports how many distinct submitters mention each term and
  how much of the volume the largest one supplies.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from .clinvar import classify_significance, map_unique, parse_date, read_tsv

SCV_COLUMNS = [
    "VariationID", "ClinicalSignificance", "DateLastEvaluated", "Description",
    "ExplanationOfInterpretation", "ReviewStatus", "CollectionMethod", "Submitter", "SCV",
]

# name -> (regex, role). Roles: "endorsed" = predictors/codes calibrated for ACMG PP3/BP4;
# "control" = evidence types the PP3/BP4 calibration did not touch; "context" = shares
# that normalise the rest.
# Word edges are written as (^|[^a-z0-9]) ... ([^a-z0-9]|$) rather than \b: criteria are often glued to their
# strength with an underscore ("PP3_Moderate"), and underscore is a word character, so \b sees no edge there.
# (The regex engine behind polars has no lookaround, hence the consuming form.)
_L, _R = r"(?:^|[^a-z0-9])", r"(?:$|[^a-z0-9])"

PATTERNS: dict[str, tuple[str, str]] = {
    "pp3_bp4": (rf"(?i){_L}(?:pp3|bp4){_R}", "endorsed"),
    "revel": (rf"(?i){_L}revel{_R}", "endorsed"),
    "alphamissense": (r"(?i)alpha[\s-]?missense", "endorsed"),
    "cadd": (rf"(?i){_L}cadd{_R}", "endorsed"),
    "bayesdel": (r"(?i)bayes[\s-]?del", "endorsed"),
    "legacy_tools": (
        rf"(?i){_L}(?:sift|polyphen|mutation[\s-]?taster|provean|fathmm|metasvm|metalr){_R}",
        "control",
    ),
    "conservation": (rf"(?i)conserv|{_L}(?:phylop|phastcons|gerp){_R}", "control"),
    "population": (
        rf"(?i)gnomad|{_L}exac{_R}|1000 genomes|population (?:database|frequency)|allele frequency",
        "control",
    ),
    "esm": (rf"(?i){_L}esm[\s-]?(?:1b|2)?{_R}", "control"),
    # The largest labs (e.g. Ambry) say "predicted to be deleterious by in silico analysis"
    # without naming a tool, so tool names alone undercount computational evidence.
    "insilico_generic": (
        r"(?i)in[\s-]?silico|computational (?:evidence|prediction|analysis|tool|model|algorithm|method)"
        r"|prediction (?:tool|algorithm|software|program)"
        r"|predicted (?:to be|as) (?:deleterious|damaging|benign|tolerated|pathogenic|disease)",
        "context",
    ),
    "acmg_code": (rf"(?i){_L}(?:pvs1|ps[1-4]|pm[1-6]|pp[1-5]|ba1|bs[1-4]|bp[1-7]){_R}", "context"),
}
FLAGS = list(PATTERNS)
MIN_TEXT = 30  # characters; shorter descriptions carry no evidence statement

# Milestones from docs/brainstorm.md section 2. The exact dates are to be verified in week 1.
MILESTONES = {
    "2022-12-01": "Pejaver et al. PP3/BP4 calibration",
    "2023-09-19": "AlphaMissense release",
}


def scv_flags(submission_summary: str | Path, variation_ids: pl.DataFrame) -> pl.DataFrame:
    """One row per SCV on a variant in `variation_ids` (column `VariationID`, Int64)."""
    text = pl.concat_str(
        pl.col("Description").fill_null(""), pl.lit(" "),
        pl.col("ExplanationOfInterpretation").fill_null(""),
    )

    def keep(df: pl.DataFrame) -> pl.DataFrame:
        df = df.with_columns(pl.col("VariationID").cast(pl.Int64, strict=False)).join(
            variation_ids, on="VariationID", how="semi"
        )
        return df.select(
            "VariationID",
            pl.col("SCV").str.extract(r"^(SCV\d+)").alias("scv"),
            pl.col("ClinicalSignificance"),
            parse_date("DateLastEvaluated").alias("date"),
            pl.col("ReviewStatus"),
            pl.col("CollectionMethod").alias("method"),
            pl.col("Submitter").alias("submitter"),
            text.str.len_chars().alias("text_len"),
            *[text.str.contains(rx).alias(name) for name, (rx, _) in PATTERNS.items()],
        )

    df = read_tsv(submission_summary, SCV_COLUMNS, keep, chunk_lines=200_000)
    return df.select(
        pl.col("VariationID").alias("variation_id"),
        "scv",
        map_unique("ClinicalSignificance", classify_significance, pl.Utf8).alias("label"),
        "date",
        pl.col("ReviewStatus").str.starts_with("criteria provided").alias("criteria_provided"),
        "method",
        "submitter",
        "text_len",
        *FLAGS,
    )


def _wilson(k: pl.Expr, n: pl.Expr, z: float = 1.96) -> tuple[pl.Expr, pl.Expr]:
    p = k / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    half = z * ((p * (1 - p) / n + z**2 / (4 * n**2)).sqrt()) / denom
    return centre - half, centre + half


def aggregate(scv: pl.DataFrame, where: pl.Expr | None = None, min_text: int = MIN_TEXT) -> pl.DataFrame:
    """Share of SCVs-with-text mentioning each flag, by quarter, long format.

    Columns: quarter, term, n_text (denominator), k, share, lo, hi (95% Wilson),
    n_submitters (distinct submitters that mention the term), top_submitter_share
    (fraction of the term's mentions supplied by the largest single submitter).
    """
    d = scv.filter(pl.col("date").is_not_null() & (pl.col("text_len") >= min_text))
    if where is not None:
        d = d.filter(where)
    d = d.with_columns(pl.col("date").dt.truncate("3mo").alias("quarter"))
    denom = d.group_by("quarter").agg(pl.len().alias("n_text"))
    parts = []
    for term in FLAGS:
        hit = d.filter(pl.col(term))
        per_sub = hit.group_by("quarter", "submitter").agg(pl.len().alias("m"))
        by_q = per_sub.group_by("quarter").agg(
            pl.col("m").sum().alias("k"),
            pl.col("submitter").n_unique().alias("n_submitters"),
            (pl.col("m").max() / pl.col("m").sum()).alias("top_submitter_share"),
        )
        parts.append(by_q.with_columns(pl.lit(term).alias("term")))
    long = pl.concat(parts)
    grid = denom.join(pl.DataFrame({"term": FLAGS}), how="cross")
    out = (
        grid.join(long, on=["quarter", "term"], how="left")
        .with_columns(
            pl.col("k").fill_null(0),
            pl.col("n_submitters").fill_null(0),
        )
    )
    lo, hi = _wilson(pl.col("k").cast(pl.Float64), pl.col("n_text").cast(pl.Float64))
    return (
        out.with_columns(
            (pl.col("k") / pl.col("n_text")).alias("share"), lo.alias("lo"), hi.alias("hi")
        )
        .select("quarter", "term", "n_text", "k", "share", "lo", "hi", "n_submitters",
                "top_submitter_share")
        .sort("term", "quarter")
    )


# --- separating "labs changed what they write" from "different labs dominate" -----


def _text_bearing(scv: pl.DataFrame, where: pl.Expr | None, min_text: int = MIN_TEXT) -> pl.DataFrame:
    d = scv.filter(pl.col("date").is_not_null() & (pl.col("text_len") >= min_text))
    return d.filter(where) if where is not None else d


def _in_window(window: tuple[str, str]) -> pl.Expr:
    lo, hi = (pl.lit(w).str.to_date() for w in window)
    return (pl.col("date") >= lo) & (pl.col("date") < hi)


def composition(scv: pl.DataFrame, where: pl.Expr | None = None, top: int = 3) -> pl.DataFrame:
    """Share of text-bearing SCVs per quarter supplied by the `top` largest labs vs the rest."""
    d = _text_bearing(scv, where)
    names = d.group_by("submitter").len().sort("len", descending=True).head(top)["submitter"].to_list()
    d = d.with_columns(
        pl.when(pl.col("submitter").is_in(names))
        .then(pl.col("submitter").str.split(",").list.first().str.replace(r"\s*\(.*$", ""))
        .otherwise(pl.lit("All other labs"))
        .alias("group"),
        pl.col("date").dt.truncate("3mo").alias("quarter"),
    )
    out = d.group_by("quarter", "group").agg(pl.len().alias("n"))
    return out.with_columns((pl.col("n") / pl.col("n").sum().over("quarter")).alias("share")).sort(
        "quarter", "group"
    )


def standardized_share(
    scv: pl.DataFrame,
    term: str,
    ref: tuple[str, str],
    where: pl.Expr | None = None,
    min_n: int = 20,
) -> pl.DataFrame:
    """Quarterly mention share with every lab held at its weight in the `ref` window.

    Labs absent from a quarter (fewer than `min_n` texts) are dropped and the rest are
    re-weighted; `coverage` is the fraction of the reference weight still represented.
    """
    d = _text_bearing(scv, where).with_columns(pl.col("date").dt.truncate("3mo").alias("quarter"))
    ref_w = d.filter(_in_window(ref)).group_by("submitter").agg(pl.len().alias("w"))
    per = (
        d.group_by("quarter", "submitter")
        .agg(pl.col(term).mean().alias("s"), pl.len().alias("n"))
        .filter(pl.col("n") >= min_n)
        .join(ref_w, on="submitter")
    )
    total = ref_w["w"].sum()
    return (
        per.group_by("quarter")
        .agg(
            ((pl.col("s") * pl.col("w")).sum() / pl.col("w").sum()).alias("share_std"),
            (pl.col("w").sum() / total).alias("coverage"),
            pl.len().alias("n_labs"),
        )
        .sort("quarter")
    )


def decompose(
    scv: pl.DataFrame,
    term: str,
    pre: tuple[str, str],
    post: tuple[str, str],
    where: pl.Expr | None = None,
    min_n: int = 50,
) -> dict:
    """Split the pre->post change in a mention share into within-lab and between-lab parts.

    Over labs with >= `min_n` texts in both windows, with w = lab's share of volume and
    s = its mention rate (symmetric Kitagawa decomposition, exact):
        change = sum( (w_pre+w_post)/2 * (s_post-s_pre) )     within: labs changed wording
               + sum( (s_pre+s_post)/2 * (w_post-w_pre) )     between: the mix changed
    `pooled_*` use every text in the window, including labs that enter or leave, so
        pooled_delta = within_lab + between_lab + entry_exit
    where `entry_exit` is whatever the continuing labs do not explain: labs that appear
    or disappear between the windows, or are too small to count.
    """
    d = _text_bearing(scv, where)

    def window(win):
        return d.filter(_in_window(win)).group_by("submitter").agg(
            pl.len().alias("n"), pl.col(term).mean().alias("s")
        )

    a, b = window(pre), window(post)
    def pooled(w):
        return float((w["n"] * w["s"]).sum() / w["n"].sum())

    c = a.join(b, on="submitter", suffix="_post").filter(
        (pl.col("n") >= min_n) & (pl.col("n_post") >= min_n)
    )
    w0, w1 = c["n"] / c["n"].sum(), c["n_post"] / c["n_post"].sum()
    s0, s1 = c["s"], c["s_post"]
    within = float(((w0 + w1) / 2 * (s1 - s0)).sum())
    between = float(((s0 + s1) / 2 * (w1 - w0)).sum())
    cont_pre, cont_post = float((w0 * s0).sum()), float((w1 * s1).sum())
    pooled_pre, pooled_post = pooled(a), pooled(b)
    return {
        "term": term,
        "pooled_pre": pooled_pre,
        "pooled_post": pooled_post,
        "pooled_delta": pooled_post - pooled_pre,
        "continuing_pre": cont_pre,
        "continuing_post": cont_post,
        "within_lab": within,
        "between_lab": between,
        "entry_exit": (pooled_post - pooled_pre) - (cont_post - cont_pre),
        "n_labs": c.height,
        "coverage_pre": float(c["n"].sum() / a["n"].sum()),
        "coverage_post": float(c["n_post"].sum() / b["n"].sum()),
    }


def by_lab(scv: pl.DataFrame, where: pl.Expr | None = None, top: int = 15) -> pl.DataFrame:
    """Mention rate per year for the `top` labs by volume, wide in terms."""
    d = _text_bearing(scv, where).with_columns(pl.col("date").dt.year().alias("year"))
    names = d.group_by("submitter").len().sort("len", descending=True).head(top)["submitter"]
    return (
        d.filter(pl.col("submitter").is_in(names))
        .group_by("submitter", "year")
        .agg(pl.len().alias("n_text"), *[pl.col(t).mean().alias(t) for t in FLAGS])
        .sort("submitter", "year")
    )
