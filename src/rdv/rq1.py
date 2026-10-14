"""RQ1 analysis exactly as frozen in docs/protocol.md.

Everything that the protocol fixes is a constant or a default argument here. A change to any of
them is a deviation and belongs in the protocol's Deviations log.
"""

from __future__ import annotations

from dataclasses import dataclass

import polars as pl

from .metrics import StratifiedAUC, minimum_detectable_effect, per_gene_stats

# name -> (column in scores.parquet, sign) so that larger means more pathogenic after orientation
ORIENTATION = {
    "revel": ("revel", 1.0),
    "am": ("am_score", 1.0),
    "esm1b": ("esm1b_llr", -1.0),
    "grantham": ("grantham", 1.0),
    "blosum62": ("blosum62", -1.0),
}
CAP = 500
CAP_SENSITIVITY = (250, 1000, None)
REPLICATES = 2000
SEED = 20261012
PLACEBO_SHIFTS = (-24, -36)  # months


@dataclass(frozen=True)
class Comparison:
    name: str
    treated: str
    control: str
    pre: tuple[str, str]
    post: tuple[str, str]


# Primary tests (protocol section 4)
P1 = Comparison("P1_pejaver_2022-12", "revel", "esm1b", ("2021-01", "2022-10"), ("2023-01", "2024-10"))
P2 = Comparison("P2_alphamissense_2023-09", "am", "esm1b", ("2022-01", "2023-07"), ("2023-10", "2025-01"))
PRIMARY = (P1, P2)


def shift_snapshot(snapshot: str, months: int) -> str:
    """'2021-01' shifted by -24 months is '2019-01'."""
    y, m = map(int, snapshot.split("-"))
    total = y * 12 + (m - 1) + months
    return f"{total // 12:04d}-{total % 12 + 1:02d}"


def shifted(test: Comparison, months: int) -> Comparison:
    sh = lambda w: (shift_snapshot(w[0], months), shift_snapshot(w[1], months))  # noqa: E731
    return Comparison(f"{test.name}_shift{months:+d}m", test.treated, test.control, sh(test.pre), sh(test.post))


def build_cohort(first_confident_submitter: pl.DataFrame, scores: pl.DataFrame) -> pl.DataFrame:
    """One row per variant that entered the confident set after the first snapshot.

    Columns: allele_id, gene, snapshot, y (1 = P/LP), submitter ('unknown' if unattributed),
    route, stars, and one oriented column per predictor (null where unscored).
    """
    oriented = [(pl.col(col) * sign).alias(name) for name, (col, sign) in ORIENTATION.items() if col in scores.columns]
    sc = scores.select("allele_id", *oriented)
    return (
        first_confident_submitter.filter(pl.col("label").cast(pl.Utf8).is_in(["PLP", "BLB"]))
        .select(
            "allele_id", "gene", "snapshot",
            (pl.col("label").cast(pl.Utf8) == "PLP").cast(pl.Int8).alias("y"),
            pl.col("label").cast(pl.Utf8).alias("label"),
            pl.col("submitter").fill_null("unknown"),
            "route", "stars",
        )
        .join(sc, on="allele_id", how="left")
    )


def add_weights(cohort: pl.DataFrame, cap: int | None) -> pl.DataFrame:
    """Weight min(1, cap / n) for the n entries of each (snapshot, submitter, label) cell."""
    if cap is None:
        return cohort.with_columns(pl.lit(1.0).alias("w"))
    n = pl.len().over("snapshot", "submitter", "label")
    return cohort.with_columns(pl.min_horizontal(pl.lit(1.0), cap / n).alias("w"))


def run_did(
    cohort: pl.DataFrame,
    test: Comparison,
    cap: int | None = CAP,
    exclude_submitters: tuple[str, ...] = (),
    include_submitters: tuple[str, ...] | None = None,
    replicates: int = REPLICATES,
    seed: int = SEED,
    min_stars: int = 0,
    routes: tuple[str, ...] | None = None,
) -> dict:
    """One difference-in-differences with its bootstrap; returns a flat dict for a results table."""
    d = add_weights(cohort, cap)  # weights from the full cohort, before any filtering
    d = d.filter(pl.col(test.treated).is_not_null() & pl.col(test.control).is_not_null())  # common support
    if exclude_submitters:
        d = d.filter(~pl.col("submitter").is_in(list(exclude_submitters)))
    if include_submitters is not None:
        d = d.filter(pl.col("submitter").is_in(list(include_submitters)))
    if min_stars:
        d = d.filter(pl.col("stars") >= min_stars)
    if routes is not None:
        d = d.filter(pl.col("route").is_in(list(routes)))
    stats = per_gene_stats(d, [test.treated, test.control], {"pre": test.pre, "post": test.post})
    if stats.is_empty() or stats["window"].n_unique() < 2:
        return {"test": test.name, "n_variants": d.height, "estimate": None}
    m = StratifiedAUC(stats)
    r = m.did(test.treated, test.control, "pre", "post", replicates=replicates, seed=seed)
    return {
        "test": test.name, "treated": test.treated, "control": test.control, "cap": cap,
        "pre": f"{test.pre[0]}..{test.pre[1]}", "post": f"{test.post[0]}..{test.post[1]}",
        "estimate": r.estimate, "se": r.se, "lo": r.lo, "hi": r.hi, "mde": minimum_detectable_effect(r.se),
        "p_boot": two_sided_p(r.replicates),
        "treated_pre": r.treated_pre, "treated_post": r.treated_post,
        "control_pre": r.control_pre, "control_post": r.control_post,
        "pairs_pre": m.n_pairs("pre", test.treated), "pairs_post": m.n_pairs("post", test.treated),
        "genes_pre": m.n_genes("pre"), "genes_post": m.n_genes("post"),
        "variants_pre": d.filter((pl.col("snapshot") >= test.pre[0]) & (pl.col("snapshot") <= test.pre[1])).height,
        "variants_post": d.filter((pl.col("snapshot") >= test.post[0]) & (pl.col("snapshot") <= test.post[1])).height,
    }


def two_sided_p(replicates) -> float:
    """Bootstrap two-sided p-value for 'estimate = 0': twice the smaller tail, never below 1/(B+1)."""
    import numpy as np

    reps = np.asarray(replicates)
    tail = min((reps <= 0).mean(), (reps >= 0).mean())
    return float(max(2 * tail, 1.0 / (len(reps) + 1)))


def holm(pvalues: list[float]) -> list[float]:
    """Holm-adjusted p-values, in the input order."""
    order = sorted(range(len(pvalues)), key=lambda i: pvalues[i])
    adj, running = [0.0] * len(pvalues), 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (len(pvalues) - rank) * pvalues[i]))
        adj[i] = running
    return adj


def half_year_bins(first: str = "2019-01", last: str = "2026-10") -> dict[str, tuple[str, str]]:
    """Half-year bins of entry snapshot: Jan+Apr and Jul+Oct, keyed '2019H1' etc."""
    bins = {}
    y = int(first[:4])
    while f"{y}-01" <= last:
        bins[f"{y}H1"] = (f"{y}-01", f"{y}-04")
        bins[f"{y}H2"] = (f"{y}-07", f"{y}-10")
        y += 1
    return {k: v for k, v in bins.items() if v[0] <= last}


def event_study(
    cohort: pl.DataFrame,
    preds: list[str],
    diffs: tuple[tuple[str, str], ...] = (("revel", "esm1b"), ("am", "esm1b")),
    cap: int | None = CAP,
    replicates: int = REPLICATES,
    seed: int = SEED,
) -> pl.DataFrame:
    """WG-AUC per predictor and half-year bin, on common support of all `preds`.

    Long table: kind 'auc' (pred = name) or 'diff' (pred = 'a-b'), each with a gene-cluster
    bootstrap 95% interval, plus the number of within-gene pairs and genes in the bin.
    """
    d = add_weights(cohort, cap)
    for p in preds:
        d = d.filter(pl.col(p).is_not_null())
    bins = half_year_bins()
    m = StratifiedAUC(per_gene_stats(d, preds, bins))
    counts = m.bootstrap_counts(replicates, seed)
    rows = []
    for window in bins:
        if not any((window, p) in m._k for p in preds):
            continue
        base = {"bin": window, "pairs": m.n_pairs(window, preds[0]), "genes": m.n_genes(window)}
        for p in preds:
            reps = m.auc(window, p, counts)
            rows.append({**base, "kind": "auc", "pred": p, "value": float(m.auc(window, p)),
                         "lo": float(pl.Series(reps).quantile(0.025)), "hi": float(pl.Series(reps).quantile(0.975))})
        for a, b in diffs:
            if a in preds and b in preds:
                reps = m.auc(window, a, counts) - m.auc(window, b, counts)
                rows.append({**base, "kind": "diff", "pred": f"{a}-{b}",
                             "value": float(m.auc(window, a) - m.auc(window, b)),
                             "lo": float(pl.Series(reps).quantile(0.025)), "hi": float(pl.Series(reps).quantile(0.975))})
    return pl.DataFrame(rows)


# --- post-hoc additions, made after the primary results were seen (see protocol Deviations log) -----


def clean_boundary(test: Comparison) -> Comparison:
    """Same test with the first post snapshot dropped.

    An entry at snapshot S was labelled during the quarter before S, so the first post snapshot
    straddles the milestone (for AlphaMissense about 80% of it predates the release).
    """
    return Comparison(f"{test.name}_clean", test.treated, test.control, test.pre,
                      (shift_snapshot(test.post[0], 3), test.post[1]))


def attach_labelling_scv_flags(cohort: pl.DataFrame, first_confident_submitter: pl.DataFrame,
                               scv_flags: pl.DataFrame) -> pl.DataFrame:
    """Add the text flags (pp3_bp4, insilico_generic, ..., text_len) of each variant's labelling SCV."""
    flags = scv_flags.select("scv", "text_len", *[c for c in scv_flags.columns if c in
                                                   ("revel", "pp3_bp4", "insilico_generic", "alphamissense")])
    return cohort.join(first_confident_submitter.select("allele_id", "scv").join(flags, on="scv", how="left"),
                       on="allele_id", how="left")


def group_contrast(
    cohort: pl.DataFrame, flag: str, treated: str = "revel", control: str = "esm1b", cap: int | None = CAP,
    since: str = "2019-01", min_text: int = 30, replicates: int = REPLICATES, seed: int = SEED,
) -> dict:
    """Positive control: do labels whose own submission text cites something align differently?

    Group A: the labelling SCV's text has `flag`. Group B: its text (at least `min_text` characters)
    does not. Reports within-gene AUCs per group and the 'gap difference'
    [AUC(T, A) - AUC(T, B)] - [AUC(C, A) - AUC(C, B)] with a gene-cluster bootstrap interval.
    A design that can see predictor-shaped labels should find it positive for labels that cite the
    predictor family; case mix shifts all predictors together and cancels.
    """
    d = add_weights(cohort, cap).filter(pl.col("snapshot") >= since)
    preds = [c for c in ORIENTATION if c in d.columns]
    for p in preds:
        d = d.filter(pl.col(p).is_not_null())
    window = {"x": (since, "9999-12")}
    a = d.filter(pl.col(flag) == True)  # noqa: E712
    b = d.filter((pl.col(flag) == False) & (pl.col("text_len") >= min_text))  # noqa: E712
    stats = pl.concat([
        per_gene_stats(a, preds, window).with_columns(pl.lit("A").alias("window")),
        per_gene_stats(b, preds, window).with_columns(pl.lit("B").alias("window")),
    ])
    m = StratifiedAUC(stats)
    r = m.did(treated, control, "B", "A", replicates=replicates, seed=seed)
    out = {"flag": flag, "treated": treated, "control": control, "variants_A": a.height, "variants_B": b.height,
           "pairs_A": m.n_pairs("A", treated), "pairs_B": m.n_pairs("B", treated),
           "genes_A": m.n_genes("A"), "genes_B": m.n_genes("B"),
           "gap_difference": r.estimate, "se": r.se, "lo": r.lo, "hi": r.hi}
    for p in preds:
        out[f"auc_A_{p}"], out[f"auc_B_{p}"] = float(m.auc("A", p)), float(m.auc("B", p))
    return out
