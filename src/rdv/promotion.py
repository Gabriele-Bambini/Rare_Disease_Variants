"""VUS-promotion test (docs/protocol_promotion.md).

Between two consecutive ClinVar snapshots p and s, some variants that were VUS at p become P/LP (or
B/LB) at s. Each such *case* is compared with VUS in the same gene that stayed VUS over the same
interval (the *controls*, sampled from the gene's risk set). A weighted conditional logistic model
with REVEL, AlphaMissense, ESM1b, Grantham and a frequency proxy then asks which scores distinguish
promoted from unpromoted VUS, and whether that changes after a milestone (a period interaction).

Variants that moved VUS -> conflicting or other are neither case nor control; that is a stated limit
of the design, not something it fixes.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl

from . import clogit
from .rq1 import shift_snapshot

ML = ["revel", "am", "esm1b", "grantham"]  # first four covariates: flipped for the benign direction
COVARIATES = ML + ["log_af", "af_seen"]
M_CONTROLS = 20
CAP = 500
SEED = 20261014
DIRECTIONS = {"up": "PLP", "down": "BLB"}


@dataclass(frozen=True)
class Window:
    """Event snapshots s of the pre and post windows (inclusive). An event at s was labelled in (p, s]."""
    name: str
    focus: str  # the covariate whose change is the primary estimand
    pre: tuple[str, str]
    post: tuple[str, str]


# Primary windows. Q1: the Pejaver preprint (17 Mar 2022) falls in the interval ending 2022-04 and the
# journal version (20 Oct / 1 Dec 2022) in the intervals ending 2023-01, so those are dropped.
# Q2: AlphaMissense (19 Sep 2023) falls in the interval ending 2023-10, which is dropped.
Q1 = Window("Q1_pejaver_2022-12", "revel", ("2020-04", "2022-01"), ("2023-04", "2025-01"))
Q2 = Window("Q2_alphamissense_2023-09", "am", ("2022-01", "2023-07"), ("2024-01", "2025-01"))
PRIMARY = (Q1, Q2)


def shifted(w: Window, months: int) -> Window:
    sh = lambda x: (shift_snapshot(x[0], months), shift_snapshot(x[1], months))  # noqa: E731
    return Window(f"{w.name}_shift{months:+d}m", w.focus, sh(w.pre), sh(w.post))


def covariate_table(scores: pl.DataFrame, af: pl.DataFrame, universe: pl.DataFrame) -> pl.DataFrame:
    """z-scored covariates, oriented so that larger means more pathogenic (ESM1b is negated).

    Means and SDs come from every universe variant with all four scores. `log_af` is the base-10 log of
    the frequency proxy, z-scored over variants that have one, and 0 for those that do not;
    `af_seen` marks which is which. Rows missing any of the four scores are null.
    """
    t = (
        scores.select("allele_id", "VariationID", "revel", "am_score", pl.col("esm1b_llr").mul(-1).alias("esm1b"), "grantham")
        .rename({"am_score": "am"})
        .join(af, on="VariationID", how="left")
    )
    full = t.drop_nulls(ML)
    stats = {c: (full[c].mean(), full[c].std()) for c in ML}
    seen = full.filter(pl.col("af_proxy").is_not_null()).select(pl.col("af_proxy").clip(lower_bound=1e-6).log10().alias("l"))
    mu_l, sd_l = seen["l"].mean(), seen["l"].std()
    return t.select(
        "allele_id",
        *[((pl.col(c) - stats[c][0]) / stats[c][1]).alias(c) for c in ML],
        pl.when(pl.col("af_proxy").is_not_null())
        .then((pl.col("af_proxy").clip(lower_bound=1e-6).log10() - mu_l) / sd_l)
        .otherwise(0.0)
        .alias("log_af"),
        pl.col("af_proxy").is_not_null().cast(pl.Float64).alias("af_seen"),
    )


def snapshots(history_dir: str | Path) -> list[str]:
    return [p.stem.split("=")[1] for p in sorted(Path(history_dir).glob("snapshot=*.parquet"))]


def interval_events(history_dir, p: str, s: str, first_snapshot: pl.DataFrame) -> pl.DataFrame:
    """Outcome of every VUS at snapshot p that is still in ClinVar at s.

    Columns: allele_id, outcome in {'up', 'down', 'control', 'other'}.
      up / down   first confident label (P/LP, B/LB) appears at s: the variant has no earlier confident snapshot
      control     still VUS at s and never confident before s
      other       conflicting / other label at s, or confident before p (downgraded to VUS): not used
    """
    d = Path(history_dir)
    a = pl.read_parquet(d / f"snapshot={p}.parquet", columns=["allele_id", "label"]).filter(pl.col("label") == "VUS")
    b = pl.read_parquet(d / f"snapshot={s}.parquet", columns=["allele_id", "label"]).rename({"label": "label_s"})
    m = a.select("allele_id").join(b, on="allele_id").join(first_snapshot, on="allele_id", how="left")
    fs, ls = pl.col("first_snapshot"), pl.col("label_s").cast(pl.Utf8)
    return m.select(
        "allele_id",
        pl.when((ls == "PLP") & (fs == s)).then(pl.lit("up"))
        .when((ls == "BLB") & (fs == s)).then(pl.lit("down"))
        .when((ls == "VUS") & (fs.is_null() | (fs > s))).then(pl.lit("control"))
        .otherwise(pl.lit("other"))
        .alias("outcome"),
    )


def all_events(history_dir, first_snapshot, gene_of: pl.DataFrame, submitter_of: pl.DataFrame) -> pl.DataFrame:
    """Every promotion event (case) over every interval, with gene and attributed submitter. Labels only."""
    snaps = snapshots(history_dir)
    parts = []
    for p, s in zip(snaps[:-1], snaps[1:]):
        ev = interval_events(history_dir, p, s, first_snapshot).filter(pl.col("outcome").is_in(["up", "down"]))
        parts.append(ev.with_columns(pl.lit(s).alias("s")))
    return (
        pl.concat(parts)
        .join(gene_of, on="allele_id")
        .join(submitter_of, on="allele_id", how="left")
        .with_columns(pl.col("submitter").fill_null("unknown"))
    )


def case_weights(events: pl.DataFrame, cap: int | None) -> pl.DataFrame:
    """Add w = min(1, cap / n) for the n cases of each (s, submitter, direction) cell, from labels only."""
    if cap is None:
        return events.with_columns(pl.lit(1.0).alias("w"))
    n = pl.len().over("s", "submitter", "outcome")
    return events.with_columns(pl.min_horizontal(pl.lit(1.0), cap / n).alias("w"))


def wave_intervals(events: pl.DataFrame, share: float = 0.10) -> list[str]:
    """Intervals in the top `share` by total promotions, defined from event counts alone."""
    per = events.group_by("s").len().sort("len", descending=True)
    return per.head(max(1, int(round(share * per.height))))["s"].to_list()


@dataclass
class Strata:
    X: np.ndarray            # (S, 1 + m, P) oriented covariates; position 0 is the case
    valid: np.ndarray        # (S, 1 + m)
    meta: pl.DataFrame       # per stratum: allele_id, gene, s, outcome, submitter


def build_strata(
    history_dir, first_snapshot: pl.DataFrame, events: pl.DataFrame, gene_of: pl.DataFrame, cov: pl.DataFrame,
    direction: str, m: int = M_CONTROLS, seed: int = SEED, intervals: list[str] | None = None,
    allowed_controls=None,
) -> Strata:
    """Sample a risk set of up to `m` controls per case, within gene and interval, for one direction.

    `allowed_controls(case_row) -> sorted array of allele_ids, or None` restricts controls (used for the
    lab-matched analysis). Cases and controls need all covariates; others are dropped.
    """
    snaps = snapshots(history_dir)
    cov_ids = cov["allele_id"].to_numpy()
    order = np.argsort(cov_ids)
    cov_ids = cov_ids[order]
    cov_arr = cov.select(COVARIATES).to_numpy()[order].astype(float)
    sign = np.array([-1.0] * len(ML) + [1.0] * (len(COVARIATES) - len(ML))) if direction == "down" else np.ones(len(COVARIATES))
    complete = ~np.isnan(cov_arr[:, : len(ML)]).any(axis=1)

    def rows_of(ids):
        pos = np.searchsorted(cov_ids, ids)
        pos[pos >= len(cov_ids)] = 0
        ok = (cov_ids[pos] == ids) & complete[pos]
        return pos, ok

    rng = np.random.default_rng([seed, 1 if direction == "up" else 2])
    idx_rows, metas = [], []
    for p, s in zip(snaps[:-1], snaps[1:]):
        if intervals is not None and s not in intervals:
            continue
        ev = events.filter((pl.col("s") == s) & (pl.col("outcome") == direction))
        if ev.is_empty():
            continue
        io = interval_events(history_dir, p, s, first_snapshot).filter(pl.col("outcome") == "control")
        genes = ev["gene"].unique()
        pool = io.join(gene_of, on="allele_id").filter(pl.col("gene").is_in(genes.implode()))
        by_gene = {g: np.asarray(ids, dtype=np.int64) for g, ids in pool.group_by("gene").agg(pl.col("allele_id")).iter_rows()}
        for row in ev.iter_rows(named=True):
            ctrl = by_gene.get(row["gene"], np.empty(0, dtype=np.int64))
            if allowed_controls is not None:
                keep = allowed_controls(row)  # sorted array of allele ids, or None for no restriction
                if keep is not None:
                    at = np.searchsorted(keep, ctrl)
                    at[at >= len(keep)] = 0
                    ctrl = ctrl[(len(keep) > 0) & (keep[at] == ctrl)] if len(keep) else ctrl[:0]
            _, ok = rows_of(ctrl)
            ctrl = ctrl[ok]
            _, case_ok = rows_of(np.array([row["allele_id"]]))
            if not case_ok[0] or len(ctrl) == 0:
                continue
            chosen = ctrl if len(ctrl) <= m else rng.choice(ctrl, size=m, replace=False)
            idx_rows.append(np.concatenate([[row["allele_id"]], chosen]))
            metas.append((row["allele_id"], row["gene"], s, direction, row["submitter"]))
    S = len(idx_rows)
    P = len(COVARIATES)
    X = np.zeros((S, 1 + m, P))
    valid = np.zeros((S, 1 + m), dtype=bool)
    for k, ids in enumerate(idx_rows):
        pos, _ = rows_of(ids)
        X[k, : len(ids)] = cov_arr[pos] * sign
        valid[k, : len(ids)] = True
    meta = pl.DataFrame(metas, schema=["allele_id", "gene", "s", "outcome", "submitter"], orient="row")
    return Strata(X, valid, meta)


def fit_window(strata: Strata, window: Window, events: pl.DataFrame, cap: int | None = CAP,
               exclude_s: tuple[str, ...] = (), exclude_submitters: tuple[str, ...] = ()) -> dict:
    """Interaction fit for one window: coefficients per covariate in pre and post, and their change."""
    meta = strata.meta
    s = meta["s"]
    in_pre = (s >= window.pre[0]) & (s <= window.pre[1])
    in_post = (s >= window.post[0]) & (s <= window.post[1])
    keep = (in_pre | in_post)
    if exclude_s:
        keep = keep & ~s.is_in(pl.Series(list(exclude_s), dtype=pl.Utf8).implode())
    if exclude_submitters:
        keep = keep & ~meta["submitter"].is_in(pl.Series(list(exclude_submitters), dtype=pl.Utf8).implode())
    keep_np = keep.to_numpy()
    sub = meta.filter(keep)
    w = _weights(strata, events, cap)[keep_np]
    post = (sub["s"] >= window.post[0]).to_numpy().astype(float)
    X = strata.X[keep_np]
    Xi = np.concatenate([X, X * post[:, None, None]], axis=2)
    _, cluster = np.unique(sub["gene"].to_numpy(), return_inverse=True)
    f = clogit.fit(Xi, strata.valid[keep_np], w, cluster)
    P = len(COVARIATES)
    out = {"test": window.name, "n_strata": f.n_strata, "n_genes": f.n_clusters,
           "n_pre": int((post == 0).sum()), "n_post": int((post == 1).sum()), "converged": f.converged}
    for k, name in enumerate(COVARIATES):
        d, se = f.beta[P + k], f.se[P + k]
        out[f"pre_{name}"] = float(f.beta[k])
        out[f"post_{name}"] = float(f.beta[k] + d)
        out[f"delta_{name}"], out[f"se_{name}"] = float(d), float(se)
        out[f"lo_{name}"], out[f"hi_{name}"] = float(d - 1.96 * se), float(d + 1.96 * se)
    # change in the summed reliance on the three machine-learning predictors
    c = np.zeros(2 * P)
    c[[P + COVARIATES.index(x) for x in ("revel", "am", "esm1b")]] = 1.0
    out["delta_sum_ml"], out["se_sum_ml"] = float(c @ f.beta), float(np.sqrt(c @ f.cov @ c))
    return out


def _weights(strata: Strata, events: pl.DataFrame, cap: int | None) -> np.ndarray:
    """Per-stratum weights from the full event table, so a cap never depends on which cases could be scored."""
    wt = case_weights(events, cap).select("allele_id", "s", "outcome", pl.col("w").alias("w_new"))
    return strata.meta.join(wt, on=["allele_id", "s", "outcome"], how="left")["w_new"].to_numpy()


def pool_information(history_dir, first_snapshot, events, gene_of, cov, direction: str, window: Window,
                     cap: int | None = CAP) -> dict:
    """Expected information at beta = 0, from risk sets alone (before any outcome is modelled).

    Under no effect every member of a risk set is equally likely to be the case, so each risk set
    contributes (case weight) x (covariance of the covariates over its members, here the unpromoted VUS
    of the gene, so that no case's score enters the number). Summing over a
    window gives the model-based information; Var(change) = inv(I_pre) + inv(I_post). Returns the
    minimum detectable change (2.8 SE) per covariate. This is a lower bound: it ignores gene clustering.
    """
    snaps = snapshots(history_dir)
    cov_ids = cov["allele_id"].to_numpy()
    order = np.argsort(cov_ids)
    cov_ids, cov_arr = cov_ids[order], cov.select(COVARIATES).to_numpy()[order].astype(float)
    complete = ~np.isnan(cov_arr[:, : len(ML)]).any(axis=1)
    sign = np.array([-1.0] * len(ML) + [1.0] * (len(COVARIATES) - len(ML))) if direction == "down" else np.ones(len(COVARIATES))
    weighted = case_weights(events, cap).filter(pl.col("outcome") == direction)

    def rows_of(ids):
        pos = np.searchsorted(cov_ids, ids)
        pos[pos >= len(cov_ids)] = 0
        return pos, (cov_ids[pos] == ids) & complete[pos]

    def info(win: tuple[str, str]) -> np.ndarray:
        P = len(COVARIATES)
        total = np.zeros((P, P))
        for p, s in zip(snaps[:-1], snaps[1:]):
            if not (win[0] <= s <= win[1]):
                continue
            ev = weighted.filter(pl.col("s") == s)
            if ev.is_empty():
                continue
            ctrl = interval_events(history_dir, p, s, first_snapshot).filter(pl.col("outcome") == "control")
            ctrl = ctrl.join(gene_of, on="allele_id").filter(pl.col("gene").is_in(ev["gene"].unique().implode()))
            by_gene = {g: np.asarray(i, dtype=np.int64) for g, i in ctrl.group_by("gene").agg(pl.col("allele_id")).iter_rows()}
            for g, rows in ev.group_by("gene"):
                cases = rows["allele_id"].to_numpy()
                # outcome-blind: the covariance comes from the unpromoted VUS only, so no case's score enters it
                members = by_gene.get(g[0], np.empty(0, dtype=np.int64))
                pos, ok = rows_of(members)
                _, case_ok = rows_of(cases)
                if ok.sum() < 2 or not case_ok.any():
                    continue
                X = cov_arr[pos[ok]] * sign
                total += rows["w"].to_numpy()[case_ok].sum() * np.cov(X.T, bias=True)
        return total

    i_pre, i_post = info(window.pre), info(window.post)
    active = (np.diag(i_pre) > 1e-12) & (np.diag(i_post) > 1e-12)  # a covariate with no variation has no information
    var = np.full(len(COVARIATES), np.nan)
    sub = np.ix_(active, active)
    var[active] = np.diag(np.linalg.inv(i_pre[sub]) + np.linalg.inv(i_post[sub]))
    return {f"mde_{name}": float(2.8 * np.sqrt(v)) for name, v in zip(COVARIATES, var)}
