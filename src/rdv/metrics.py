"""Within-gene weighted AUC and its gene-cluster bootstrap (docs/protocol.md, sections 3 to 6).

The statistic is a ratio of sums over genes, so each gene is reduced once to two numbers per
(window, predictor): the weighted count of (P, B) pairs the predictor orders correctly (`num`,
ties count half) and the weighted count of all pairs (`den`). Any resampling of genes is then a
matrix product, which makes thousands of bootstrap replicates cheap.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import polars as pl


def pair_stats(y: np.ndarray, w: np.ndarray, s: np.ndarray) -> tuple[float, float]:
    """(num, den) for one stratum: y bool (True = P/LP), weights w, scores s (higher = more pathogenic)."""
    y = np.asarray(y, dtype=bool)
    if y.all() or not y.any():
        return 0.0, 0.0
    order = np.argsort(s, kind="stable")
    s, y, w = s[order], y[order], w[order]
    _, start = np.unique(s, return_index=True)  # one group per distinct score, so ties are exact
    w_neg = np.add.reduceat(np.where(y, 0.0, w), start)
    w_pos = np.add.reduceat(np.where(y, w, 0.0), start)
    below = np.cumsum(w_neg) - w_neg
    return float(np.sum(w_pos * (below + 0.5 * w_neg))), float(w_pos.sum() * w_neg.sum())


def per_gene_stats(
    df: pl.DataFrame, preds: list[str], windows: dict[str, tuple[str, str]]
) -> pl.DataFrame:
    """Per-gene (num, den) for every predictor in every window.

    `df` needs gene, y (0/1), w, snapshot ('YYYY-MM') and one column per predictor, already oriented
    so that higher means more pathogenic and already restricted to common support. A window is an
    inclusive range of entry snapshots.
    """
    rows = []
    for window, (lo, hi) in windows.items():
        d = df.filter((pl.col("snapshot") >= lo) & (pl.col("snapshot") <= hi))
        for (gene,), g in d.group_by("gene"):
            y, w = g["y"].to_numpy().astype(bool), g["w"].to_numpy().astype(float)
            if y.all() or not y.any():
                continue
            for p in preds:
                num, den = pair_stats(y, w, g[p].to_numpy().astype(float))
                rows.append((gene, window, p, num, den))
    return pl.DataFrame(rows, schema=["gene", "window", "pred", "num", "den"], orient="row")


@dataclass
class DiD:
    estimate: float
    se: float
    lo: float
    hi: float
    treated_pre: float
    treated_post: float
    control_pre: float
    control_post: float
    replicates: np.ndarray


class StratifiedAUC:
    """Gene-by-key matrices of (num, den) with point estimates and cluster-bootstrap inference."""

    def __init__(self, stats: pl.DataFrame):
        self.genes = np.array(sorted(stats["gene"].unique().to_list()))
        keys = stats.select("window", "pred").unique().sort("window", "pred").iter_rows()
        self.keys = [tuple(k) for k in keys]
        g_idx = {g: i for i, g in enumerate(self.genes)}
        k_idx = {k: i for i, k in enumerate(self.keys)}
        self.num = np.zeros((len(self.genes), len(self.keys)))
        self.den = np.zeros_like(self.num)
        for gene, window, pred, num, den in stats.iter_rows():
            self.num[g_idx[gene], k_idx[(window, pred)]] = num
            self.den[g_idx[gene], k_idx[(window, pred)]] = den
        self._k = k_idx

    def auc(self, window: str, pred: str, counts: np.ndarray | None = None):
        """WG-AUC, or one value per row of `counts` (gene multiplicities of bootstrap replicates)."""
        k = self._k[(window, pred)]
        if counts is None:
            return self.num[:, k].sum() / self.den[:, k].sum()
        return (counts @ self.num[:, k]) / (counts @ self.den[:, k])

    def n_pairs(self, window: str, pred: str) -> float:
        return float(self.den[:, self._k[(window, pred)]].sum())

    def n_genes(self, window: str) -> int:
        k = next(i for key, i in self._k.items() if key[0] == window)
        return int((self.den[:, k] > 0).sum())

    def bootstrap_counts(self, replicates: int, seed: int) -> np.ndarray:
        """Multinomial gene multiplicities: replicates x genes, rows summing to the number of genes."""
        g = len(self.genes)
        return np.random.default_rng(seed).multinomial(g, np.full(g, 1.0 / g), size=replicates)

    def did(
        self, treated: str, control: str, pre: str, post: str,
        replicates: int = 2000, seed: int = 20261012, counts: np.ndarray | None = None,
    ) -> DiD:
        """[WG-AUC(T, post) - WG-AUC(T, pre)] - [WG-AUC(C, post) - WG-AUC(C, pre)] with a gene-cluster bootstrap."""
        if counts is None:
            counts = self.bootstrap_counts(replicates, seed)

        def stat(c=None):
            return (self.auc(post, treated, c) - self.auc(pre, treated, c)) - (
                self.auc(post, control, c) - self.auc(pre, control, c)
            )

        reps = stat(counts)
        lo, hi = np.percentile(reps, [2.5, 97.5])
        return DiD(
            estimate=float(stat()), se=float(reps.std(ddof=1)), lo=float(lo), hi=float(hi),
            treated_pre=float(self.auc(pre, treated)), treated_post=float(self.auc(post, treated)),
            control_pre=float(self.auc(pre, control)), control_post=float(self.auc(post, control)),
            replicates=reps,
        )


def minimum_detectable_effect(se: float) -> float:
    """Effect detectable with 80% power at two-sided alpha 0.05, under a normal approximation."""
    return 2.8 * se
