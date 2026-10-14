import numpy as np
import polars as pl
import pytest

from rdv.metrics import StratifiedAUC, minimum_detectable_effect, pair_stats, per_gene_stats


def brute(y, w, s):
    num = den = 0.0
    for i in np.flatnonzero(y):
        for j in np.flatnonzero(~y):
            den += w[i] * w[j]
            num += w[i] * w[j] * ((s[i] > s[j]) + 0.5 * (s[i] == s[j]))
    return num, den


def test_pair_stats_matches_brute_force_with_ties_and_weights():
    rng = np.random.default_rng(0)
    for _ in range(25):
        n = int(rng.integers(2, 40))
        y = rng.random(n) < 0.4
        w = rng.choice([0.2, 0.5, 1.0, 3.0], n)
        s = rng.integers(0, 6, n).astype(float)  # many ties
        got, want = pair_stats(y, w, s), brute(y, w, s)
        assert got == pytest.approx(want)


def test_pair_stats_extremes_and_single_class():
    y = np.array([True, True, False, False])
    w = np.ones(4)
    assert pair_stats(y, w, np.array([4.0, 3.0, 2.0, 1.0])) == (4.0, 4.0)  # perfect
    assert pair_stats(y, w, np.array([1.0, 2.0, 3.0, 4.0])) == (0.0, 4.0)  # inverted
    assert pair_stats(y, w, np.zeros(4)) == (2.0, 4.0)  # all tied: AUC 0.5
    assert pair_stats(np.array([True, True]), np.ones(2), np.array([1.0, 2.0])) == (0.0, 0.0)


def test_integer_weights_equal_replication():
    y = np.array([True, False, True, False])
    s = np.array([0.9, 0.5, 0.2, 0.4])
    w = np.array([2.0, 1.0, 3.0, 2.0])
    rep = np.repeat(np.arange(4), w.astype(int))
    assert pair_stats(y, w, s) == pytest.approx(pair_stats(y[rep], np.ones(len(rep)), s[rep]))


def _frame(rows):
    return pl.DataFrame(rows, schema=["gene", "y", "w", "snapshot", "a", "b"], orient="row")


def test_per_gene_stats_pools_genes_and_ignores_single_class_genes():
    df = _frame([
        ("G1", 1, 1.0, "2022-01", 0.9, 0.1), ("G1", 0, 1.0, "2022-01", 0.1, 0.9),
        ("G2", 1, 1.0, "2022-01", 0.8, 0.8), ("G2", 0, 1.0, "2022-01", 0.9, 0.8),  # a inverted, b tied
        ("G3", 1, 1.0, "2022-01", 0.5, 0.5), ("G3", 1, 1.0, "2022-01", 0.6, 0.6),  # only P: no pairs
        ("G1", 1, 1.0, "2023-01", 0.9, 0.9), ("G1", 0, 1.0, "2023-01", 0.1, 0.1),  # outside the window
    ])
    stats = per_gene_stats(df, ["a", "b"], {"pre": ("2021-01", "2022-10")})
    assert set(stats["gene"]) == {"G1", "G2"}
    m = StratifiedAUC(stats)
    assert m.auc("pre", "a") == pytest.approx((1 + 0) / 2)  # G1 right, G2 wrong
    assert m.auc("pre", "b") == pytest.approx((0 + 0.5) / 2)  # G1 wrong, G2 tied
    assert (m.n_pairs("pre", "a"), m.n_genes("pre")) == (2.0, 2)


def _simulate(effect_post, n_genes=300, per_gene=20, seed=3):
    """Treated predictor T separates classes by 1.0 before and 1.0+effect_post after; control C stays at 1.0."""
    rng = np.random.default_rng(seed)
    rows = []
    for g in range(n_genes):
        for snap, d_t in (("2021-04", 1.0), ("2023-04", 1.0 + effect_post)):
            y = rng.random(per_gene) < 0.5
            t = rng.normal(size=per_gene) + d_t * y
            c = rng.normal(size=per_gene) + 1.0 * y
            rows += [(f"G{g}", int(yy), 1.0, snap, tt, cc) for yy, tt, cc in zip(y, t, c)]
    return _frame(rows)


WINDOWS = {"pre": ("2021-01", "2022-10"), "post": ("2023-01", "2024-10")}


def test_did_recovers_a_planted_gain_and_gives_a_null_interval_for_no_gain():
    gain = StratifiedAUC(per_gene_stats(_simulate(0.5), ["a", "b"], WINDOWS))  # a = T, b = C
    d = gain.did("a", "b", "pre", "post", replicates=500, seed=1)
    # AUC = Phi(d / sqrt 2): T goes 0.760 -> 0.856 while C stays at 0.760, so the planted DiD is 0.095
    assert d.estimate == pytest.approx(0.095, abs=0.02) and d.lo > 0
    null = StratifiedAUC(per_gene_stats(_simulate(0.0, seed=4), ["a", "b"], WINDOWS))
    d0 = null.did("a", "b", "pre", "post", replicates=500, seed=1)
    assert d0.lo < 0 < d0.hi and abs(d0.estimate) < 0.03
    assert d.replicates.shape == (500,) and d.se > 0 and minimum_detectable_effect(d.se) == pytest.approx(2.8 * d.se)


def test_bootstrap_counts_are_gene_multiplicities():
    m = StratifiedAUC(per_gene_stats(_simulate(0.0, n_genes=50), ["a"], WINDOWS))
    c = m.bootstrap_counts(20, seed=0)
    assert c.shape == (20, 50) and (c.sum(axis=1) == 50).all()
    # resampling every gene exactly once reproduces the point estimate
    assert m.auc("pre", "a", np.ones((1, 50))) == pytest.approx(m.auc("pre", "a"))
