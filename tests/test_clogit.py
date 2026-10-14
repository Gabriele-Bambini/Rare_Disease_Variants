import numpy as np
import pytest

from rdv.clogit import fit, information_at_null, loglik_grad_hess


def _simulate(rng, S, M, beta, n_clusters=None, pad_to=None):
    """S strata of M members: covariates ~ N(0, I); the case is drawn with probability ~ exp(x'beta)."""
    P = len(beta)
    X = rng.normal(size=(S, M, P))
    p = np.exp(X @ beta)
    p /= p.sum(axis=1, keepdims=True)
    case = np.array([rng.choice(M, p=pi) for pi in p])
    order = np.array([[case[s]] + [j for j in range(M) if j != case[s]] for s in range(S)])
    X = np.take_along_axis(X, order[:, :, None], axis=1)
    valid = np.ones((S, M), dtype=bool)
    if pad_to:
        X = np.concatenate([X, np.zeros((S, pad_to - M, P))], axis=1)
        valid = np.concatenate([valid, np.zeros((S, pad_to - M), dtype=bool)], axis=1)
    cluster = rng.integers(0, n_clusters or S, size=S)
    return X, valid, cluster


def test_gradient_and_hessian_match_finite_differences_with_padding_and_weights():
    rng = np.random.default_rng(0)
    X, valid, _ = _simulate(rng, 12, 5, np.array([0.5, -0.3]), pad_to=7)
    valid[3, 4:] = False  # a stratum with fewer controls
    X[3, 4:] = 0.0
    w = rng.uniform(0.2, 1.0, 12)
    b = np.array([0.2, 0.1])
    ll, g, H = loglik_grad_hess(b, X, valid, w)
    eps = 1e-6
    for k in range(2):
        d = np.zeros(2)
        d[k] = eps
        up, dn = loglik_grad_hess(b + d, X, valid, w), loglik_grad_hess(b - d, X, valid, w)
        assert (up[0] - dn[0]) / (2 * eps) == pytest.approx(g[k], rel=1e-5, abs=1e-7)
        assert (up[1] - dn[1]) / (2 * eps) == pytest.approx(H[k], rel=1e-5, abs=1e-6)


def test_loglik_equals_brute_force():
    rng = np.random.default_rng(1)
    X, valid, _ = _simulate(rng, 6, 4, np.array([0.4]))
    w = np.array([1.0, 0.5, 2.0, 1.0, 1.0, 0.25])
    b = np.array([0.7])
    brute = sum(w[s] * (X[s, 0] @ b - np.log(sum(np.exp(X[s, j] @ b) for j in range(4)))) for s in range(6))
    assert loglik_grad_hess(b, X, valid, w)[0] == pytest.approx(brute)


def test_fit_recovers_planted_coefficients_and_ignores_padding():
    rng = np.random.default_rng(2)
    beta = np.array([0.8, -0.5, 0.0])
    X, valid, cluster = _simulate(rng, 4000, 8, beta)
    f = fit(X, valid, np.ones(4000), cluster)
    assert f.converged and np.all(np.abs(f.beta - beta) < 4 * f.se)
    pad = np.concatenate([X, np.zeros((4000, 5, 3))], axis=1)
    padv = np.concatenate([valid, np.zeros((4000, 5), dtype=bool)], axis=1)
    f2 = fit(pad, padv, np.ones(4000), cluster)
    assert f2.beta == pytest.approx(f.beta, abs=1e-8)


def test_interaction_recovers_a_planted_change_between_periods():
    rng = np.random.default_rng(3)
    S = 6000
    post = np.repeat([0.0, 1.0], S // 2)
    x = rng.normal(size=(S, 8, 1))
    slope = 0.3 + 0.5 * post  # the coefficient rises by 0.5 after the milestone
    p = np.exp(x[:, :, 0] * slope[:, None])
    p /= p.sum(axis=1, keepdims=True)
    case = np.array([rng.choice(8, p=pi) for pi in p])
    order = np.array([[case[s]] + [j for j in range(8) if j != case[s]] for s in range(S)])
    x = np.take_along_axis(x, order[:, :, None], axis=1)
    X = np.concatenate([x, x * post[:, None, None]], axis=2)
    f = fit(X, np.ones((S, 8), dtype=bool), np.ones(S), rng.integers(0, 300, S))
    assert abs(f.beta[1] - 0.5) < 3 * f.se[1] and f.beta[1] / f.se[1] > 4


def test_weights_act_like_replication():
    rng = np.random.default_rng(4)
    X, valid, cluster = _simulate(rng, 60, 5, np.array([0.6]))
    w = rng.integers(1, 4, 60).astype(float)
    rep = np.repeat(np.arange(60), w.astype(int))
    a = fit(X, valid, w, cluster)
    b = fit(X[rep], valid[rep], np.ones(len(rep)), cluster[rep])
    assert a.beta == pytest.approx(b.beta, abs=1e-8)


def test_clustered_sandwich_covers_about_95_percent():
    rng = np.random.default_rng(5)
    beta, hits, n = np.array([0.5, 0.0]), 0, 300
    for _ in range(n):
        X, valid, cluster = _simulate(rng, 250, 6, beta, n_clusters=40)
        f = fit(X, valid, np.ones(250), cluster)
        hits += abs(f.beta[0] - beta[0]) < 1.96 * f.se[0]
    assert 0.90 <= hits / n <= 0.99


def test_information_at_null_matches_the_hessian_at_zero():
    rng = np.random.default_rng(6)
    X, valid, _ = _simulate(rng, 40, 6, np.array([0.0, 0.0]))
    w = rng.uniform(0.3, 1.0, 40)
    covs = np.stack([np.cov(X[s].T, bias=True) for s in range(40)])  # uniform weights within each stratum
    _, _, H = loglik_grad_hess(np.zeros(2), X, valid, w)
    assert information_at_null(covs, w) == pytest.approx(-H, rel=1e-9)
