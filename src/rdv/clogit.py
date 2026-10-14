"""Weighted conditional logistic regression for matched risk sets, with gene-clustered errors.

Each stratum holds one case and up to M-1 controls. Position 0 is the case; `valid` marks real
members (padding is zero). With weights w_s the log-likelihood is

    sum_s w_s [ x_s0' b - log sum_{j valid} exp(x_sj' b) ].

The covariance is the cluster-robust sandwich A^-1 B A^-1 with A the negative Hessian and B the sum
over clusters of the outer products of the cluster-summed scores, scaled by G / (G - 1).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Fit:
    beta: np.ndarray
    cov: np.ndarray  # cluster-robust
    model_cov: np.ndarray  # inverse information, assuming independent strata
    loglik: float
    iterations: int
    converged: bool
    n_strata: int
    n_clusters: int

    @property
    def se(self) -> np.ndarray:
        return np.sqrt(np.diag(self.cov))


def _pieces(beta, X, valid, w):
    eta = np.where(valid, X @ beta, -np.inf)
    m = eta.max(axis=1, keepdims=True)
    e = np.exp(eta - m)
    z = e.sum(axis=1, keepdims=True)
    pi = e / z
    xbar = np.einsum("sm,smp->sp", pi, X)
    loglik = float(np.sum(w * (eta[:, 0] - m[:, 0] - np.log(z[:, 0]))))
    return pi, xbar, loglik


def loglik_grad_hess(beta, X, valid, w):
    pi, xbar, ll = _pieces(beta, X, valid, w)
    grad = np.einsum("s,sp->p", w, X[:, 0, :] - xbar)
    second = np.einsum("s,sm,smp,smq->pq", w, pi, X, X) - np.einsum("s,sp,sq->pq", w, xbar, xbar)
    return ll, grad, -second  # Hessian of the log-likelihood


def fit(X, valid, w, cluster, max_iter: int = 60, tol: float = 1e-10, ridge: float = 1e-9) -> Fit:
    """Newton-Raphson fit. `cluster` holds one integer id per stratum (the gene)."""
    X = np.asarray(X, dtype=float)
    valid = np.asarray(valid, dtype=bool)
    w = np.asarray(w, dtype=float)
    S, _, P = X.shape
    beta = np.zeros(P)
    ll_old = -np.inf
    converged = False
    for it in range(1, max_iter + 1):
        ll, grad, hess = loglik_grad_hess(beta, X, valid, w)
        step = np.linalg.solve(-hess + ridge * np.eye(P), grad)
        # backtracking keeps the likelihood from decreasing on badly scaled problems
        t = 1.0
        while t > 1e-6:
            cand = beta + t * step
            if _pieces(cand, X, valid, w)[2] >= ll - 1e-12:
                break
            t /= 2
        beta = beta + t * step
        if abs(ll - ll_old) < tol * (1 + abs(ll)) and np.max(np.abs(t * step)) < 1e-8:
            converged = True
            break
        ll_old = ll
    ll, _, hess = loglik_grad_hess(beta, X, valid, w)
    pi, xbar, _ = _pieces(beta, X, valid, w)
    A = -hess + ridge * np.eye(P)
    A_inv = np.linalg.inv(A)
    scores = w[:, None] * (X[:, 0, :] - xbar)
    _, inv = np.unique(np.asarray(cluster), return_inverse=True)
    G = inv.max() + 1
    cs = np.zeros((G, P))
    np.add.at(cs, inv, scores)
    B = cs.T @ cs
    cov = (G / max(G - 1, 1)) * A_inv @ B @ A_inv
    return Fit(beta, cov, A_inv, ll, it, converged, S, int(G))


def information_at_null(cov_by_group: np.ndarray, case_weight: np.ndarray) -> np.ndarray:
    """Expected information at beta = 0 for risk sets whose members are exchangeable.

    cov_by_group: (K, P, P) covariance of the covariates within each risk set's pool.
    case_weight: (K,) total weight of the cases that fall in that pool.
    """
    return np.einsum("k,kpq->pq", case_weight, cov_by_group)
