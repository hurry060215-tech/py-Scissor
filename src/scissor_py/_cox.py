"""Independent Cox/Breslow numerical kernels for the Scissor prototype.

Provenance: mathematical behavior was checked against Scissor's ``PrepCox``,
``NetCoxC`` and ``EnetCoxC`` at upstream commit
311560a1f160f665917e7da0fcf289914ee99773 (R/net_cox.R and src/Scissor.cpp).
This Python implementation is part of the GPL-3.0-or-later prototype; see NOTICE.
It contains no R bridge and never constructs a subject-by-subject risk matrix.

The corrected solver minimizes
    -log_partial_likelihood / n_original
    + lambda * alpha * ||beta||_1
    + lambda * (1-alpha) / 2 * beta.T @ (I-W) @ beta.
``official`` reproduces the defined upstream coordinate-update/stopping behavior,
including early path termination. Graph preparation and lambda-label semantics
belong to the caller. It is not a bit-for-bit floating-point compatibility claim.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse


@dataclass(frozen=True)
class _CoxData:
    x: np.ndarray
    status: np.ndarray
    event_starts: np.ndarray
    event_counts: np.ndarray
    event_rows: np.ndarray
    event_risk_starts: np.ndarray
    n_original: int
    mean: np.ndarray


def _prepare(X, y) -> _CoxData:
    """Sort tied events first and center on records at/after the first event."""
    x = np.asarray(X, dtype=float)
    outcome = np.asarray(y, dtype=float)
    if x.ndim != 2 or min(x.shape) == 0:
        raise ValueError("X must be a nonempty two-dimensional matrix")
    if outcome.shape != (x.shape[0], 2):
        raise ValueError("Cox y must have shape (n_samples, 2): time, status")
    if not np.isfinite(x).all() or not np.isfinite(outcome).all():
        raise ValueError("X and Cox outcomes must contain only finite values")
    if np.any(outcome[:, 0] < 0):
        raise ValueError("Cox times must be nonnegative")
    if not np.isin(outcome[:, 1], (0, 1)).all():
        raise ValueError("Cox status must be 0 or 1")
    n_original = x.shape[0]
    order = np.lexsort((-outcome[:, 1], outcome[:, 0]))
    x = x[order].copy()
    time, status = outcome[order].T
    events = np.flatnonzero(status)
    # With no events, the partial likelihood is identically zero. Path fitting
    # rejects this case, but public likelihood/score evaluation remains defined.
    if events.size:
        offset = events[0]
        x, time, status = x[offset:], time[offset:], status[offset:]
    mean = x.mean(axis=0)
    x -= mean
    starts = np.r_[0, np.flatnonzero(np.diff(time) != 0) + 1]
    counts = np.add.reduceat(status, starts).astype(int)
    event_starts = starts[counts > 0]
    event_counts = counts[counts > 0]
    return _CoxData(
        x, status, event_starts, event_counts, np.flatnonzero(status),
        np.repeat(event_starts, event_counts), n_original, mean,
    )


def _eta_stats(data: _CoxData, eta: np.ndarray) -> tuple[float, np.ndarray]:
    """Breslow likelihood and eta score in O(n) time and temporary memory.

    Suffix log-sums avoid exp overflow and subtraction of nearly equal risk
    sums. Cumulative hazard is also accumulated in log space, so no product
    of an overflowing reciprocal risk sum and an underflowing exp is needed.
    """
    if not data.event_starts.size:
        return 0.0, np.zeros_like(eta)
    log_risk = np.logaddexp.accumulate(eta[::-1])[::-1]
    ll = float(np.sum(eta[data.event_rows] - log_risk[data.event_risk_starts]))
    log_increment = np.full(eta.size, -np.inf)
    log_increment[data.event_starts] = (
        np.log(data.event_counts) - log_risk[data.event_starts]
    )
    log_hazard = np.logaddexp.accumulate(log_increment)
    expected = np.exp(eta + log_hazard)
    return ll, data.status - expected


def _beta(beta, p: int) -> np.ndarray:
    value = np.asarray(beta, dtype=float)
    if value.shape != (p,) or not np.isfinite(value).all():
        raise ValueError("beta must be a finite vector with one value per feature")
    return value


def cox_log_likelihood(X, y, beta) -> float:
    """Unnormalized Breslow log partial likelihood, using all risk-set records.

    Censoring at an event time remains in that event's risk set. Records
    censored strictly before the first event contribute zero. No intercept
    or baseline hazard is estimated.
    """
    data = _prepare(X, y)
    value = _beta(beta, data.x.shape[1])
    return _eta_stats(data, data.x @ value)[0]


def cox_score(X, y, beta) -> np.ndarray:
    """Gradient of ``cox_log_likelihood`` with respect to beta (not / n)."""
    data = _prepare(X, y)
    value = _beta(beta, data.x.shape[1])
    return data.x.T @ _eta_stats(data, data.x @ value)[1]


def cox_lambda_max(X, y, alpha: float) -> float:
    """Null-score maximum / original n / alpha; use 0.001 when alpha is zero."""
    if not np.isfinite(alpha) or not 0 <= alpha <= 1:
        raise ValueError("alpha must lie in [0, 1]")
    data = _prepare(X, y)
    if not data.event_starts.size:
        raise ValueError("At least one observed event is required to fit a Cox model")
    score = data.x.T @ _eta_stats(data, np.zeros(data.x.shape[0]))[1]
    return float(np.max(np.abs(score)) / data.n_original / (alpha or 0.001))


def _network(network, p: int) -> sparse.csr_matrix:
    if network is None:
        return sparse.csr_matrix((p, p), dtype=float)
    if not sparse.issparse(network):
        raise ValueError("network must be a pre-normalized sparse matrix or None")
    result = sparse.csr_matrix(network, dtype=float, copy=True)
    result.sum_duplicates()
    result.eliminate_zeros()
    result.sort_indices()
    if result.shape != (p, p):
        raise ValueError("network must have shape (n_features, n_features)")
    if not np.isfinite(result.data).all():
        raise ValueError("network weights must be finite")
    if np.any(result.diagonal() != 0):
        raise ValueError("network must have zero diagonal")
    asymmetry = result - result.T
    if asymmetry.nnz and np.max(np.abs(asymmetry.data)) > 1e-12:
        raise ValueError("network must be symmetric")
    return result


def _soft_threshold(value: np.ndarray, amount: float) -> np.ndarray:
    return np.sign(value) * np.maximum(np.abs(value) - amount, 0.0)


def _kkt_residual(gradient: np.ndarray, beta: np.ndarray, l1: float) -> float:
    residual = np.maximum(np.abs(gradient) - l1, 0.0)
    nonzero = beta != 0
    residual[nonzero] = np.abs(gradient[nonzero] + l1 * np.sign(beta[nonzero]))
    return float(np.max(residual))


def _corrected_path(data, graph, lambdas, alpha, thresh, maxit):
    x, n, p = data.x, data.n_original, data.x.shape[1]
    beta = np.zeros(p)
    coefficients = np.zeros((p, len(lambdas)))
    converged = np.zeros(len(lambdas), dtype=bool)
    iterations = np.zeros(len(lambdas), dtype=int)
    residuals = np.full(len(lambdas), np.inf)
    objectives = np.full(len(lambdas), np.nan)
    reasons = []
    # Reused between lambdas, but backtracking always verifies local curvature.
    lipschitz = max(float(np.max(np.sum(x * x, axis=0))) / n, 1e-6)

    for column, lam in enumerate(lambdas):
        l1, l2 = float(lam * alpha), float(lam * (1 - alpha))

        def smooth(value):
            eta = x @ value
            if not np.isfinite(eta).all():
                return np.inf, np.full(p, np.nan)
            ll, score_eta = _eta_stats(data, eta)
            lap_value = value - graph @ value
            loss = -ll / n + 0.5 * l2 * float(value @ lap_value)
            gradient = -(x.T @ score_eta) / n + l2 * lap_value
            return loss, gradient

        smooth_beta, grad_beta = smooth(beta)
        objective = smooth_beta + l1 * float(np.abs(beta).sum())
        extrapolated = beta.copy()
        momentum = 1.0
        reason = "maxit"
        for iteration in range(1, maxit + 1):
            residual = _kkt_residual(grad_beta, beta, l1)
            if residual <= thresh:
                converged[column] = True
                reason = "kkt"
                # A warm start already satisfying KKT performs no update.
                iterations[column] = iteration - 1
                break
            base_loss, base_gradient = smooth(extrapolated)
            accepted = False
            # Monotone accelerated proximal gradient. If acceleration raises
            # the objective, restart at the last accepted coefficient vector.
            for restart in range(2):
                trial_lipschitz = max(lipschitz * 0.8, 1e-12)
                for _ in range(80):
                    proposal = _soft_threshold(
                        extrapolated - base_gradient / trial_lipschitz,
                        l1 / trial_lipschitz,
                    )
                    delta = proposal - extrapolated
                    proposal_loss, proposal_gradient = smooth(proposal)
                    bound = (base_loss + float(base_gradient @ delta)
                             + 0.5 * trial_lipschitz * float(delta @ delta))
                    slack = 1e-13 * max(1.0, abs(base_loss))
                    if np.isfinite(proposal_loss) and proposal_loss <= bound + slack:
                        accepted = True
                        break
                    trial_lipschitz *= 2
                if not accepted:
                    break
                proposal_objective = proposal_loss + l1 * float(np.abs(proposal).sum())
                if proposal_objective <= objective + 1e-13 * max(1.0, abs(objective)):
                    break
                if restart == 0:
                    extrapolated = beta.copy()
                    base_loss, base_gradient = smooth_beta, grad_beta
                    momentum = 1.0
                    accepted = False
                else:
                    accepted = False
            if not accepted:
                reason = "line_search_failed"
                iterations[column] = iteration - 1
                break
            previous = beta
            beta = proposal
            smooth_beta, grad_beta = proposal_loss, proposal_gradient
            objective = proposal_objective
            lipschitz = trial_lipschitz
            next_momentum = (1 + np.sqrt(1 + 4 * momentum * momentum)) / 2
            extrapolated = beta + ((momentum - 1) / next_momentum) * (beta - previous)
            momentum = next_momentum
            iterations[column] = iteration
        residuals[column] = _kkt_residual(grad_beta, beta, l1)
        if residuals[column] <= thresh:
            converged[column] = True
            reason = "kkt"
        coefficients[:, column] = beta
        objectives[column] = objective
        reasons.append(reason)
    return coefficients, converged, iterations, {
        "solver": "monotone_accelerated_proximal_gradient",
        "kkt_residual": residuals,
        "objective": objectives,
        "termination_reason": reasons,
        "effective_path_length": len(lambdas),
    }


def _official_stats(data, eta, ifast):
    """Source-style unshifted exponentials and diagonal eta information.

    ``ifast`` uses sequential subtraction, and flags any nonpositive diagonal
    curvature/risk denominator just as the upstream early-termination route.
    The slow route sums surviving records rather than subtracting prior ones.
    """
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        weights = np.exp(eta)
        starts = data.event_starts
        ends = np.r_[starts[1:], len(eta)]
        denominators = np.empty(len(starts))
        if ifast:
            running = float(weights.sum())
            for k, (start, end) in enumerate(zip(starts, ends)):
                denominators[k] = running
                for row in range(start, end):
                    running -= weights[row]
        else:
            running = 0.0
            for k in range(len(starts) - 1, -1, -1):
                for row in range(starts[k], ends[k]):
                    running += weights[row]
                denominators[k] = running
        ll = 0.0
        score = np.empty(len(eta))
        curvature = np.empty(len(eta))
        cumulative, cumulative2 = 0.0, 0.0
        for k, (start, end, count) in enumerate(zip(starts, ends, data.event_counts)):
            denominator = denominators[k]
            # numpy scalar operations intentionally turn invalid arithmetic
            # into flags instead of raising ZeroDivisionError.
            cumulative += count / denominator
            cumulative2 += (count / denominator) / denominator
            segment = weights[start:end]
            score[start:end] = data.status[start:end] - segment * cumulative
            curvature[start:end] = segment * (cumulative - segment * cumulative2)
            for row in range(start, start + count):
                ll += eta[row]
            ll -= count * np.log(denominator)
    invalid = not (np.isfinite(ll) and np.isfinite(score).all()
                   and np.isfinite(curvature).all())
    fast_failure = bool(ifast and (np.any(denominators <= 0)
                                  or np.any(curvature <= 0)))
    return float(ll), score, curvature, invalid, fast_failure


def _official_path(data, graph, lambdas, alpha, thresh, maxit, ifast):
    x, n, p = data.x, data.n_original, data.x.shape[1]
    beta = np.zeros(p)
    eta = np.zeros(x.shape[0])
    active: list[int] = []
    is_active = np.zeros(p, dtype=bool)
    coefficients = np.full((p, len(lambdas)), np.nan)
    converged = np.zeros(len(lambdas), dtype=bool)
    iterations = np.zeros(len(lambdas), dtype=int)
    flags = np.full(len(lambdas), -1, dtype=int)
    residuals = np.full(len(lambdas), np.nan)
    objectives = np.full(len(lambdas), np.nan)
    ll, score, curvature, invalid, fast_failure = _official_stats(data, eta, ifast)
    objective, quadratic = -ll / n, 0.0
    completed = 0
    failure = None

    def admit(l1_scaled, l2_scaled):
        added = False
        for j in range(p):
            if not is_active[j]:
                start, end = graph.indptr[j:j + 2]
                neighbor = float(graph.data[start:end] @ beta[graph.indices[start:end]])
                value = float(score @ x[:, j]) + l2_scaled * neighbor
                if abs(value) > l1_scaled:
                    is_active[j] = True
                    active.append(j)
                    added = True
        return added

    for column, lam in enumerate(lambdas):
        if invalid or fast_failure:
            failure = "nonfinite" if invalid else "upstream_ifast_termination"
            break
        l1, l2 = float(lam * alpha), float(lam * (1 - alpha))
        l1_scaled, l2_scaled = l1 * n, l2 * n
        admit(l1_scaled, l2_scaled)
        stopped = False
        count_iterations = 0
        while not stopped:
            count_iterations += 1
            for j in active:
                col = x[:, j]
                diagonal = float(curvature @ (col * col))
                start, end = graph.indptr[j:j + 2]
                neighbor = float(graph.data[start:end] @ beta[graph.indices[start:end]])
                z = beta[j] * diagonal + float(score @ col) + l2_scaled * neighbor
                denominator = l2_scaled + diagonal
                if abs(z) > l1_scaled:
                    if denominator <= 0 or not np.isfinite(denominator):
                        invalid = True
                        break
                    new = (z - np.sign(z) * l1_scaled) / denominator
                else:
                    new = 0.0
                difference = beta[j] - new
                if difference:
                    quadratic -= difference * (beta[j] + new - 2 * neighbor)
                    beta[j] = new
                    score += curvature * col * difference
                    eta -= difference * col
            if invalid:
                failure = "nonfinite_or_nonpositive_coordinate_curvature"
                break
            old_ll, old_objective = ll, objective
            ll, new_score, new_curvature, invalid, fast_failure = _official_stats(data, eta, ifast)
            objective = -ll / n + l1 * float(np.abs(beta).sum()) + l2 * quadratic / 2
            iterations[column] = count_iterations
            if invalid or fast_failure:
                failure = "nonfinite" if invalid else "upstream_ifast_termination"
                break
            # Upstream tests likelihood OR objective change before maxit.
            stopped = (abs(old_ll - ll) < abs(thresh * old_ll)
                       or abs(old_objective - objective) < abs(thresh * old_objective))
            if not stopped and count_iterations >= maxit:
                failure = "maxit"
                flags[column] = 1
                break
            score, curvature = new_score, new_curvature
            if stopped and admit(l1_scaled, l2_scaled):
                stopped = False
        if failure is not None:
            if flags[column] < 0:
                flags[column] = 2
            # Upstream additionally drops the preceding completed lambda on a
            # fast-risk numerical termination. Preserve that returned prefix.
            if fast_failure and column > 0:
                completed = max(completed - 1, 0)
                coefficients[:, completed:] = np.nan
                converged[completed:] = False
            break
        coefficients[:, column] = beta
        converged[column] = True
        flags[column] = 0
        gradient = -(x.T @ score) / n + l2 * (beta - graph @ beta)
        residuals[column] = _kkt_residual(gradient, beta, l1)
        objectives[column] = objective
        completed += 1
    return coefficients, converged, iterations, {
        "solver": "source_style_coordinate_descent",
        "convergence_definition": "upstream likelihood/objective change, not KKT",
        "flags": flags,
        "kkt_residual": residuals,
        "objective": objectives,
        "effective_path_length": completed,
        "termination_reason": failure or "upstream_tolerance",
    }


def cox_path(X, y, network, lambdas, alpha, *, semantics="corrected",
             thresh=1e-7, maxit=100000, ifast=True):
    """Fit a warm-started Cox path with a pre-normalized sparse graph.

    Returns ``(coefficients, intercepts, converged, iterations, metadata)``.
    Coefficients have shape (p, len(lambdas)); intercepts are always zero.
    ``thresh`` is an absolute KKT tolerance for ``corrected`` and the original
    relative change tolerance for ``official``. The official mode leaves NaN
    columns outside its returned ``effective_path_length``. Corrected mode
    always uses stable suffix risks; ``ifast`` only affects official mode.
    """
    data = _prepare(X, y)
    if not data.event_starts.size:
        raise ValueError("At least one observed event is required to fit a Cox model")
    if semantics not in ("corrected", "official"):
        raise ValueError("semantics must be 'corrected' or 'official'")
    if not np.isfinite(alpha) or not 0 <= alpha <= 1:
        raise ValueError("alpha must lie in [0, 1]")
    penalties = np.asarray(lambdas, dtype=float)
    if (penalties.ndim != 1 or penalties.size == 0
            or not np.isfinite(penalties).all() or np.any(penalties < 0)):
        raise ValueError("lambdas must be a nonempty finite nonnegative vector")
    if not np.isfinite(thresh) or thresh <= 0:
        raise ValueError("thresh must be finite and positive")
    if isinstance(maxit, bool) or int(maxit) != maxit or maxit < 1:
        raise ValueError("maxit must be a positive integer")
    graph = _network(network, data.x.shape[1])
    if semantics == "corrected":
        result = _corrected_path(data, graph, penalties, alpha, float(thresh), int(maxit))
    else:
        result = _official_path(data, graph, penalties, alpha, float(thresh), int(maxit), bool(ifast))
    coefficients, converged, iterations, metadata = result
    metadata.update({
        "semantics": semantics,
        "ties": "breslow",
        "n_original": data.n_original,
        "n_retained": data.x.shape[0],
        "center": data.mean,
        "ifast": bool(ifast),
    })
    return coefficients, np.zeros(penalties.size), converged, iterations, metadata
