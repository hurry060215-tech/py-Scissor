"""Gaussian and binomial network paths.

Official algorithms translated/adapted from Duanchen Sun's Scissor, fixed SHA
311560a1f160f665917e7da0fcf289914ee99773, src/Scissor.cpp (GPL; see NOTICE).
2026 modifications: sparse CSR access, input guards, explicit semantics,
diagnostics, and an independently converged corrected Gaussian path.
"""
import numpy as np
from scipy import sparse
from scipy.special import expit
from ._types import ScissorError


def soft(value, threshold):
    return np.sign(value) * max(abs(value) - threshold, 0.0)


def neighbor(W, beta, j):
    a, b = W.indptr[j:j+2]
    return float(W.data[a:b] @ beta[W.indices[a:b]])


def network_quadratic(beta, W):
    return float(beta @ beta - beta @ (W @ beta))


def kkt_residual(gradient, beta, l1):
    active = beta != 0
    residual = np.maximum(np.abs(gradient) - l1, 0)
    residual[active] = np.abs(gradient[active] + l1 * np.sign(beta[active]))
    return float(np.max(residual, initial=0))


def wrapper_network(W, family, semantics, alpha):
    """Represent the published wrapper, refusing its undefined memory reads."""
    p = W.shape[0]
    if semantics == "corrected" or family == "binomial" or alpha == 1:
        return W
    if p == 1 and family == "gaussian":
        return W
    unsafe = 2 if family == "gaussian" else 1
    if W[max(0, p-unsafe):].nnz:
        raise ScissorError(
            f"official {family} wrapper prepends an intercept graph node but its kernel has p coefficients; "
            f"this graph can trigger out-of-bounds access (last {unsafe} graph nodes must be isolated for the audited safe subset). "
            "Use semantics='corrected' explicitly for aligned graph behavior; it is not strict official equivalence.")
    return sparse.block_diag((sparse.csr_matrix((1, 1)), W), format="csr")[:p, :p].tocsr()


def gaussian_path(X, y, W, lambdas, alpha, *, semantics, thresh, maxit, cv=False):
    n, p = X.shape
    meanX, meany = X.mean(axis=0), y.mean()
    Z, residual = X-meanX, y-meany
    beta = np.zeros(p)
    output, intercepts, flags, its, errors, rss_track = [], [], [], [], [], []
    initial_rss = float(residual @ residual)
    rss = initial_rss
    objective_old = rss/n/2
    quadratic = 0.
    active, active_mask = [], np.zeros(p, dtype=bool)
    correlations = np.abs(Z.T @ residual/n)
    curvature = np.mean(Z*Z, axis=0)
    early_truncated = False
    for lam in lambdas:
        l1, l2 = lam*alpha, lam*(1-alpha)
        if semantics == "corrected":
            converged = False
            gradient = -Z.T @ residual/n + l2*(beta-W@beta)
            for iteration in range(1, maxit+1):
                # Zero coordinates satisfying their current KKT inequality
                # need no update. Recheck *all* coordinates after every sweep,
                # including nodes newly affected by their graph neighbors.
                # This is an active-set optimization, not a screening rule:
                # no cell is removed, and the same full KKT stop is required.
                working = np.flatnonzero((beta != 0) | (np.abs(gradient) > l1))
                for j in working:
                    denominator = curvature[j]+l2
                    new = (soft(Z[:, j] @ residual/n + curvature[j]*beta[j] + l2*neighbor(W, beta, j), l1)/denominator
                           if denominator > 0 else 0.)
                    residual += (beta[j]-new)*Z[:, j]
                    beta[j] = new
                gradient = -Z.T @ residual/n + l2*(beta-W@beta)
                error = kkt_residual(gradient, beta, l1)
                if not np.isfinite(beta).all() or not np.isfinite(error):
                    break
                if error <= thresh:
                    converged = True
                    break
            rss = float(residual@residual)
        else:
            for j in range(p):
                if not active_mask[j] and correlations[j] > l1:
                    active.append(j); active_mask[j] = True
            iteration = 0
            converged = False
            while True:
                converged = False
                while True:
                    iteration += 1
                    old_rss, old_obj = rss, objective_old
                    absolute, max_change = 0., 0.
                    for j in active:
                        xr = float(residual@Z[:, j])
                        adj = neighbor(W, beta, j)
                        new = soft(xr/n+beta[j]+l2*adj, l1)/(1+l2)
                        change = beta[j]-new
                        quadratic -= change*(beta[j]+new-2*adj)
                        beta[j] = new
                        residual += change*Z[:, j]
                        absolute += abs(new)*l1
                        rss += change*(change*n+2*xr)
                        max_change = max(max_change, change*change)
                    objective_old = rss/n/2+absolute+l2*quadratic/2
                    if not np.isfinite(objective_old) or not np.isfinite(beta).all():
                        break
                    if (abs(old_rss-rss) < abs(1e-5*old_rss) or rss < initial_rss*.001
                            or max_change < thresh or abs(old_obj-objective_old) < abs(1e-5*old_obj)):
                        converged = True
                        break
                    if iteration >= maxit:
                        break
                added = False
                for j in range(p):
                    if not active_mask[j]:
                        # Published NetLmC uses loc without -1 here. Preserve
                        # only within the bounds-checked safe graph subset.
                        adj = 0.
                        start, end = W.indptr[j:j+2]
                        for m0 in W.indices[start:end]:
                            m = m0+1
                            if m < p and active_mask[m]:
                                adj += beta[m]*W[m, j]
                        correlations[j] = abs(residual@Z[:, j]/n+l2*adj)
                        if correlations[j] > l1:
                            active.append(j); active_mask[j] = True; added = True
                if not added or iteration >= maxit or not np.isfinite(objective_old):
                    if added:
                        converged = False
                    break
            gradient = -Z.T@residual/n+l2*(beta-W@beta)
            error = kkt_residual(gradient, beta, l1)
            if initial_rss > 0 and 1-rss/initial_rss > .999 and not cv:
                # C++ returns il, thereby excluding the just-computed column.
                early_truncated = True
                break
        output.append(beta.copy()); intercepts.append(float(meany-meanX@beta))
        flags.append(converged); its.append(iteration); errors.append(error); rss_track.append(rss)
    return (np.asarray(output).T if output else np.empty((p, 0)), np.asarray(intercepts),
            np.asarray(flags, bool), np.asarray(its, int),
            {"kkt_residual": errors, "reported_rss": rss_track,
             "early_truncated": early_truncated,
             "convergence_definition": "KKT residual" if semantics == "corrected" else "published stopping tests; inspect KKT diagnostics"})


def binomial_path(X, y, W, lambdas, alpha, *, semantics, thresh, maxit):
    if semantics == "corrected":
        return _corrected_binomial_path(X, y, W, lambdas, alpha, thresh, maxit)
    n, p0 = X.shape
    means = X.mean(axis=0)
    Z = np.column_stack([np.ones(n), X-means])
    p = p0+1
    W = sparse.block_diag((sparse.csr_matrix((1, 1)), W), format="csr")
    beta = np.zeros(p); beta[0] = np.log(y.mean()/(1-y.mean()))
    eta = np.full(n, beta[0])
    prob = np.clip(expit(eta), 1e-5, 1-1e-5)
    correlations = np.abs(Z.T@(y-prob)/n)
    active, mask = [0], np.zeros(p, bool); mask[0] = True
    output, intercepts, flags, its, errors = [], [], [], [], []
    square = Z*Z
    for lam in lambdas:
        l1, l2 = np.full(p, lam*alpha), np.full(p, lam*(1-alpha))
        l1[0] = l2[0] = 0
        if semantics == "corrected":
            active, mask = list(range(p)), np.ones(p, bool)
        else:
            for j in range(1, p):
                if not mask[j] and correlations[j] > l1[j]:
                    active.append(j); mask[j] = True
        iteration, converged = 0, False
        while True:
            converged = False
            while True:
                iteration += 1
                weights = prob*(1-prob)
                working = (y-prob)/weights
                largest = 0.
                for j in active:
                    h = float(weights@square[:, j]/n)
                    partial = float((working*weights)@Z[:, j]/n+h*beta[j])
                    partial += l2[j]*neighbor(W, beta, j)
                    new = soft(partial, l1[j])/(h+l2[j]) if h+l2[j] > 0 else 0.
                    change = new-beta[j]
                    beta[j] = new
                    working -= change*Z[:, j]
                    eta += change*Z[:, j]
                    largest = max(largest, change*change)
                prob = np.clip(expit(eta), 1e-5, 1-1e-5)
                gradient = Z.T@(prob-y)/n+l2*(beta-W@beta)
                # Intercept has l1=l2=0; include its score in convergence.
                active_beta = beta != 0
                r = np.maximum(np.abs(gradient)-l1, 0)
                r[active_beta] = np.abs(gradient[active_beta]+l1[active_beta]*np.sign(beta[active_beta]))
                error = float(np.max(r))
                if not np.isfinite(beta).all() or not np.isfinite(error):
                    break
                if (largest < thresh if semantics == "official" else error <= thresh):
                    converged = True
                    break
                if iteration >= maxit:
                    break
            if semantics == "corrected":
                break
            added = False
            for j in range(1, p):
                if not mask[j]:
                    # Intentional compatibility with the upstream signed term
                    # outside abs(); corrected mode checks every coordinate.
                    correlations[j] = abs(Z[:, j]@(y-prob)/n)+l2[j]*neighbor(W, beta, j)
                    if correlations[j] > l1[j]:
                        active.append(j); mask[j] = True; added = True
            if not added or iteration >= maxit or not np.isfinite(error):
                if added:
                    converged = False
                break
        output.append(beta[1:].copy()); intercepts.append(float(beta[0]-means@beta[1:]))
        flags.append(converged); its.append(iteration); errors.append(error)
    return (np.asarray(output).T, np.asarray(intercepts), np.asarray(flags, bool),
            np.asarray(its, int), {"kkt_residual": errors, "probability_clip": 1e-5,
            "convergence_definition": "KKT residual of clipped score" if semantics == "corrected" else "published squared coefficient change"})


def _corrected_binomial_path(X, y, W, lambdas, alpha, thresh, maxit):
    """Exact logistic objective with monotone accelerated proximal gradient.

    Unlike official mode, probabilities are not clipped for fitting. A stable
    logaddexp objective, line search, and a KKT check determine convergence.
    Published clipped deviance is still used by cross-validation in both modes.
    """
    n, p = X.shape
    means = X.mean(axis=0)
    Z = np.column_stack([np.ones(n), X-means])
    beta = np.r_[np.log(y.mean()/(1-y.mean())), np.zeros(p)]
    output, intercepts, flags, its, errors, objectives, reasons = [], [], [], [], [], [], []
    lipschitz = max(np.max(np.sum(Z*Z, axis=0))/n*.25, 1e-6)
    for lam in lambdas:
        l1, l2 = lam*alpha, lam*(1-alpha)
        def smooth(b):
            eta = Z@b
            lap = b[1:]-W@b[1:]
            value = np.mean(np.logaddexp(0, (1-2*y)*eta))+l2*(b[1:]@lap)/2
            gradient = Z.T@(expit(eta)-y)/n
            gradient[1:] += l2*lap
            return float(value), gradient
        value, gradient = smooth(beta)
        objective = value+l1*np.abs(beta[1:]).sum()
        extrapolated, momentum = beta.copy(), 1.
        converged, reason = False, "maxit"
        for iteration in range(1, maxit+1):
            error = max(abs(gradient[0]), kkt_residual(gradient[1:], beta[1:], l1))
            if error <= thresh:
                converged, reason = True, "kkt"
                break
            base_value, base_gradient = smooth(extrapolated)
            accepted = False
            for restart in range(2):
                stepL = max(lipschitz*.8, 1e-12)
                for _ in range(80):
                    proposal = extrapolated-base_gradient/stepL
                    proposal[1:] = np.sign(proposal[1:])*np.maximum(abs(proposal[1:])-l1/stepL, 0)
                    delta = proposal-extrapolated
                    next_value, next_gradient = smooth(proposal)
                    bound = base_value+base_gradient@delta+stepL*(delta@delta)/2
                    if np.isfinite(next_value) and next_value <= bound+1e-13*max(1, abs(base_value)):
                        accepted = True
                        break
                    stepL *= 2
                if not accepted:
                    break
                next_objective = next_value+l1*np.abs(proposal[1:]).sum()
                if next_objective <= objective+1e-13*max(1, abs(objective)):
                    break
                accepted = False
                if restart == 0:
                    extrapolated = beta.copy(); momentum = 1.
                    base_value, base_gradient = value, gradient
            if not accepted:
                reason = "line_search_failed"
                break
            old = beta
            beta = proposal
            value, gradient, objective = next_value, next_gradient, next_objective
            lipschitz = stepL
            next_momentum = (1+np.sqrt(1+4*momentum*momentum))/2
            extrapolated = beta+(momentum-1)/next_momentum*(beta-old)
            momentum = next_momentum
        error = max(abs(gradient[0]), kkt_residual(gradient[1:], beta[1:], l1))
        if error <= thresh:
            converged, reason = True, "kkt"
        output.append(beta[1:].copy()); intercepts.append(float(beta[0]-means@beta[1:]))
        flags.append(converged); its.append(iteration); errors.append(error)
        objectives.append(float(objective)); reasons.append(reason)
    return (np.asarray(output).T, np.asarray(intercepts), np.asarray(flags, bool), np.asarray(its, int),
            {"kkt_residual": errors, "objective": objectives, "termination_reason": reasons,
             "solver": "monotone_accelerated_proximal_gradient", "probability_clip": None,
             "convergence_definition": "absolute KKT residual of exact logistic objective"})
