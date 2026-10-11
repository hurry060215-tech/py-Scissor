"""Optional float64 acceleration of the corrected Gaussian coordinate loop.

Imported only when engine='numba' is explicitly requested. No fastmath,
parallel JIT, automatic fallback, R bridge, or extra screening semantics.
The reference NumPy implementation remains available and is the default.
"""
import numpy as np
import numba
from numba import njit


@njit(cache=False, fastmath=False, parallel=False)
def _kernel(Z, residual, meanX, meany, curvature, values, indices, indptr,
            lambdas, alpha, thresh, maxit):
    n, p = Z.shape
    size = len(lambdas)
    coefficients = np.empty((p, size), dtype=np.float64)
    intercepts = np.empty(size, dtype=np.float64)
    flags = np.zeros(size, dtype=np.bool_)
    iterations = np.zeros(size, dtype=np.int64)
    errors = np.empty(size, dtype=np.float64)
    rss = np.empty(size, dtype=np.float64)
    beta = np.zeros(p, dtype=np.float64)
    gradient = np.empty(p, dtype=np.float64)
    for k in range(size):
        l1, l2 = lambdas[k]*alpha, lambdas[k]*(1-alpha)
        score = Z.T @ residual
        for j in range(p):
            adjacent = 0.
            for edge in range(indptr[j], indptr[j+1]):
                adjacent += values[edge]*beta[indices[edge]]
            gradient[j] = -score[j]/n+l2*(beta[j]-adjacent)
        error = np.inf
        for iteration in range(1, maxit+1):
            for j in range(p):
                if beta[j] == 0 and abs(gradient[j]) <= l1:
                    continue
                adjacent = 0.
                for edge in range(indptr[j], indptr[j+1]):
                    adjacent += values[edge]*beta[indices[edge]]
                value = np.dot(Z[:, j], residual)/n+curvature[j]*beta[j]+l2*adjacent
                denominator = curvature[j]+l2
                new = 0.
                if denominator > 0:
                    if value > l1:
                        new = (value-l1)/denominator
                    elif value < -l1:
                        new = (value+l1)/denominator
                change = beta[j]-new
                for i in range(n):
                    residual[i] += change*Z[i, j]
                beta[j] = new
            score = Z.T @ residual
            error = 0.
            finite = True
            for j in range(p):
                adjacent = 0.
                for edge in range(indptr[j], indptr[j+1]):
                    adjacent += values[edge]*beta[indices[edge]]
                gradient[j] = -score[j]/n+l2*(beta[j]-adjacent)
                if not np.isfinite(beta[j]) or not np.isfinite(gradient[j]):
                    finite = False
                if beta[j] > 0:
                    deviation = abs(gradient[j]+l1)
                elif beta[j] < 0:
                    deviation = abs(gradient[j]-l1)
                else:
                    deviation = max(abs(gradient[j])-l1, 0.)
                error = max(error, deviation)
            if not finite:
                error = np.inf
                break
            if error <= thresh:
                flags[k] = True
                break
        coefficients[:, k] = beta
        intercepts[k] = meany-np.dot(meanX, beta)
        iterations[k] = iteration
        errors[k] = error
        rss[k] = np.dot(residual, residual)
    return coefficients, intercepts, flags, iterations, errors, rss


def gaussian_numba_path(X, y, W, lambdas, alpha, *, thresh, maxit):
    meanX, meany = X.mean(axis=0), y.mean()
    Z = np.asfortranarray(X-meanX)
    residual = y-meany
    curvature = np.mean(Z*Z, axis=0)
    coefficients, intercepts, flags, iterations, errors, rss = _kernel(
        Z, residual, meanX, meany, curvature, W.data, W.indices, W.indptr,
        lambdas, alpha, thresh, maxit)
    return coefficients, intercepts, flags, iterations, {
        "kkt_residual": errors.tolist(), "reported_rss": rss.tolist(),
        "early_truncated": False, "convergence_definition": "KKT residual",
        "engine": "numba", "numba_version": numba.__version__,
        "jit_fastmath": False, "jit_parallel": False, "jit_disk_cache": False,
    }
