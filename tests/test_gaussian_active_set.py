"""Corrected active-set acceleration retains the full mathematical problem."""
import numpy as np
import pytest
from scipy import sparse

from scissor_py import fit_path, normalize_graph
from scissor_py._solvers import kkt_residual, neighbor, soft


def exhaustive(X, y, W, lambdas, alpha):
    """Previous exhaustive corrected coordinate algorithm, test oracle only."""
    Z, residual = X-X.mean(axis=0), y-y.mean()
    curvature = np.mean(Z*Z, axis=0)
    beta = np.zeros(X.shape[1])
    output = []
    for lam in lambdas:
        l1, l2 = lam*alpha, lam*(1-alpha)
        for _ in range(100000):
            for j in range(len(beta)):
                denominator = curvature[j]+l2
                new = soft(Z[:, j]@residual/len(y)+curvature[j]*beta[j]+l2*neighbor(W, beta, j), l1)/denominator if denominator else 0.
                residual += (beta[j]-new)*Z[:, j]
                beta[j] = new
            gradient = -Z.T@residual/len(y)+l2*(beta-W@beta)
            if kkt_residual(gradient, beta, l1) <= 1e-11:
                break
        else:
            raise AssertionError("exhaustive test oracle did not converge")
        output.append(beta.copy())
    return np.array(output).T


@pytest.mark.parametrize("seed", [7, 8, 9])
@pytest.mark.parametrize("alpha", [0., .05, .6, 1.])
def test_matches_exhaustive_coordinates(seed, alpha):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(30, 14))*.4
    # Connected graph means previously inactive cells can become active due to
    # changed neighbors; the KKT scan must still consider every cell each pass.
    graph = sparse.diags([np.ones(13), np.ones(13)], [-1, 1], format="csr")
    y = X[:, 2]*2-X[:, 7]+.1*rng.normal(size=30)
    lambdas = [.3, .1, .03, .01]
    actual = fit_path(X, y, graph, family="gaussian", semantics="corrected",
                      alpha=alpha, lambdas=lambdas, thresh=1e-11)
    expected = exhaustive(X, y, normalize_graph(graph), lambdas, alpha)
    assert actual.converged.all()
    np.testing.assert_allclose(actual.coefficients, expected, atol=1e-8, rtol=1e-7)
    assert max(actual.metadata["kkt_residual"]) <= 1e-11
