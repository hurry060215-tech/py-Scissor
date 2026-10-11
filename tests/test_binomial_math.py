import numpy as np
from scipy import sparse
from scipy.optimize import minimize
from scipy.special import expit
from scissor_py import fit_path, normalize_graph


def test_corrected_binomial_matches_independent_optimizer():
    rng = np.random.default_rng(82)
    X = rng.normal(size=(40, 5))
    y = (rng.random(40) < expit(X[:, 0]-X[:, 1]-.8)).astype(float)
    graph = sparse.diags([np.ones(4), np.ones(4)], [-1, 1], format="csr")
    W = normalize_graph(graph)
    result = fit_path(X, y, graph, family="binomial", alpha=.4, lambdas=[.1, .03], semantics="corrected", thresh=1e-9)
    assert result.converged.all()
    for k, lam in enumerate(result.lambdas):
        def objective(v):
            b = v[1:6]-v[6:]
            eta = v[0]+X@b
            lap = b-W@b
            loss = np.logaddexp(0, (1-2*y)*eta).mean()+lam*.4*v[1:].sum()+lam*.6/2*(b@lap)
            score = expit(eta)-y
            gradient = X.T@score/len(y)+lam*.6*lap
            return loss, np.r_[score.mean(), gradient+lam*.4, -gradient+lam*.4]
        oracle = minimize(objective, np.zeros(11), jac=True, bounds=[(None,None)]+[(0,None)]*10,
                           method="L-BFGS-B", options={"ftol":1e-14,"gtol":1e-11,"maxiter":10000})
        assert oracle.success
        np.testing.assert_allclose(result.coefficients[:,k], oracle.x[1:6]-oracle.x[6:], atol=1e-6)
        np.testing.assert_allclose(result.intercepts[k], oracle.x[0], atol=1e-6)
        assert result.metadata["kkt_residual"][k] <= 1e-9


def test_unbalanced_classes_and_extreme_predictors():
    rng = np.random.default_rng(82)
    X = rng.normal(size=(40, 4))*100
    y = np.r_[1.,1.,np.zeros(38)]
    out = fit_path(X, y, sparse.csr_matrix((4,4)), family="binomial", alpha=.5,
                    lambdas=[.1], semantics="corrected", thresh=1e-8, maxit=30000)
    assert out.converged.all(), out.metadata
    assert np.isfinite(out.coefficients).all()
    assert out.metadata["probability_clip"] is None
