"""Independent mathematical checks of the linear-memory Cox implementation."""
import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy import sparse
from scipy.optimize import minimize
from scipy.special import logsumexp

from scissor_py._cox import cox_lambda_max, cox_log_likelihood, cox_path, cox_score


def _example():
    rng = np.random.default_rng(106)
    x = rng.normal(size=(18, 4))
    y = np.column_stack((np.repeat(np.arange(1, 10), 2),
                         [0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1]))
    graph = sparse.csr_matrix(np.array([[0, .5, 0, .5], [.5, 0, .5, 0],
                                        [0, .5, 0, .5], [.5, 0, .5, 0]]))
    return x, y, graph


def _reference(x, y, beta):
    eta = x @ beta
    result = 0.0
    score = np.zeros(x.shape[1])
    for time in np.unique(y[y[:, 1] == 1, 0]):
        event = (y[:, 0] == time) & (y[:, 1] == 1)
        risk = y[:, 0] >= time
        count = int(event.sum())
        denominator = logsumexp(eta[risk])
        result += np.sum(eta[event] - denominator)
        score += x[event].sum(axis=0) - count * (np.exp(eta[risk] - denominator) @ x[risk])
    return result, score


def test_likelihood_and_score_match_independent_breslow_reference():
    x, y, _ = _example()
    beta = np.array([.2, -.5, .1, .9])
    ll, score = _reference(x, y, beta)
    assert_allclose(cox_log_likelihood(x, y, beta), ll, atol=1e-12)
    assert_allclose(cox_score(x, y, beta), score, atol=1e-12)
    h = 1e-6
    finite_difference = [
        (cox_log_likelihood(x, y, beta + h * direction)
         - cox_log_likelihood(x, y, beta - h * direction)) / (2 * h)
        for direction in np.eye(x.shape[1])
    ]
    assert_allclose(cox_score(x, y, beta), finite_difference, rtol=2e-8, atol=2e-8)


def test_tied_censoring_stays_in_risk_set():
    x = np.array([[0.0], [1.0], [2.0], [3.0]])
    y = np.array([[1, 0], [2, 1], [2, 0], [3, 1]])
    beta = np.array([0.0])
    assert_allclose(cox_log_likelihood(x, y, beta), -np.log(3))
    assert_allclose(cox_score(x, y, beta), [-1])
    assert_allclose(cox_lambda_max(x, y, .5), .5)
    assert_allclose(cox_lambda_max(x, y, 0), 250)


def test_permutation_location_and_leading_censor_invariance():
    x, y, _ = _example()
    beta = np.array([.2, -.5, .1, .9])
    order = np.random.default_rng(99).permutation(len(x))
    assert_allclose(cox_log_likelihood(x[order] + 10000, y[order], beta),
                    cox_log_likelihood(x, y, beta), atol=1e-11)
    assert_allclose(cox_score(x[order] + 10000, y[order], beta),
                    cox_score(x, y, beta), atol=1e-10)
    assert_allclose(cox_log_likelihood(x[2:], y[2:], beta), cox_log_likelihood(x, y, beta))
    assert_allclose(cox_lambda_max(x, y, .5) * len(x),
                    cox_lambda_max(x[2:], y[2:], .5) * (len(x) - 2))


def test_extreme_linear_predictors_remain_finite():
    x = np.array([[2000.0], [1000.0], [0.0], [-1000.0], [-2000.0]])
    y = np.column_stack((np.arange(5), np.ones(5)))
    for beta in (np.array([1.0]), np.array([-1.0])):
        ll, score = _reference(x, y, beta)
        assert np.isfinite(cox_log_likelihood(x, y, beta))
        assert np.isfinite(cox_score(x, y, beta)).all()
        assert_allclose(cox_log_likelihood(x, y, beta), ll, atol=1e-10)
        assert_allclose(cox_score(x, y, beta), score, atol=1e-9)


def test_corrected_path_matches_independent_bound_constrained_optimizer():
    x, y, graph = _example()
    alpha = .37
    lambdas = np.array([.18, .06, .012])
    betas, intercepts, converged, _, metadata = cox_path(
        x, y, graph, lambdas, alpha, thresh=1e-9, maxit=10000,
    )
    assert converged.all(), metadata
    assert_allclose(intercepts, 0)
    assert metadata["n_original"] == 18
    assert metadata["n_retained"] == 16
    assert np.max(metadata["kkt_residual"]) <= 1e-9
    for k, lam in enumerate(lambdas):
        l1, l2 = lam * alpha, lam * (1 - alpha)

        def objective(uv):
            beta = uv[:4] - uv[4:]
            ll, score = _reference(x, y, beta)
            lap = beta - graph @ beta
            value = -ll / len(x) + l1 * uv.sum() + .5 * l2 * (beta @ lap)
            grad = -score / len(x) + l2 * lap
            return value, np.r_[grad + l1, -grad + l1]

        fit = minimize(objective, np.zeros(8), jac=True, bounds=[(0, None)] * 8,
                       method="L-BFGS-B", options={"ftol": 1e-14, "gtol": 1e-11, "maxiter": 10000})
        assert fit.success
        assert_allclose(betas[:, k], fit.x[:4] - fit.x[4:], rtol=3e-5, atol=3e-7)
        assert_allclose(metadata["objective"][k], fit.fun, atol=1e-10)


def test_null_lambda_is_zero_and_no_dense_graph_conversion(monkeypatch):
    x, y, graph = _example()
    def forbid(*args, **kwargs):
        raise AssertionError("The Cox kernel must not densify its sparse graph")
    monkeypatch.setattr(sparse.csr_matrix, "toarray", forbid)
    largest = cox_lambda_max(x, y, .5)
    result = cox_path(x, y, graph, [largest * 1.0001, largest * .7], .5)
    assert_allclose(result[0][:, 0], 0)
    assert result[2].all()


def test_report_maxit_without_claiming_convergence():
    x, y, graph = _example()
    result = cox_path(x, y, graph, [.001], .5, maxit=1, thresh=1e-14)
    assert not result[2][0]
    assert result[3][0] == 1
    assert result[4]["termination_reason"] == ["maxit"]


def test_official_coordinate_mode_near_same_optimum_with_tight_threshold():
    x, y, graph = _example()
    corrected = cox_path(x, y, graph, [.2, .08], .6, thresh=1e-10)
    official = cox_path(x, y, graph, [.2, .08], .6, semantics="official",
                        thresh=1e-13, ifast=False)
    assert official[2].all()
    assert official[4]["effective_path_length"] == 2
    assert_allclose(official[0], corrected[0], atol=2e-6)


def test_no_events_likelihood_is_zero_but_fit_is_rejected():
    x, y, graph = _example()
    y[:, 1] = 0
    assert cox_log_likelihood(x, y, np.ones(4)) == 0
    assert_allclose(cox_score(x, y, np.ones(4)), 0)
    with pytest.raises(ValueError, match="observed event"):
        cox_path(x, y, graph, [.1], .5)


@pytest.mark.parametrize("bad_y", [np.array([[1, .2], [2, 1]]),
                                     np.array([[-1, 1], [2, 1]]),
                                     np.array([[np.nan, 1], [2, 1]])])
def test_invalid_outcomes_are_rejected(bad_y):
    with pytest.raises(ValueError):
        cox_log_likelihood(np.ones((2, 1)), bad_y, [0])


@pytest.mark.parametrize("tolerance", ["default", "tight"])
def test_official_path_matches_pinned_actual_r_fixture(tolerance):
    import json
    from pathlib import Path

    path = (Path(__file__).parents[1] / "reference" / "fixtures"
            / f"cox_{tolerance}_explicit_cv.json")
    with path.open() as stream:
        fixture = json.load(stream)
    inputs = fixture["inputs"]
    adjacency = np.pad(np.asarray(inputs["Omega"], dtype=float), ((1, 0), (1, 0)))
    degree = adjacency.sum(axis=0)
    inverse_sqrt = np.zeros_like(degree)
    np.divide(1, np.sqrt(degree), out=inverse_sqrt, where=degree > 0)
    normalized = inverse_sqrt[:, None] * adjacency * inverse_sqrt[None, :]
    # The caller, not the low-level Cox kernel, applies upstream's intercept
    # padding/index-alignment semantics. These fixtures deliberately avoid the
    # source's out-of-range final-node adjacency case.
    graph = sparse.csr_matrix(normalized[:-1, :-1])
    result = cox_path(inputs["x"], inputs["y"], graph, inputs["lambda"],
                      inputs["alpha"], semantics="official", thresh=inputs["thresh"],
                      maxit=inputs["maxit"])
    assert fixture["provenance"]["modified_upstream"] is False
    assert result[2].all()
    assert_allclose(result[0], fixture["output"]["result"]["Beta"], atol=2e-11, rtol=2e-11)
    assert_allclose(cox_lambda_max(inputs["x"], inputs["y"], inputs["alpha"]),
                    .2368093421544794, atol=1e-14)


def test_score_uses_linear_memory_for_many_risk_sets():
    import tracemalloc

    rng = np.random.default_rng(712)
    n = 6000
    x = rng.normal(size=(n, 3))
    y = np.column_stack((np.arange(n), np.ones(n)))
    tracemalloc.start()
    try:
        score = cox_score(x, y, np.array([.1, -.3, .2]))
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert np.isfinite(score).all()
    # An n-by-n float risk matrix alone would require 288 MB here. This
    # generous limit permits implementation bookkeeping and BLAS temporaries.
    assert peak < 15_000_000


def test_only_last_subject_event_has_zero_information():
    x = np.arange(12, dtype=float).reshape(4, 3)
    y = np.array([[1, 0], [2, 0], [3, 0], [4, 1]])
    result = cox_path(x, y, None, [.1, 0], .5)
    assert_allclose(result[0], 0)
    assert result[2].all()
    assert result[4]["n_retained"] == 1
    assert cox_log_likelihood(x, y, np.ones(3)) == 0


def test_all_tied_events_have_breslow_count_factor():
    x = np.array([[-2.], [-1.], [0.], [1.], [2.]])
    y = np.column_stack((np.ones(5), [1, 1, 0, 1, 0]))
    assert_allclose(cox_log_likelihood(x, y, [0]), -3 * np.log(5))
    assert_allclose(cox_score(x, y, [0]), [-2])
