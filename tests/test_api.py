import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scissor_py import fit_path, cross_validate, scissor, scissor_from_anndata, ScissorError, ConvergenceError
from scissor_py import prepare


def problem():
    rng = np.random.default_rng(92)
    X = rng.normal(size=(36, 6))*.3
    y = 1+1.2*X[:, 0]-.6*X[:, 1]+rng.normal(size=36)*.03
    graph = sparse.diags([np.ones(5), np.ones(5)], [-1, 1], format="csr")
    return X, y, graph


def test_corrected_gaussian_permutation_and_kkt():
    X, y, graph = problem()
    kwargs = dict(alpha=.4, lambdas=[.1, .03, .01], semantics="corrected", thresh=1e-10)
    a = fit_path(X, y, graph, **kwargs)
    order = [3, 1, 5, 0, 4, 2]
    b = fit_path(X[:, order], y, graph[order, :][:, order], **kwargs)
    np.testing.assert_allclose(a.coefficients[order], b.coefficients, atol=1e-8)
    assert max(a.metadata["kkt_residual"]) <= 1e-10
    assert a.converged.all()


def test_all_zero_selection_and_nonconvergence():
    X, y, graph = problem()
    zero = fit_path(X, y, graph, semantics="corrected", lambdas=[100])
    assert not np.count_nonzero(zero.coefficients)
    bad = fit_path(X, y, graph, semantics="corrected", lambdas=[.001], thresh=1e-15, maxit=1)
    assert not bad.converged[0]
    with pytest.raises(ConvergenceError):
        cross_validate(X, y, graph, semantics="corrected", lambdas=[.001], thresh=1e-15, maxit=1, nfolds=3)


def test_zero_score_auto_lambda_path():
    X, _, graph = problem()
    fit = fit_path(X, np.ones(len(X)), graph, semantics="corrected", nlambda=3)
    assert np.all(fit.lambdas == 0)
    assert np.all(fit.coefficients == 0)
    assert np.all(fit.intercepts == 1)


@pytest.mark.parametrize("family", ["gaussian", "binomial", "cox"])
def test_invalid_folds(family):
    X, y, graph = problem()
    if family == "binomial": y = np.r_[np.zeros(18), np.ones(18)]
    if family == "cox": y = np.column_stack([np.arange(1, 37), np.r_[np.zeros(18), np.ones(18)]])
    fold = np.r_[np.zeros(18, int), np.ones(18, int)]
    if family == "gaussian": fold[-1] = 4
    with pytest.raises(ScissorError):
        cross_validate(X, y, graph, family=family, semantics="corrected", foldid=fold)


def test_high_level_ids_labels_cutoff_and_no_mutation():
    rng = np.random.default_rng(8)
    bulk = pd.DataFrame(rng.poisson(6, size=(30, 18)), index=[f"g{i}" for i in range(30)], columns=[f"s{i}" for i in range(18)])
    sc = pd.DataFrame(rng.poisson(6, size=(30, 5)), index=bulk.index, columns=[f"c{i}" for i in range(5)])
    graph = pd.DataFrame(np.zeros((5, 5)), index=sc.columns, columns=sc.columns)
    y = pd.Series(["case", "control"]*9, index=bulk.columns).iloc[::-1]
    with pytest.raises(ScissorError, match="positive_label"):
        scissor(bulk, sc, y, graph, family="binomial", alpha=.5, lambdas=[1], nfolds=3)
    result = scissor(bulk, sc, y, graph, family="binomial", alpha=.5, lambdas=[1], nfolds=3, positive_label="case")
    assert result.positive_cells == result.negative_cells == []
    assert result.cutoff_reached
    assert result.metadata["positive_class"] == {"positive": "case", "negative": "control"}
    assert tuple(result.coefficients.index) == tuple(sc.columns)
    assert all(result.labels == "Background")
    class FakeAnnData:
        X = sc.T.to_numpy(); layers = {"counts": sparse.csr_matrix(sc.T)}
        obsp = {"network": sparse.csr_matrix(graph)}; var_names = sc.index; obs_names = sc.columns
    adapted = scissor_from_anndata(FakeAnnData(), bulk, y, layer="counts", graph_key="network", family="binomial", alpha=.5, lambdas=[1], nfolds=3, positive_label="case")
    pd.testing.assert_series_equal(result.coefficients, adapted.coefficients)


def test_cutoff_unmet_explicit_flag():
    rng = np.random.default_rng(8)
    bulk = pd.DataFrame(rng.normal(size=(50, 20)), columns=[f"s{i}" for i in range(20)])
    sc = pd.DataFrame(rng.normal(size=(50, 3)), columns=[f"c{i}" for i in range(3)])
    y = pd.Series(np.arange(20), index=bulk.columns)
    graph = pd.DataFrame(np.zeros((3, 3)), index=sc.columns, columns=sc.columns)
    with pytest.warns(RuntimeWarning, match="cutoff"):
        r = scissor(bulk, sc, y, graph, alpha=[.3, .5], lambdas=[1e-5], cutoff=.01, semantics="corrected", nfolds=4)
    assert not r.cutoff_reached
    assert r.alpha == .5


@pytest.mark.parametrize("family,seed,maxit", [("binomial",0,10),("gaussian",23,30)])
def test_active_set_restart_does_not_reuse_stale_convergence(family, seed, maxit):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(20,5)) * (1 if family == "binomial" else .3)
    X[:,1] = .6*X[:,0]+.4*X[:,1]
    y = rng.binomial(1,.5,20).astype(float) if family == "binomial" else rng.normal(size=20)
    result = fit_path(X,y,np.zeros((5,5)),family=family,alpha=.9,lambdas=[.03],thresh=1e-5,maxit=maxit)
    assert not result.converged[0]


def test_nonunique_lasso_compares_predictions_not_individual_coefficients():
    x = np.linspace(-1,1,30)
    X = np.column_stack([x,x])
    y = 2*x
    a = fit_path(X,y,np.zeros((2,2)),alpha=1,lambdas=[.05],semantics="corrected",thresh=1e-12)
    b = fit_path(X[:,::-1],y,np.zeros((2,2)),alpha=1,lambdas=[.05],semantics="corrected",thresh=1e-12)
    np.testing.assert_allclose(X@a.coefficients, X[:,::-1]@b.coefficients, atol=1e-12)
    assert np.max(abs(a.coefficients-b.coefficients[::-1])) > 1
    assert a.converged.all() and b.converged.all()


def test_official_single_cell_cox_wrapper_failure_boundary():
    X = np.arange(12).reshape(12,1)/10
    y = np.column_stack([np.arange(1,13),np.ones(12)])
    with pytest.raises(ScissorError,match="p=1"):
        fit_path(X,y,[[0]],family="cox",lambdas=[.1])
    assert fit_path(X,y,[[0]],family="cox",lambdas=[.1],semantics="corrected").converged.all()


@pytest.mark.parametrize("random_state", [None, np.random.default_rng(123)])
def test_folds_reused_across_alphas_with_unseeded_or_generator_rng(random_state):
    rng = np.random.default_rng(8)
    bulk = pd.DataFrame(rng.normal(size=(50,20)),columns=[f"s{i}" for i in range(20)])
    sc = pd.DataFrame(rng.normal(size=(50,3)),columns=[f"c{i}" for i in range(3)])
    y = pd.Series(np.arange(20),index=bulk.columns)
    graph = pd.DataFrame(np.zeros((3,3)),index=sc.columns,columns=sc.columns)
    with pytest.warns(RuntimeWarning,match="cutoff"):
        r = scissor(bulk,sc,y,graph,alpha=[.3,.5],lambdas=[1e-5],cutoff=.01,
                    semantics="corrected",nfolds=4,random_state=random_state)
    np.testing.assert_array_equal(r.cv_results[0].foldid,r.cv_results[1].foldid)


def test_near_zero_coefficients_remain_strict_sign_labels():
    rng = np.random.default_rng(18)
    bulk = pd.DataFrame(rng.normal(size=(30,20)),columns=[f"s{i}" for i in range(20)])
    sc = pd.DataFrame(rng.normal(size=(30,3)),columns=[f"c{i}" for i in range(3)])
    graph = pd.DataFrame(np.zeros((3,3)),index=sc.columns,columns=sc.columns)
    X = prepare(bulk,sc,graph).similarity
    y = pd.Series(1e-7*X[:,0],index=bulk.columns)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",RuntimeWarning)
        r = scissor(bulk,sc,y,graph,alpha=1,lambdas=[0],cutoff=.99,
                    semantics="corrected",thresh=1e-14,nfolds=4)
    assert 0 < r.coefficients["c0"] < 1e-5
    assert r.labels["c0"] == "Scissor+"
    assert "c0" in r.metadata["near_zero_cells"]
