"""Public Scissor API. Explicit numerical semantics; R-free runtime."""
import warnings
import numpy as np
import pandas as pd
from scipy.special import expit
from ._types import PathResult, CVResult, ScissorResult, ScissorError, ConvergenceError, REFERENCE_COMMIT
from ._preprocessing import prepare, normalize_graph, _validate_dense_budget
from ._validation import validate_problem, align_phenotype, make_folds
from ._solvers import gaussian_path, binomial_path, wrapper_network

DEFAULT_ALPHAS = (.005, .01, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9)


def _lambda_path(X, y, family, alpha, lambdas, nlambda, rlambda):
    if lambdas is not None:
        if np.iscomplexobj(lambdas):
            raise ScissorError("lambdas must be real-valued")
        path = np.atleast_1d(np.asarray(lambdas, dtype=float))
        if path.ndim != 1 or len(path) == 0 or not np.isfinite(path).all() or np.any(path < 0) or np.any(np.diff(path) > 0):
            raise ScissorError("lambdas must be a nonempty finite nonnegative descending vector")
        return path, path.copy(), False
    if not isinstance(nlambda, (int, np.integer)) or nlambda < 2:
        raise ScissorError("automatic lambda path requires nlambda >= 2")
    ratio = (1e-4 if X.shape[0] > X.shape[1] else 1e-2) if rlambda is None else rlambda
    if np.ndim(ratio) != 0 or np.iscomplexobj(ratio) or not np.isfinite(ratio) or not 0 < ratio < 1:
        raise ScissorError("rlambda must be between 0 and 1")
    ratios = ratio**(np.arange(nlambda)/(nlambda-1))
    a = alpha if alpha > 0 else .001
    if family in ("gaussian", "binomial"):
        maximum = np.max(np.abs((X-X.mean(axis=0)).T@(y-y.mean())))/len(y)/a
    else:
        from ._cox import cox_lambda_max
        maximum = cox_lambda_max(X, y, alpha)
    actual = ratios*maximum
    return actual, ratios, True


def _fit_path(X, y, network, *, family, alpha, semantics, lambdas=None,
              nlambda=100, rlambda=None, thresh=1e-7, maxit=100000, ifast=True, cv=False, engine="numpy"):
    X, y, W = validate_problem(X, y, network, family, alpha, semantics, thresh, maxit)
    actual, ratios, automatic = _lambda_path(X, y, family, alpha, lambdas, nlambda, rlambda)
    W = wrapper_network(W, family, semantics, alpha)
    if engine not in ("numpy", "numba"):
        raise ScissorError("engine must be 'numpy' or 'numba'")
    if engine == "numba":
        if family != "gaussian" or semantics != "corrected":
            raise ScissorError("engine='numba' currently supports corrected Gaussian only")
        try:
            from ._gaussian_numba import gaussian_numba_path
        except ModuleNotFoundError as error:
            if error.name == "numba":
                raise ScissorError("engine='numba' requires the optional scissor-py[acceleration] dependency") from error
            raise
        out = gaussian_numba_path(X, y, W, actual, alpha, thresh=thresh, maxit=maxit)
    elif family == "gaussian":
        out = gaussian_path(X, y, W, actual, alpha, semantics=semantics, thresh=thresh, maxit=maxit, cv=cv)
    elif family == "binomial":
        out = binomial_path(X, y, W, actual, alpha, semantics=semantics, thresh=thresh, maxit=maxit)
    else:
        from ._cox import cox_path
        out = cox_path(X, y, W, actual, alpha, semantics=semantics, thresh=thresh, maxit=maxit, ifast=ifast)
    coefficients, intercepts, converged, iterations, metadata = out
    length = min(coefficients.shape[1], metadata.get("effective_path_length", coefficients.shape[1]))
    reported = ratios if family == "cox" and semantics == "official" and automatic else actual
    finite = np.isfinite(coefficients[:, :length]).all(axis=0) & np.isfinite(intercepts[:length])
    metadata.update({"semantics": semantics, "reference_commit": REFERENCE_COMMIT, "engine": engine,
                     "solver_lambdas": actual[:length].tolist(), "automatic_lambda": automatic,
                     "thresh": thresh, "maxit": maxit, "ifast": ifast,
                     "graph_alignment": "published padded graph" if family != "binomial" and semantics == "official" and alpha != 1 else "cell aligned",
                     "cox_lambda_label_mismatch": bool(family == "cox" and semantics == "official" and automatic),
                     "unused_L0_trimming": "not implemented; Scissor uses untrimmed Beta and lambda.min"})
    return PathResult(family, alpha, reported[:length].copy(), coefficients[:, :length],
                      intercepts[:length], converged[:length] & finite, iterations[:length], metadata)


def fit_path(X, y, graph, *, family="gaussian", alpha=.1, semantics="official",
             lambdas=None, nlambda=100, rlambda=None, thresh=1e-7, maxit=100000, ifast=True, engine="numpy"):
    """Fit an APML1-Net path from samples × cells similarity and a raw cell graph.

    y is numeric n-vector, or n×2 [time, status] for Cox. Inputs are positional
    here; use scissor() for ID-aware expression inputs. Graph is binarized and
    degree normalized. No expression preprocessing is applied at this layer.

    official reproduces defined published behavior, including documented
    path/stopping quirks, and rejects unsafe wrapper graph indexing. corrected
    explicitly changes these behaviors; it is not strict official equivalence.
    """
    return _fit_path(X, y, normalize_graph(graph), family=family, alpha=alpha, semantics=semantics,
                     lambdas=lambdas, nlambda=nlambda, rlambda=rlambda, thresh=thresh, maxit=maxit, ifast=ifast, engine=engine)


def _cross_validate(X, y, network, *, family, alpha, semantics, foldid=None,
                     nfolds=10, random_state=123, lambdas=None, nlambda=100,
                     rlambda=None, thresh=1e-7, maxit=100000, engine="numpy"):
    X, y, network = validate_problem(X, y, network, family, alpha, semantics, thresh, maxit)
    folds = make_folds(y, family, foldid, nfolds, random_state)
    full = _fit_path(X, y, network, family=family, alpha=alpha, semantics=semantics,
                     lambdas=lambdas, nlambda=nlambda, rlambda=rlambda, thresh=thresh, maxit=maxit, engine=engine)
    levels = np.unique(folds)
    count = len(full.lambdas)
    if not count:
        raise ConvergenceError("published path truncation returned no fit")
    losses = np.full((len(levels), count), np.nan)
    valid = full.converged.copy()
    weights = np.empty(len(levels))
    fold_diagnostics = []
    for k, fold in enumerate(levels):
        train, held = folds != fold, folds == fold
        sub = _fit_path(X[train], y[train], network, family=family, alpha=alpha, semantics=semantics,
                        lambdas=full.lambdas, thresh=thresh, maxit=maxit, ifast=False, cv=True, engine=engine)
        q = len(sub.lambdas)
        good = np.zeros(count, bool); good[:q] = sub.converged
        valid &= good
        if family == "gaussian":
            prediction = X[held]@sub.coefficients+sub.intercepts
            losses[k, :q] = np.mean((y[held, None]-prediction)**2, axis=0)
            weights[k] = held.sum()
        elif family == "binomial":
            probability = np.clip(expit(X[held]@sub.coefficients+sub.intercepts), 1e-5, 1-1e-5)
            losses[k, :q] = -2*np.mean(y[held, None]*np.log(probability)+(1-y[held, None])*np.log1p(-probability), axis=0)
            weights[k] = held.sum()
        else:
            from ._cox import cox_log_likelihood
            weights[k] = y[held, 1].sum()
            for j in range(q):
                b = sub.coefficients[:, j]
                losses[k, j] = -(cox_log_likelihood(X, y, b)-cox_log_likelihood(X[train], y[train], b))/weights[k]
        fold_diagnostics.append({"fold": int(fold), "converged": sub.converged.tolist(),
                                 "iterations": sub.iterations.tolist()})
    valid &= np.isfinite(losses).all(axis=0)
    means = np.average(losses, axis=0, weights=weights)
    se = np.sqrt(np.average((losses-means)**2, axis=0, weights=weights)/(len(levels)-1))
    if not np.any(valid):
        raise ConvergenceError("no lambda converged on the full data and every validation fold")
    best = int(np.argmin(np.where(valid, means, np.inf)))
    one_se = int(np.flatnonzero(valid & (means <= means[best]+se[best]))[0])
    return CVResult(full, folds, means, se, losses, float(full.lambdas[best]),
                    float(full.lambdas[one_se]), valid,
                    {"semantics": semantics, "fold_weights": weights.tolist(), "folds": fold_diagnostics,
                     "loss": {"gaussian": "held-out mean squared error", "binomial": "held-out binomial deviance",
                               "cox": "negative full-minus-training Breslow partial log likelihood per held-out event"}[family],
                     "convergence_policy": "reject a candidate if full fit or any fold did not converge",
                     "rng": "NumPy default_rng; equal R and Python seeds do not imply identical folds"})


def cross_validate(X, y, graph, *, family="gaussian", alpha=.1, semantics="official", foldid=None,
                    nfolds=10, random_state=123, lambdas=None, nlambda=100,
                    rlambda=None, thresh=1e-7, maxit=100000, engine="numpy"):
    """Fixed-graph CV with explicit folds and family-specific published losses."""
    return _cross_validate(X, y, normalize_graph(graph), family=family, alpha=alpha,
                             semantics=semantics, foldid=foldid, nfolds=nfolds,
                             random_state=random_state, lambdas=lambdas, nlambda=nlambda,
                             rlambda=rlambda, thresh=thresh, maxit=maxit, engine=engine)


def scissor(bulk, single_cell, phenotype, graph, *, family="gaussian", alpha=None,
            cutoff=.2, semantics="official", graph_cell_ids=None, positive_label=None,
            foldid=None, nfolds=None, random_state=123, lambdas=None, nlambda=100,
            rlambda=None, thresh=1e-7, maxit=100000, max_dense_bytes=1_000_000_000, engine="numpy"):
    """Select phenotype-associated cells from genes × observations DataFrames.

    Joint quantile normalization and Pearson similarity follow Scissor. Graph
    construction, count/log transforms, and batch correction remain explicit
    upstream steps. The first alpha with selected_fraction < cutoff is used,
    otherwise the last successful candidate is returned with cutoff_reached=False.
    At each alpha, CV lambda.min is refitted from zero as in Scissor.R.
    """
    if np.ndim(cutoff) != 0 or np.iscomplexobj(cutoff) or not np.isfinite(cutoff) or not 0 < cutoff <= 1:
        raise ScissorError("cutoff must lie in (0, 1]")
    if np.iscomplexobj(alpha):
        raise ScissorError("alpha must be real-valued")
    alphas = np.asarray(DEFAULT_ALPHAS) if alpha is None else np.atleast_1d(alpha).astype(float)
    if alphas.ndim != 1 or len(alphas) == 0 or not np.isfinite(alphas).all() or np.any(alphas < 0) or np.any(alphas > 1):
        raise ScissorError("alpha must be a nonempty scalar/vector within [0, 1]")
    prepared = prepare(bulk, single_cell, graph, graph_cell_ids=graph_cell_ids, max_dense_bytes=max_dense_bytes)
    y, class_mapping = align_phenotype(phenotype, prepared.sample_ids, family, positive_label)
    if isinstance(foldid, pd.Series):
        if not foldid.index.is_unique or set(foldid.index) != set(prepared.sample_ids):
            raise ScissorError("foldid Series IDs must equal bulk sample IDs")
        foldid = foldid.loc[list(prepared.sample_ids)].to_numpy()
    nfolds = min(10, len(y)) if nfolds is None else nfolds
    fixed_folds = make_folds(y, family, foldid, nfolds, random_state)
    X, network = prepared.similarity, prepared.network
    results, failures, chosen = [], [], None
    reached = False
    for a in alphas:
        try:
            cv = _cross_validate(X, y, network, family=family, alpha=float(a), semantics=semantics,
                                  foldid=fixed_folds, nfolds=nfolds, random_state=random_state,
                                  lambdas=lambdas, nlambda=nlambda, rlambda=rlambda, thresh=thresh, maxit=maxit, engine=engine)
            refit = _fit_path(X, y, network, family=family, alpha=float(a), semantics=semantics,
                              lambdas=[cv.lambda_min], thresh=thresh, maxit=maxit, engine=engine)
            if len(refit.lambdas) != 1 or not refit.converged[0]:
                raise ConvergenceError("lambda.min refit did not converge")
        except ConvergenceError as error:
            failures.append({"alpha": float(a), "error": str(error)})
            continue
        results.append(cv)
        coefficients = refit.coefficients[:, 0]
        reached = bool(np.count_nonzero(coefficients)/len(coefficients) < cutoff)
        chosen = (float(a), cv, refit, coefficients.copy())
        if reached:
            break
    if chosen is None:
        raise ConvergenceError(f"all alpha candidates failed: {failures}")
    a, cv, refit, coefficients = chosen
    labels = np.where(coefficients > 0, "Scissor+", np.where(coefficients < 0, "Scissor-", "Background"))
    near_zero = [prepared.cell_ids[i] for i in np.flatnonzero((coefficients != 0) & (np.abs(coefficients) < 1e-5))]
    if not reached:
        warnings.warn("No successful alpha met cutoff; returning last successful alpha with cutoff_reached=False", RuntimeWarning, stacklevel=2)
    return ScissorResult(pd.Series(coefficients, index=prepared.cell_ids, name="coefficient"),
                         pd.Series(labels, index=prepared.cell_ids, name="scissor_label"),
                         a, cv.lambda_min, float(refit.intercepts[0]), family, reached,
                         results, prepared,
                         {"semantics": semantics, "reference_commit": REFERENCE_COMMIT,
                          "positive_class": class_mapping, "cutoff": cutoff,
                          "requested_alphas": [float(value) for value in alphas],
                          "lambda_mode": "automatic" if lambdas is None else "explicit",
                          "requested_nlambda": nlambda if lambdas is None else None,
                          "requested_rlambda": rlambda,
                          "fold_source": "generated" if foldid is None else "provided",
                          "random_state": random_state,
                          "selected_fraction": float(np.count_nonzero(coefficients)/len(coefficients)),
                          "failed_alpha_candidates": failures, "near_zero_cells": near_zero,
                          "near_zero_threshold_diagnostic_only": 1e-5,
                          "label_rule": "strict sign; no tolerance-based relabeling",
                          "non_unique_solutions": "not automatically certified; compare predictions/objectives under permutations",
                          "refit": refit.metadata})


def scissor_from_anndata(adata, bulk, phenotype, *, layer, graph_key,
                         writeback=False, key_added="scissor", overwrite=False, **kwargs):
    """Optional duck-typed AnnData adapter, with explicit expression layer/graph.

    layer='X' means adata.X; all other names index adata.layers. This adapter
    is nonmutating by default and does not call Omicverse, Scanpy, or R.
    Explicit writeback=True annotates an actual in-memory AnnData only after
    fitting/serialization succeed; collisions require overwrite=True.
    Expression is materialized only after checking the configured byte budget.
    """
    from scipy import sparse
    if not isinstance(writeback, (bool, np.bool_)):
        raise ScissorError("writeback must be boolean")
    if writeback:
        from ._anndata import _check_writeback
        _check_writeback(adata, key_added, overwrite)
    if not isinstance(layer, str) or not isinstance(graph_key, str):
        raise ScissorError("layer and graph_key must be explicitly supplied strings")
    if graph_key not in adata.obsp:
        raise ScissorError(f"graph_key {graph_key!r} is missing from adata.obsp")
    if layer != "X" and layer not in adata.layers:
        raise ScissorError(f"layer {layer!r} is missing from adata.layers")
    values = adata.X if layer == "X" else adata.layers[layer]
    estimated = 8*values.shape[0]*values.shape[1]
    limit = kwargs.get("max_dense_bytes", 1_000_000_000)
    _validate_dense_budget(limit)
    if estimated > limit:
        raise ScissorError(f"dense expression alone needs {estimated} bytes, above max_dense_bytes={limit}")
    # AnnData's backed CSR/CSC dataset is not a scipy sparse matrix. Its public
    # to_memory() method materializes sparse storage; only then densify under
    # the already checked byte budget. No AnnData import or file write occurs.
    if not sparse.issparse(values) and callable(getattr(values, "to_memory", None)):
        values = values.to_memory()
    if sparse.issparse(values):
        values = values.toarray()
    if np.asarray(values).shape != (len(adata.obs_names), len(adata.var_names)):
        raise ScissorError("selected expression layer shape must match AnnData obs/var IDs")
    expression = pd.DataFrame(np.asarray(values).T, index=adata.var_names, columns=adata.obs_names)
    result = scissor(bulk, expression, phenotype, adata.obsp[graph_key], graph_cell_ids=adata.obs_names, **kwargs)
    result.metadata["anndata_layer"] = layer
    result.metadata["anndata_graph_key"] = graph_key
    if writeback:
        from ._anndata import annotate_anndata
        annotate_anndata(adata, result, key_added=key_added, overwrite=overwrite)
    return result
