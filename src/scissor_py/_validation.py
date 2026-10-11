import numpy as np
import pandas as pd
from scipy import sparse
from ._types import ScissorError


def validate_problem(X, y, network, family, alpha, semantics, thresh, maxit):
    if np.iscomplexobj(X) or np.iscomplexobj(y):
        raise ScissorError("X and numeric phenotype must be real-valued")
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    if X.ndim != 2 or min(X.shape) == 0 or X.shape[0] < 2 or not np.isfinite(X).all():
        raise ScissorError("X must be finite samples × cells, with at least two samples")
    if family not in ("gaussian", "binomial", "cox"):
        raise ScissorError("family must be gaussian, binomial, or cox")
    expected = (X.shape[0], 2) if family == "cox" else (X.shape[0],)
    if y.shape != expected or not np.isfinite(y).all():
        raise ScissorError(f"{family} phenotype must be finite with shape {expected}")
    if family == "binomial" and not np.array_equal(np.unique(y), [0., 1.]):
        raise ScissorError("binomial phenotype must contain both 0 and 1")
    if family == "cox":
        if np.any(y[:, 0] < 0) or not np.isin(y[:, 1], [0, 1]).all() or not y[:, 1].sum():
            raise ScissorError("Cox requires nonnegative times, status 0/1, and at least one event")
        if semantics == "official" and X.shape[1] == 1:
            raise ScissorError("official Cox wrapper loses matrix dimensions for p=1; use corrected semantics for one cell")
    if np.ndim(alpha) != 0 or np.iscomplexobj(alpha) or not np.isfinite(alpha) or not 0 <= alpha <= 1:
        raise ScissorError("alpha must lie in [0, 1]")
    if semantics not in ("official", "corrected"):
        raise ScissorError("semantics must explicitly be 'official' or 'corrected'")
    if np.ndim(thresh) != 0 or np.iscomplexobj(thresh) or not np.isfinite(thresh) or thresh <= 0 or isinstance(maxit, (bool, np.bool_)) or not isinstance(maxit, (int, np.integer)) or maxit < 1:
        raise ScissorError("thresh must be positive and maxit a positive integer")
    W = sparse.csr_matrix(network, dtype=float, copy=True)
    if W.shape != (X.shape[1], X.shape[1]):
        raise ScissorError("network shape must match the number of cells")
    if not np.isfinite(W.data).all() or np.any(W.data < 0) or np.any(W.diagonal() != 0) or (abs(W-W.T) > 1e-12).nnz:
        raise ScissorError("normalized network must be finite, symmetric, nonnegative, with zero diagonal")
    W.sort_indices()
    if semantics == "official" and family == "binomial" and np.any(np.ptp(X, axis=0) == 0):
        raise ScissorError("official binomial wrapper drops constant columns without updating its graph; use corrected semantics or remove the column and graph node explicitly")
    return X, y, W


def align_phenotype(phenotype, sample_ids, family, positive_label=None):
    if family == "cox":
        if not isinstance(phenotype, pd.DataFrame) or set(phenotype.columns) != {"time", "status"}:
            raise ScissorError("Cox phenotype must be an ID-indexed DataFrame with exactly time and status columns")
    elif not isinstance(phenotype, pd.Series):
        raise ScissorError("phenotype must be a sample-ID-indexed pandas Series")
    if not phenotype.index.is_unique or phenotype.index.hasnans or set(phenotype.index) != set(sample_ids):
        raise ScissorError("phenotype sample IDs must equal bulk sample IDs, without duplicates or missing IDs")
    ordered = phenotype.loc[list(sample_ids)]
    if family != "binomial" and np.iscomplexobj(ordered.to_numpy()):
        raise ScissorError("numeric phenotype must be real-valued")
    mapping = None
    if family == "binomial":
        if ordered.isna().any() or ordered.nunique() != 2:
            raise ScissorError("binomial phenotype must contain exactly two nonmissing classes")
        if positive_label is None:
            if not set(ordered.unique()) == {0, 1}:
                raise ScissorError("positive_label is required unless classes are already coded 0 and 1")
            positive_label = 1
        if positive_label not in ordered.unique():
            raise ScissorError("positive_label must be one of the observed classes")
        negative = next(v for v in ordered.unique() if v != positive_label)
        mapping = {"positive": positive_label, "negative": negative}
        values = (ordered == positive_label).to_numpy(dtype=float)
    elif family == "cox":
        values = ordered[["time", "status"]].to_numpy(dtype=float)
    else:
        values = ordered.to_numpy(dtype=float)
    return values, mapping


def make_folds(y, family, foldid, nfolds, random_state):
    n = len(y)
    if foldid is None:
        if not isinstance(nfolds, (int, np.integer)) or not 2 <= nfolds <= n:
            raise ScissorError("nfolds must be an integer between 2 and the sample count")
        folds = np.arange(n) % nfolds
        np.random.default_rng(random_state).shuffle(folds)
    else:
        folds = np.asarray(foldid)
        if folds.shape != (n,) or not np.issubdtype(folds.dtype, np.integer):
            raise ScissorError("foldid must contain one integer per sample in X order")
        levels = np.unique(folds)
        if len(levels) < 2 or not np.array_equal(levels, np.arange(levels[0], levels[0]+len(levels))) or levels[0] not in (0, 1):
            raise ScissorError("foldid must use contiguous fold numbers starting at 0 or 1, with at least two folds")
        folds = folds.copy()
    for k in np.unique(folds):
        train, valid = folds != k, folds == k
        if train.sum() < 2:
            raise ScissorError(f"fold {k} has fewer than two training samples")
        if family == "binomial" and len(np.unique(y[train])) != 2:
            raise ScissorError(f"fold {k} has only one training class; provide valid foldid")
        if family == "cox" and (not y[train, 1].sum() or not y[valid, 1].sum()):
            raise ScissorError(f"fold {k} has no training or validation events; provide event-balanced foldid")
    return folds
