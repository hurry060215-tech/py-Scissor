"""ID-aware Scissor preprocessing, with an explicitly supplied sparse cell graph.

Algorithm based on sunduanchen/Scissor R/Scissor.R at REFERENCE_COMMIT.
See NOTICE for upstream authors, licensing, and modifications.
"""
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata
from ._types import PreparedData, ScissorError, REFERENCE_COMMIT

_BLOCK_COLUMNS = 64


def _validate_dense_budget(value):
    if (isinstance(value, (bool, np.bool_)) or
            not isinstance(value, (int, float, np.integer, np.floating)) or
            value <= 0 or
            (isinstance(value, (float, np.floating)) and not np.isfinite(value))):
        raise ScissorError("max_dense_bytes must be a finite positive byte budget")


def quantile_normalize(matrix: np.ndarray) -> np.ndarray:
    """Joint column quantile normalization, with average-rank interpolation.

    Rows are genes; columns are bulk samples and single cells together.
    Finite float64 data are required. Tied observations receive the target
    quantile at their average rank (linear interpolation for half ranks).
    """
    raw = np.asarray(matrix)
    if np.iscomplexobj(raw):
        raise ScissorError("expression must be real-valued; complex components are not discarded")
    a = np.asarray(raw, dtype=np.float64)
    if a.ndim != 2 or not a.size:
        raise ScissorError("expression must be a nonempty finite two-dimensional matrix")
    # Two column-block passes retain every gene/observation without allocating
    # full-sized ranks, lower/upper indices, and target-indexed copies together.
    target = np.zeros(a.shape[0])
    for start in range(0, a.shape[1], _BLOCK_COLUMNS):
        block = a[:, start:start+_BLOCK_COLUMNS]
        if not np.isfinite(block).all():
            raise ScissorError("expression must be a nonempty finite two-dimensional matrix")
        ordered = np.sort(block, axis=0)
        ordered /= a.shape[1]
        target += ordered.sum(axis=1)
    if not np.isfinite(target).all():
        raise ScissorError("quantile target exceeds finite float64 arithmetic")
    output = np.empty(a.shape, dtype=np.float64, order="F")
    for start in range(0, a.shape[1], _BLOCK_COLUMNS):
        ranks = rankdata(a[:, start:start+_BLOCK_COLUMNS], method="average", axis=0)-1
        lower = np.floor(ranks).astype(np.intp)
        upper = np.ceil(ranks).astype(np.intp)
        left, right = target[lower], target[upper]
        # Stable midpoint: a same-sign sum may overflow, while an opposite-
        # sign difference may overflow. Exact equal ranks remain unchanged.
        values = np.empty_like(left)
        same_sign = np.signbit(left) == np.signbit(right)
        values[same_sign] = left[same_sign] + (right[same_sign]-left[same_sign])/2
        values[~same_sign] = (left[~same_sign]+right[~same_sign])/2
        output[:, start:start+_BLOCK_COLUMNS] = values
    return output


def normalize_graph(graph, *, cell_ids=None, graph_cell_ids=None) -> sparse.csr_matrix:
    """Remove loops, binarize, and return D^-1/2 A D^-1/2 in CSR format.

    The graph must be undirected with nonnegative finite weights. Isolated
    nodes remain zero rows. Positional graphs require explicit graph_cell_ids
    when cell_ids are supplied, preventing accidental cell-order mismatch.
    """
    if np.iscomplexobj(graph):
        raise ScissorError("cell graph must be real-valued")
    if cell_ids is not None:
        _unique_ids(pd.Index(cell_ids), "cell_ids")
    if isinstance(graph, pd.DataFrame):
        _unique_ids(graph.index, "graph rows")
        _unique_ids(graph.columns, "graph columns")
        if set(graph.index) != set(graph.columns):
            raise ScissorError("graph row and column cell IDs must match")
        order = list(graph.index if cell_ids is None else cell_ids)
        if set(order) != set(graph.index):
            raise ScissorError("graph IDs must equal single-cell IDs")
        a = sparse.csr_matrix(graph.loc[order, order].to_numpy(dtype=float))
    else:
        if sparse.issparse(graph):
            # COO->CSR normally sums duplicates before returning. Validate and
            # binarize individual edges first so huge finite positive duplicate
            # weights cannot overflow during that implicit sum.
            edges = graph.tocoo(copy=True).astype(np.float64)
            if not np.isfinite(edges.data).all() or np.any(edges.data < 0):
                raise ScissorError("cell graph weights must be finite and nonnegative")
            edges.data[edges.data != 0] = 1.0
            a = edges.tocsr()
        else:
            a = sparse.csr_matrix(graph, dtype=np.float64, copy=True)
        if a.shape[0] != a.shape[1] or a.shape[0] == 0:
            raise ScissorError("cell graph must be a nonempty square matrix before ID alignment")
        if cell_ids is not None:
            if graph_cell_ids is None:
                raise ScissorError("graph_cell_ids are required for a positional sparse/array graph")
            ids = pd.Index(graph_cell_ids)
            _unique_ids(ids, "graph_cell_ids")
            if set(ids) != set(cell_ids) or len(ids) != a.shape[0]:
                raise ScissorError("graph IDs must equal single-cell IDs")
            permutation = ids.get_indexer(cell_ids)
            a = a[permutation, :][:, permutation].tocsr()
    if a.ndim != 2 or a.shape[0] != a.shape[1] or a.shape[0] == 0:
        raise ScissorError("cell graph must be a nonempty square matrix")
    if not np.isfinite(a.data).all() or np.any(a.data < 0):
        raise ScissorError("cell graph weights must be finite and nonnegative")
    # Binarize before canonicalizing positive duplicate edges so their sum
    # cannot overflow before the final binary adjacency is constructed.
    a.data[a.data != 0] = 1.0
    a.sum_duplicates()
    a.setdiag(0)
    a.eliminate_zeros()
    a.data[:] = 1.0
    if (a != a.T).nnz:
        raise ScissorError("cell graph must be undirected after binarization")
    degree = np.asarray(a.sum(axis=1)).ravel()
    scale = np.zeros_like(degree)
    np.divide(1.0, np.sqrt(degree), out=scale, where=degree > 0)
    network = (sparse.diags(scale) @ a @ sparse.diags(scale)).tocsr()
    network.sort_indices()
    return network


def _unique_ids(ids, name):
    if len(ids) == 0 or not ids.is_unique or ids.hasnans:
        raise ScissorError(f"{name} must be nonempty, unique and nonmissing")


def _pearson_columns(values):
    """Center/normalize without overflow or underflow of extreme profiles."""
    maxima = np.max(np.abs(values), axis=0)
    scale = np.where((maxima > 1e150) | ((maxima > 0) & (maxima < 1e-150)), maxima, 1.)
    centered = values/scale
    centered -= centered.mean(axis=0)
    norms = np.linalg.norm(centered, axis=0)
    if np.any(norms == 0) or not np.isfinite(norms).all():
        raise ScissorError("normalized expression profile has zero variance or nonfinite norm")
    centered /= norms
    return centered


def prepare(bulk: pd.DataFrame, single_cell: pd.DataFrame, graph, *,
            graph_cell_ids=None, max_dense_bytes=1_000_000_000) -> PreparedData:
    """Prepare genes × observations expression with ID-safe alignment.

    Gene intersection keeps bulk order. No log transform, count normalization,
    imputation, HVG selection, or graph construction is performed. The memory
    bound conservatively budgets normalization and the samples × cells output;
    cell graph storage remains sparse throughout.
    """
    if not isinstance(bulk, pd.DataFrame) or not isinstance(single_cell, pd.DataFrame):
        raise ScissorError("bulk and single_cell must be pandas DataFrames, genes × observations")
    _validate_dense_budget(max_dense_bytes)
    for name, frame in (("bulk", bulk), ("single_cell", single_cell)):
        _unique_ids(frame.index, name + " gene IDs")
        _unique_ids(frame.columns, name + " observation IDs")
    genes = bulk.index[bulk.index.isin(single_cell.index)]
    if len(genes) < 2:
        raise ScissorError("at least two shared genes are required")
    n, p, g = bulk.shape[1], single_cell.shape[1], len(genes)
    network = normalize_graph(graph, cell_ids=single_cell.columns,
                              graph_cell_ids=graph_cell_ids)
    width = min(n+p, _BLOCK_COLUMNS)
    # Estimate incremental workspace, not the caller-owned input DataFrames:
    # joint + normalized + bounded rank temporaries, or normalized + centered
    # bulk + bounded centered cells + similarity. 25% headroom covers indexing.
    normalization_peak = 2*g*(n+p)+10*g*width
    correlation_peak = g*(n+p)+2*g*n+3*g*width+2*n*p
    estimate = int(8*1.25*max(normalization_peak, correlation_peak))
    if estimate > max_dense_bytes:
        raise ScissorError(f"estimated dense preprocessing workspace {estimate} bytes exceeds max_dense_bytes={max_dense_bytes}")
    joint = np.empty((g, n+p), dtype=np.float64, order="F")
    for frame, offset in ((bulk, 0), (single_cell, n)):
        positions = frame.index.get_indexer(genes)
        for start in range(0, frame.shape[1], _BLOCK_COLUMNS):
            stop = min(start+_BLOCK_COLUMNS, frame.shape[1])
            block = frame.iloc[positions, start:stop].to_numpy()
            if np.iscomplexobj(block):
                raise ScissorError("expression must be real-valued; complex components are not discarded")
            block = np.asarray(block, dtype=float)
            if not np.isfinite(block).all():
                raise ScissorError("expression must be finite; missing values are not imputed")
            joint[:, offset+start:offset+stop] = block
    # A constant raw column is not a meaningful expression profile, even if
    # normalization of ties could otherwise give it an artificial variance.
    if np.any(np.min(joint, axis=0) == np.max(joint, axis=0)):
        raise ScissorError("constant expression profile cannot define Pearson similarity")
    normalized = quantile_normalize(joint)
    del joint
    centered_bulk = _pearson_columns(normalized[:, :n])
    similarity = np.empty((n,p), dtype=np.float64)
    for start in range(0, p, _BLOCK_COLUMNS):
        stop = min(start+_BLOCK_COLUMNS, p)
        centered_cells = _pearson_columns(normalized[:, n+start:n+stop])
        similarity[:, start:stop] = centered_bulk.T@centered_cells
    np.clip(similarity, -1, 1, out=similarity)
    return PreparedData(similarity, network, tuple(genes), tuple(bulk.columns),
                        tuple(single_cell.columns), normalized[:, :n], normalized[:, n:],
                        {"reference_commit": REFERENCE_COMMIT,
                         "estimated_dense_bytes": estimate,
                         "workspace_budget_excludes_caller_inputs": True,
                         "column_block_size": _BLOCK_COLUMNS,
                         "graph_nnz": network.nnz,
                         "isolated_cells": int(np.count_nonzero(np.diff(network.indptr) == 0)),
                         "dropped_bulk_genes": len(bulk.index) - g,
                         "dropped_single_cell_genes": len(single_cell.index) - g,
                         "orientation": "genes_by_observations"})
