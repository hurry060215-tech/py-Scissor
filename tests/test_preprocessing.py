import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scissor_py import prepare, normalize_graph, quantile_normalize, ScissorError


def expressions():
    rng = np.random.default_rng(2)
    genes = [f"g{i}" for i in range(8)]
    bulk = pd.DataFrame(rng.normal(size=(8, 5)), index=genes, columns=[f"s{i}" for i in range(5)])
    sc = pd.DataFrame(rng.normal(size=(8, 4)), index=genes, columns=[f"c{i}" for i in range(4)])
    graph = pd.DataFrame([[9, 2, 1, 0], [2, 0, 0, 0], [1, 0, 2, 0], [0, 0, 0, 1]], index=sc.columns, columns=sc.columns)
    return bulk, sc, graph


def test_id_alignment_permutation_and_nonmutation():
    bulk, sc, graph = expressions()
    originals = [v.copy() for v in (bulk, sc, graph)]
    a = prepare(bulk, sc, graph)
    b = prepare(bulk, sc.iloc[::-1], graph.iloc[::-1, ::-1])
    np.testing.assert_allclose(a.similarity, b.similarity)
    np.testing.assert_allclose(a.network.toarray(), b.network.toarray())
    for old, new in zip(originals, (bulk, sc, graph)):
        pd.testing.assert_frame_equal(old, new)
    assert a.gene_ids == tuple(bulk.index)
    assert a.metadata["isolated_cells"] == 1


def test_graph_sparse_binarization_and_isolates():
    _, sc, graph = expressions()
    normalized = normalize_graph(graph)
    expected = np.zeros((4, 4)); expected[0, 1:3] = expected[1:3, 0] = 1/np.sqrt(2)
    assert sparse.isspmatrix_csr(normalized)
    np.testing.assert_allclose(normalized.toarray(), expected)
    reordered = normalize_graph(sparse.csr_matrix(graph.iloc[::-1, ::-1]), cell_ids=sc.columns, graph_cell_ids=sc.columns[::-1])
    np.testing.assert_allclose(reordered.toarray(), expected)
    with pytest.raises(ScissorError, match="graph_cell_ids"):
        prepare(*expressions()[:2], sparse.csr_matrix(graph))


def test_sparse_graph_large_no_dense_conversion(monkeypatch):
    n = 100_000
    graph = sparse.diags([np.ones(n-1), np.ones(n-1)], [-1, 1], format="csr")
    def fail(*args, **kwargs):
        raise AssertionError("graph was densified")
    monkeypatch.setattr(sparse.csr_matrix, "toarray", fail)
    out = normalize_graph(graph)
    assert out.nnz == 2*(n-1)
    assert out.data.nbytes+out.indices.nbytes+out.indptr.nbytes < 4_000_000


@pytest.mark.parametrize("issue", ["duplicate_gene", "duplicate_cell", "missing_gene", "nan", "constant", "asymmetric", "negative", "memory"])
def test_bad_inputs(issue):
    bulk, sc, graph = expressions()
    if issue == "duplicate_gene": bulk.index = ["x"]*len(bulk)
    if issue == "duplicate_cell": sc.columns = ["x"]*sc.shape[1]
    if issue == "missing_gene": sc.index = ["x"+s for s in sc.index]
    if issue == "nan": bulk.iloc[0, 0] = np.nan
    if issue == "constant": sc.iloc[:, 0] = 1
    if issue == "asymmetric": graph.iloc[0, 1] = 0
    if issue == "negative": graph.iloc[0, 1] = -1
    kwargs = {"max_dense_bytes": 1} if issue == "memory" else {}
    with pytest.raises(ScissorError): prepare(bulk, sc, graph, **kwargs)


def test_ties_midpoint_not_whole_block_average():
    # For a tie spanning ranks 1..3, midpoint is target[1], not mean(target[:3]).
    a = np.array([[0, 0], [0, 1], [0, 10], [8, 20]], dtype=float)
    q = quantile_normalize(a)
    target = np.sort(a, axis=0).mean(axis=1)
    assert q[0, 0] == target[1]
    assert q[0, 0] != np.mean(target[:3])
    np.testing.assert_allclose(q[:, 1], target)

