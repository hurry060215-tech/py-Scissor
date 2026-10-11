"""Graph ID validation must precede every potentially lossy alignment."""
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scissor_py import normalize_graph, prepare, ScissorError


def test_rectangular_graph_cannot_be_silently_cropped_by_cell_ids():
    graph = sparse.csr_matrix([[0,1,1],[1,0,1]])
    with pytest.raises(ScissorError,match="square"):
        normalize_graph(graph,cell_ids=["a","b"],graph_cell_ids=["a","b"])


def test_duplicate_output_ids_cannot_clone_graph_nodes():
    graph = sparse.csr_matrix([[0,1],[1,0]])
    with pytest.raises(ScissorError,match="unique"):
        normalize_graph(graph,cell_ids=["a","a","b"],graph_cell_ids=["a","b"])


def test_huge_duplicate_positive_edges_are_binary_without_overflow():
    graph = sparse.csr_matrix((np.array([1e308,1e308,1e308,1e308]),
                               np.array([1,1,0,0]),np.array([0,2,4])),shape=(2,2))
    with np.errstate(over="raise"):
        out = normalize_graph(graph)
    np.testing.assert_array_equal(out.toarray(),[[0,1],[1,0]])


def test_invalid_graph_fails_before_expensive_expression_normalization(monkeypatch):
    bulk=pd.DataFrame([[1,2],[2,1]],columns=["s0","s1"])
    cells=pd.DataFrame([[1,3],[3,1]],columns=["c0","c1"])
    def forbid(*args,**kwargs):
        raise AssertionError("normalization was reached before graph rejection")
    monkeypatch.setattr("scissor_py._preprocessing.quantile_normalize",forbid)
    with pytest.raises(ScissorError,match="square"):
        prepare(bulk,cells,sparse.csr_matrix((2,3)),graph_cell_ids=cells.columns)
