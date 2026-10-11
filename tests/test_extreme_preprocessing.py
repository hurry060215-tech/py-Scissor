import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from scissor_py import normalize_graph, prepare, quantile_normalize, ScissorError


@pytest.mark.parametrize("factor", [1e308, 1e-300])
def test_extreme_finite_quantiles_and_pearson(factor):
    values = np.array([[.1, .7, .9, .1], [.8, .4, .6, .7],
                       [.5, .9, .3, .6], [1., .2, .1, .9]])
    expected = quantile_normalize(values)
    with np.errstate(over="raise", invalid="raise"):
        actual = quantile_normalize(values*factor)
    assert np.isfinite(actual).all()
    np.testing.assert_allclose(actual/factor, expected, atol=1e-15)
    bulk = pd.DataFrame(values[:, :2]*factor, index=list("abcd"), columns=["s1", "s2"])
    cells = pd.DataFrame(values[:, 2:]*factor, index=list("abcd"), columns=["c1", "c2"])
    result = prepare(bulk, cells, sparse.csr_matrix((2, 2)), graph_cell_ids=cells.columns)
    np.testing.assert_allclose(result.similarity, np.corrcoef(expected.T)[:2, 2:], atol=1e-14)


def test_extreme_same_and_opposite_sign_midpoints():
    data = np.array([[-1.6e308, -1.4e308], [0., -.5e308],
                     [0., .5e308], [1.6e308, 1.4e308]])
    with np.errstate(over="raise", invalid="raise"):
        actual = quantile_normalize(data)
    np.testing.assert_allclose(actual/1e308, [[-1.5, -1.5], [0, -.25], [0, .25], [1.5, 1.5]])


def test_coo_huge_duplicate_edges_are_binary_before_sum():
    graph = sparse.coo_matrix(([1e308]*4, ([0, 0, 1, 1], [1, 1, 0, 0])), shape=(2, 2))
    with np.errstate(over="raise", invalid="raise"):
        actual = normalize_graph(graph)
    np.testing.assert_array_equal(actual.toarray(), [[0., 1.], [1., 0.]])
    assert np.all(graph.data == 1e308)


def test_negative_duplicate_edge_cannot_be_hidden_by_coalescing():
    graph = sparse.coo_matrix(([2., -1., 1.], ([0, 0, 1], [1, 1, 0])), shape=(2, 2))
    with pytest.raises(ValueError, match="nonnegative"):
        normalize_graph(graph)
@pytest.mark.parametrize("representation", ["dense", "sparse", "dataframe"])
def test_complex_graph_is_not_silently_projected_to_real(representation):
    graph = np.array([[0, 1+2j], [1-2j, 0]])
    if representation == "sparse":
        graph = sparse.csr_matrix(graph)
    elif representation == "dataframe":
        graph = pd.DataFrame(graph)
    with pytest.raises(ScissorError, match="real-valued"):
        normalize_graph(graph)


def test_complex_expression_is_not_silently_projected_to_real():
    values = np.array([[1, 2], [3, 4]], dtype=complex)
    values[0, 0] += 1j
    with pytest.raises(ScissorError, match="real-valued"):
        quantile_normalize(values)
    with pytest.raises(ScissorError, match="real-valued"):
        prepare(pd.DataFrame(values), pd.DataFrame(values.real), sparse.csr_matrix((2, 2)), graph_cell_ids=[0, 1])


def test_quantile_target_float64_overflow_is_an_explicit_error():
    # Dividing before summation greatly extends the usable range but does not
    # promise support for every extreme float64 matrix. Never return silent inf.
    with np.errstate(over='ignore'):
        with pytest.raises(ScissorError, match='finite float64 arithmetic'):
            quantile_normalize(np.full((3, 3), np.finfo(float).max))
