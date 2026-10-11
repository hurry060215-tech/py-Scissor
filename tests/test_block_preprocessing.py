"""Full-data blocking preserves the reference definition; it is not sampling."""
import tracemalloc
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.stats import rankdata
from scissor_py import quantile_normalize, prepare


def dense_reference(a):
    target=np.sort(a,axis=0).mean(axis=1)
    ranks=rankdata(a,method="average",axis=0)-1
    return (target[np.floor(ranks).astype(int)]+target[np.ceil(ranks).astype(int)])/2


def test_quantile_multiple_blocks_and_zero_heavy_ties():
    rng=np.random.default_rng(328)
    a=rng.integers(0,4,size=(131,259)).astype(float)
    a[rng.random(a.shape)<.75]=0
    original=a.copy()
    expected=dense_reference(a)
    np.testing.assert_allclose(quantile_normalize(a),expected,atol=1e-12,rtol=1e-12)
    np.testing.assert_array_equal(a,original)


def test_full_preprocessing_multiple_blocks_matches_dense_reference():
    rng=np.random.default_rng(22)
    bulk=pd.DataFrame(rng.normal(size=(97,73)),index=[f"g{i}" for i in range(97)],columns=[f"s{i}" for i in range(73)])
    cells=pd.DataFrame(rng.normal(size=(97,137)),index=bulk.index,columns=[f"c{i}" for i in range(137)])
    graph=sparse.csr_matrix((137,137))
    result=prepare(bulk,cells.iloc[::-1],graph,graph_cell_ids=cells.columns)
    normalized=dense_reference(np.column_stack([bulk,cells]))
    centered=normalized-normalized.mean(axis=0)
    centered/=np.linalg.norm(centered,axis=0)
    np.testing.assert_allclose(result.normalized_bulk,normalized[:,:73],atol=1e-12)
    np.testing.assert_allclose(result.normalized_single_cell,normalized[:,73:],atol=1e-12)
    np.testing.assert_allclose(result.similarity,centered[:,:73].T@centered[:,73:],atol=1e-12)
    assert result.gene_ids==tuple(bulk.index)
    assert result.cell_ids==tuple(cells.columns)


def test_quantile_peak_memory_is_bounded_beyond_output():
    rng=np.random.default_rng(102)
    a=rng.normal(size=(1000,1024))
    tracemalloc.start()
    out=quantile_normalize(a)
    _,peak=tracemalloc.get_traced_memory()
    tracemalloc.stop()
    # Output is 8.2 MB; ranks/indices are bounded to 64 columns. The former
    # implementation required >6 whole matrix temporaries on this case.
    assert out.shape==a.shape
    assert peak < a.nbytes*3, peak
@pytest.mark.parametrize("budget", [np.nan, np.inf, -np.inf, 0, -1, True, None, "1 GB"])
def test_invalid_memory_budget_is_rejected_before_allocation(budget):
    from scissor_py import prepare, ScissorError
    bulk = pd.DataFrame([[1., 2.], [3., 4.]])
    cells = pd.DataFrame([[2., 3.], [4., 5.]])
    with pytest.raises(ScissorError, match="finite positive byte budget"):
        prepare(bulk, cells, sparse.csr_matrix((2, 2)), graph_cell_ids=[0, 1], max_dense_bytes=budget)
