"""Actual optional AnnData/H5AD input checks, separate from future writeback."""
import copy
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scissor_py import scissor_from_anndata, ScissorError

ad=pytest.importorskip("anndata")


def inputs():
    rng=np.random.default_rng(600)
    genes=pd.Index([f"g{i}" for i in range(30)])
    cells=pd.Index([f"c{i}" for i in range(6)])
    bulk=pd.DataFrame(rng.poisson(8,size=(30,18)),index=genes,columns=[f"s{i}" for i in range(18)])
    counts=rng.poisson(7,size=(6,30)).astype(np.float64)
    data=ad.AnnData(sparse.csr_matrix(np.log1p(counts)),obs=pd.DataFrame(index=cells),var=pd.DataFrame(index=genes))
    data.layers["counts"]=sparse.csr_matrix(counts)
    data.obsp["network"]=sparse.diags([np.ones(5),np.ones(5)],[-1,1],format="csr")
    data.obs["original"]=pd.Categorical(["a","a","b","b","b","a"])
    data.uns["original"]={"source":"synthetic"}
    phenotype=pd.Series(np.arange(18)%2,index=bulk.columns)
    return data,bulk,phenotype


def run(data,bulk,y,**kwargs):
    return scissor_from_anndata(data,bulk,y,layer=kwargs.pop("layer","counts"),graph_key="network",
                                family="binomial",alpha=.4,lambdas=[1],foldid=np.arange(18)%3,**kwargs)


def test_explicit_layer_and_real_anndata_remain_unchanged():
    data,bulk,y=inputs()
    original=data.copy()
    result=run(data,bulk,y)
    assert result.metadata["anndata_layer"]=="counts"
    assert result.prepared.cell_ids==tuple(data.obs_names)
    np.testing.assert_array_equal(data.X.toarray(),original.X.toarray())
    np.testing.assert_array_equal(data.layers["counts"].toarray(),original.layers["counts"].toarray())
    pd.testing.assert_frame_equal(data.obs,original.obs)
    assert data.uns==original.uns
    assert result.prepared.normalized_single_cell.shape==(30,6)


def test_h5ad_round_trip_preserves_inputs_and_results(tmp_path):
    data,bulk,y=inputs()
    before=run(data,bulk,y)
    path=tmp_path/"input.h5ad"
    data.write_h5ad(path)
    recovered=ad.read_h5ad(path)
    after=run(recovered,bulk,y)
    pd.testing.assert_series_equal(before.coefficients,after.coefficients)
    pd.testing.assert_series_equal(before.labels,after.labels)
    np.testing.assert_allclose(before.prepared.similarity,after.prepared.similarity)
    assert recovered.uns==data.uns


@pytest.mark.parametrize("failure",["missing_layer","missing_graph","duplicate_ids","invalid_folds","memory"])
def test_failure_has_no_partial_anndata_mutation(failure):
    data,bulk,y=inputs()
    kwargs={}
    if failure=="missing_layer": kwargs["layer"]="missing"
    if failure=="missing_graph": del data.obsp["network"]
    if failure=="duplicate_ids": data.obs_names=["same"]*6
    if failure=="invalid_folds": y.iloc[:]=0
    if failure=="memory":kwargs["max_dense_bytes"]=1
    old_obs=data.obs.copy(deep=True);old_uns=copy.deepcopy(data.uns)
    old_X=data.X.copy();old_graphs=set(data.obsp)
    with pytest.raises((ScissorError,ValueError)):
        run(data,bulk,y,**kwargs)
    pd.testing.assert_frame_equal(data.obs,old_obs)
    assert data.uns==old_uns
    assert set(data.obsp)==old_graphs
    assert (data.X!=old_X).nnz==0


def test_backed_sparse_X_is_read_only_and_can_be_budgeted(tmp_path):
    data,bulk,y=inputs()
    path=tmp_path/"backed.h5ad"
    data.write_h5ad(path)
    backed=ad.read_h5ad(path,backed="r")
    try:
        result=run(backed,bulk,y,layer="X")
        expected=run(data,bulk,y,layer="X")
        pd.testing.assert_series_equal(result.coefficients,expected.coefficients)
        np.testing.assert_allclose(result.prepared.similarity,expected.prepared.similarity)
        assert backed.isbacked
    finally:backed.file.close()
@pytest.mark.parametrize("budget", [np.nan, np.inf, None, True])
def test_invalid_adapter_budget_never_materializes_backed_matrix(budget):
    data, bulk, y = inputs()
    class BackedSentinel:
        shape = (6, 30)
        def to_memory(self):
            raise AssertionError("invalid budget must be rejected before materialization")
    class AdapterInput:
        obsp = data.obsp
        layers = {"counts": BackedSentinel()}
        obs_names = data.obs_names
        var_names = data.var_names
    with pytest.raises(ScissorError, match="finite positive byte budget"):
        run(AdapterInput(), bulk, y, max_dense_bytes=budget)


@pytest.mark.parametrize('family', ['gaussian', 'cox'])
def test_complex_numeric_phenotypes_rejected_before_float_conversion(family):
    data, bulk, y = inputs()
    if family == 'gaussian':
        phenotype = y.astype(complex)+1j
    else:
        phenotype = pd.DataFrame({'time':np.arange(1,19).astype(complex)+1j, 'status':y}, index=bulk.columns)
    with pytest.raises(ScissorError, match='real-valued'):
        scissor_from_anndata(data, bulk, phenotype, layer='counts', graph_key='network',
                            family=family, semantics='corrected', alpha=.4, lambdas=[1],
                            foldid=np.arange(18)%3, writeback=True)
    assert 'scissor' not in data.obs and 'scissor' not in data.uns
