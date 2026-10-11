"""Tests against actual installed, minimally patched OmicVerse, with no network."""
import copy
import json
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from offline_guard import install
NETWORK_ATTEMPTS = install()
import omicverse as ov
from scissor_py import scissor_from_anndata, ScissorError, write_h5ad_atomic
from test_anndata_adapter import inputs


def setup_model(family='binomial', layer='counts', **kwargs):
    data, bulk, y = inputs()
    if family == 'gaussian': y = pd.Series(np.sin(np.arange(18)), index=bulk.columns)
    elif family == 'cox': y = pd.DataFrame({'time':np.ceil(np.arange(1,19)/2), 'status':np.ones(18)}, index=bulk.columns)
    model = ov.single.Scissor(data, bulk=bulk, phenotype=y, family=family, layer=layer,
                             graph_key='network', semantics='corrected', **kwargs)
    return data, bulk, y, model


@pytest.mark.parametrize('family', ['gaussian', 'binomial', 'cox'])
@pytest.mark.parametrize('layer', [None, 'counts'])
def test_actual_public_class_equals_backend_and_writes_schema(family, layer):
    data, bulk, y, model = setup_model(family, layer)
    original = data.copy()
    expected = scissor_from_anndata(original, bulk, y, family=family, layer='X' if layer is None else layer,
                                   graph_key='network', semantics='corrected', alpha=.4,
                                   lambdas=[1.], foldid=np.arange(18)%3)
    result = model.run(alpha=.4, lambdas=[1.], foldid=np.arange(18)%3)
    assert model.result_ is result
    pd.testing.assert_series_equal(result.coefficients, expected.coefficients)
    pd.testing.assert_series_equal(result.labels, expected.labels)
    assert data.obs['scissor'].cat.categories.tolist() == ['Scissor+', 'Scissor-', 'background']
    assert 'scissor_coef' in data.obs
    assert set(['params','encoding','sharedgenes','source','CV','convergence','versions']).issubset(data.uns['scissor'])
    assert json.loads(data.uns['scissor']['sharedgenes']) == list(result.prepared.gene_ids)
    assert json.loads(data.uns['scissor']['CV'])[0]['foldid']['values'] == (np.arange(18)%3).tolist()
    assert (data.X != original.X).nnz == 0
    assert (data.layers['counts'] != original.layers['counts']).nnz == 0
    assert (data.obsp['network'] != original.obsp['network']).nnz == 0
    pd.testing.assert_series_equal(data.obs['original'], original.obs['original'])
    assert data.uns['original'] == original.uns['original']


def test_string_encoding_reordered_ids_and_h5ad(tmp_path):
    data, bulk, y = inputs()
    labels = y.map({0:'wildtype', 1:'mutant'}).iloc[::-1]
    folds = pd.Series(np.arange(18)%3, index=bulk.columns).iloc[::-1]
    model = ov.single.Scissor(data, bulk=bulk, phenotype=labels, family='binomial',
                             layer='counts', graph_key='network', positive_label='mutant',
                             semantics='corrected', key_added='TP53')
    result = model.run(alpha=.4, lambdas=[1], foldid=folds)
    assert result.metadata['positive_class'] == {'positive':'mutant', 'negative':'wildtype'}
    assert json.loads(data.uns['TP53']['encoding']) == result.metadata['positive_class']
    assert 'TP53_coef' in data.obs
    path = tmp_path/'result.h5ad'; write_h5ad_atomic(data, path)
    import anndata
    restored = anndata.read_h5ad(path)
    pd.testing.assert_frame_equal(restored.obs, data.obs)
    assert restored.uns['TP53'] == data.uns['TP53']


def test_collision_and_failed_second_run_preserve_first_result():
    data, _, _, model = setup_model()
    first = model.run(alpha=.4, lambdas=[1], foldid=np.arange(18)%3)
    old_obs, old_uns = data.obs.copy(deep=True), copy.deepcopy(data.uns)
    with pytest.raises(ScissorError, match='already exists'):
        model.run(alpha=.4, lambdas=[1], foldid=np.arange(18)%3)
    assert model.result_ is first
    pd.testing.assert_frame_equal(data.obs, old_obs); assert data.uns == old_uns
    model.overwrite = True
    assert model.run(alpha=.4, lambdas=[1], foldid=np.arange(18)%3) is model.result_


@pytest.mark.parametrize('failure', ['graph','layer','ambiguous_X','classes','folds','duplicate_ids','budget','maxit'])
def test_invalid_inputs_do_not_mutate_object_or_result(failure):
    data, bulk, y, model = setup_model()
    kwargs = {}
    if failure == 'graph': model.graph_key = 'absent'
    elif failure == 'layer': model.layer = 'absent'
    elif failure == 'ambiguous_X':
        data.layers['X'] = data.layers['counts'].copy(); model.layer = 'X'
    elif failure == 'classes': model.phenotype = y.map({0:'WT',1:'MUT'})
    elif failure == 'folds': kwargs['foldid'] = np.zeros(18, dtype=int)
    elif failure == 'duplicate_ids': data.obs_names = ['same']*6
    elif failure == 'budget': kwargs['max_dense_bytes'] = np.nan
    else: kwargs['maxit'] = 0
    old_obs, old_uns = data.obs.copy(deep=True), copy.deepcopy(data.uns)
    args = dict(alpha=.4, lambdas=[1], foldid=np.arange(18)%3);args.update(kwargs)
    with pytest.raises((ValueError, ScissorError)): model.run(**args)
    assert model.result_ is None
    pd.testing.assert_frame_equal(data.obs, old_obs);assert data.uns == old_uns


def test_module_adds_no_extra_backend_import_at_module_scope():
    import ast
    source = Path(importlib.util.find_spec('omicverse.single._scissor').origin).read_text()
    tree = ast.parse(source)
    imports = [node for node in tree.body if isinstance(node, (ast.Import, ast.ImportFrom))]
    assert len(imports) == 1 and imports[0].module == '_registry'
    assert 'Scissor' in ov.single.__all__
    from omicverse._registry import _global_registry
    assert any(value['full_name'] == 'omicverse.single._scissor.Scissor' for value in _global_registry._registry.values())


def test_actual_installed_packages_not_source_stubs():
    import scissor_py
    assert 'site-packages' in str(Path(ov.__file__))
    assert 'site-packages' in str(Path(scissor_py.__file__))
    assert ov.__version__ == '2.3.2'
