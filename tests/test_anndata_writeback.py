"""Real optional AnnData annotation, audit JSON, and atomic file boundaries."""
import copy
import json
import os
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scissor_py import (ScissorError, SUMMARY_SCHEMA, annotate_anndata,
                        result_summary, write_h5ad_atomic)
from test_anndata_adapter import inputs, run

ad = pytest.importorskip('anndata')

@pytest.fixture
def fitted():
    data, bulk, y = inputs()
    return data, bulk, y, run(data, bulk, y)


def test_standard_json_and_complete_audit_contract(fitted):
    data, bulk, y, result = fitted
    summary = json.loads(json.dumps(result_summary(result), allow_nan=False))
    assert summary['schema_version'] == SUMMARY_SCHEMA
    assert summary['ids']['cell_ids'] == list(data.obs_names)
    assert summary['ids']['sample_ids'] == list(bulk.columns)
    assert summary['coefficients']['dtype'] == 'float64'
    assert summary['coefficients']['shape'] == [data.n_obs]
    assert summary['metadata']['positive_class'] == {'positive': 1, 'negative': 0}
    assert summary['metadata']['requested_alphas'] == [.4]
    assert summary['metadata']['fold_source'] == 'provided'
    assert summary['metadata']['lambda_mode'] == 'explicit'
    assert summary['network']['nnz'] == data.obsp['network'].nnz
    assert len(summary['network']['normalized_network_sha256']) == 64
    cv = summary['cv_results'][0]
    assert cv['foldid']['values'] == (np.arange(18) % 3).tolist()
    assert cv['fold_loss']['shape'] == [3, 1]
    assert cv['valid']['dtype'] == 'bool'
    assert 'similarity' not in summary and 'normalized_bulk' not in summary


def test_tagged_nonfinite_diagnostics_do_not_emit_nonstandard_json(fitted):
    result = fitted[-1]
    result.cv_results[0].mean_loss[0] = np.nan
    result.cv_results[0].standard_error[0] = np.inf
    result.metadata['diagnostics'] = [np.float64(-np.inf), np.int32(7), np.bool_(True)]
    summary = json.loads(json.dumps(result_summary(result), allow_nan=False))
    marker = '__scissor_nonfinite_float__'
    assert summary['cv_results'][0]['mean_loss']['values'] == [{marker: 'nan'}]
    assert summary['cv_results'][0]['standard_error']['values'] == [{marker: '+inf'}]
    assert summary['metadata']['diagnostics'] == [{marker: '-inf'}, 7, True]


@pytest.mark.parametrize('failure', ['custom', 'key', 'marker', 'coefficient', 'labels', 'ids', 'cv_shape'])
def test_invalid_summary_never_partially_annotates(fitted, failure):
    data, bulk, y, result = fitted
    before_obs, before_uns = data.obs.copy(deep=True), copy.deepcopy(data.uns)
    if failure == 'custom': result.metadata['bad'] = object()
    elif failure == 'key': result.metadata[123] = 'bad'
    elif failure == 'marker': result.metadata['bad'] = {'__scissor_nonfinite_float__': 'nan'}
    elif failure == 'coefficient': result.coefficients.iloc[0] = np.nan
    elif failure == 'labels': result.labels.iloc[0] = 'invalid'
    elif failure == 'ids': result.coefficients.index = ['x']*len(result.coefficients)
    else: result.cv_results[0].fold_loss = np.zeros((2, 1))
    with pytest.raises(ScissorError): annotate_anndata(data, result)
    pd.testing.assert_frame_equal(data.obs, before_obs)
    assert data.uns == before_uns


def test_numeric_ids_and_class_types_are_preserved(fitted):
    result = fitted[-1]
    result.prepared.cell_ids = tuple(range(10, 16))
    result.prepared.sample_ids = tuple(range(100, 118))
    result.prepared.gene_ids = tuple(range(30))
    result.coefficients.index = result.prepared.cell_ids
    result.labels.index = result.prepared.cell_ids
    summary = json.loads(json.dumps(result_summary(result), allow_nan=False))
    assert summary['ids']['cell_ids'] == list(range(10, 16))
    assert isinstance(summary['metadata']['positive_class']['positive'], int)


def test_annotation_aligns_cells_without_touching_inputs(fitted):
    data, bulk, y, result = fitted
    data = data[::-1].copy()
    original = data.copy()
    annotate_anndata(data, result, key_added='phenotype_cells')
    np.testing.assert_array_equal(data.obs['phenotype_cells_coef'], result.coefficients.loc[data.obs_names])
    assert data.obs['phenotype_cells'].tolist() == result.labels.loc[data.obs_names].replace({'Background': 'background'}).tolist()
    assert data.obs['phenotype_cells'].cat.categories.tolist() == ['Scissor+', 'Scissor-', 'background']
    assert data.uns['original'] == original.uns['original']
    pd.testing.assert_series_equal(data.obs['original'], original.obs['original'])
    assert (data.X != original.X).nnz == 0
    assert (data.layers['counts'] != original.layers['counts']).nnz == 0
    assert (data.obsp['network'] != original.obsp['network']).nnz == 0
    summary = json.loads(data.uns['phenotype_cells']['result_json'])
    assert summary['annotation']['target_cell_order_matches_model'] is False


@pytest.mark.parametrize('collision', ['label', 'coefficient', 'uns'])
def test_collision_requires_explicit_overwrite(fitted, collision):
    data, bulk, y, result = fitted
    if collision == 'uns': data.uns['scissor'] = {'old': 7}
    else: data.obs['scissor' if collision == 'label' else 'scissor_coef'] = 'old'
    original = data.copy()
    with pytest.raises(ScissorError, match='already exists'): annotate_anndata(data, result)
    pd.testing.assert_frame_equal(data.obs, original.obs)
    assert data.uns == original.uns
    annotate_anndata(data, result, overwrite=True)
    assert data.uns['scissor']['schema_version'] == SUMMARY_SCHEMA


@pytest.mark.parametrize('key', ['', ' ', '.', '..', 'a/b', 'a\x00b', None, 12])
def test_invalid_namespace_rejected_before_fit(fitted, key, monkeypatch):
    data, bulk, y, _ = fitted
    import scissor_py._api as api
    monkeypatch.setattr(api, 'scissor', lambda *args, **kwargs: pytest.fail('fit must not start'))
    with pytest.raises(ScissorError, match='key_added'): run(data, bulk, y, writeback=True, key_added=key)


@pytest.mark.parametrize('value', [None, 0, 'yes'])
def test_writeback_boolean_required(fitted, value):
    with pytest.raises(ScissorError, match='boolean'): run(*fitted[:3], writeback=value)


def test_view_backed_and_mismatched_cells_rejected(fitted, tmp_path):
    data, bulk, y, result = fitted
    with pytest.raises(ScissorError, match='non-view'): annotate_anndata(data[:], result)
    different = data.copy(); different.obs_names = ['different'] + list(data.obs_names[1:])
    with pytest.raises(ScissorError, match='exactly equal'): annotate_anndata(different, result)
    path = tmp_path/'backed.h5ad'; data.write_h5ad(path)
    backed = ad.read_h5ad(path, backed='r')
    try:
        with pytest.raises(ScissorError, match='non-view'): annotate_anndata(backed, result)
        with pytest.raises(ScissorError, match='non-view'): run(backed, bulk, y, writeback=True)
    finally: backed.file.close()


def test_setter_failure_rolls_back_obs_and_uns(fitted, monkeypatch):
    data, _, _, result = fitted
    old_obs, old_uns = data.obs.copy(deep=True), copy.deepcopy(data.uns)
    prop = ad.AnnData.uns
    calls = []
    def fail_once(self, value):
        calls.append(1)
        if len(calls) == 1: raise RuntimeError('injected uns setter failure')
        return prop.fset(self, value)
    monkeypatch.setattr(ad.AnnData, 'uns', property(prop.fget, fail_once, prop.fdel))
    with pytest.raises(RuntimeError, match='injected'): annotate_anndata(data, result)
    pd.testing.assert_frame_equal(data.obs, old_obs)
    assert data.uns == old_uns


def test_fit_failure_has_no_half_writeback(fitted):
    data, bulk, y, _ = fitted
    old_obs, old_uns = data.obs.copy(deep=True), copy.deepcopy(data.uns)
    with pytest.raises(ScissorError): run(data, bulk, y*0, writeback=True)
    pd.testing.assert_frame_equal(data.obs, old_obs)
    assert data.uns == old_uns


def test_explicit_writeback_and_atomic_h5ad_roundtrip(fitted, tmp_path):
    data, bulk, y, _ = fitted
    # Object string columns are intentionally not converted to categoricals.
    data.obs['free_text'] = ['x','x','y','y','y','x']
    data.var['free_text'] = ['gene']*data.n_vars
    result = run(data, bulk, y.map({0:'WT', 1:'mutant'}), positive_label='mutant', writeback=True)
    before = data.copy()
    path = tmp_path/'result.h5ad'
    write_h5ad_atomic(data, path)
    restored = ad.read_h5ad(path)
    pd.testing.assert_frame_equal(restored.obs, data.obs)
    pd.testing.assert_frame_equal(data.obs, before.obs)
    pd.testing.assert_frame_equal(data.var, before.var)
    assert restored.uns['scissor'] == data.uns['scissor']
    assert json.loads(restored.uns['scissor']['result_json'])['metadata']['positive_class'] == {'positive':'mutant', 'negative':'WT'}
    assert restored.obsp['network'].format == data.obsp['network'].format
    assert restored.obsp['network'].dtype == data.obsp['network'].dtype
    assert restored.obsp['network'].indices.dtype == data.obsp['network'].indices.dtype
    assert (restored.obsp['network'] != data.obsp['network']).nnz == 0
    assert (restored.layers['counts'] != data.layers['counts']).nnz == 0
    assert (restored.X != data.X).nnz == 0
    assert path.is_file() and list(tmp_path.iterdir()) == [path]
    assert result.prepared.cell_ids == tuple(restored.obs_names)
    previous = path.read_bytes()
    with pytest.raises(FileExistsError): write_h5ad_atomic(data, path)
    assert path.read_bytes() == previous
    write_h5ad_atomic(data, path, overwrite=True)
    assert list(tmp_path.iterdir()) == [path]


def test_atomic_write_failure_preserves_old_file_and_metadata(fitted, tmp_path):
    data = fitted[0]
    data.obs['plain_strings'] = ['a']*len(data)
    before_obs = data.obs.copy(deep=True)
    path = tmp_path/'old.h5ad'; path.write_bytes(b'prior content')
    data.uns['unsupported'] = object()
    with pytest.raises(Exception): write_h5ad_atomic(data, path, overwrite=True)
    assert path.read_bytes() == b'prior content'
    assert list(tmp_path.iterdir()) == [path]
    pd.testing.assert_frame_equal(data.obs, before_obs)


def test_atomic_no_clobber_race(fitted, tmp_path, monkeypatch):
    data = fitted[0]
    path = tmp_path/'race.h5ad'
    old_writer = ad.AnnData.write_h5ad
    def writer(self, filename, **kwargs):
        old_writer(self, filename, **kwargs)
        path.write_bytes(b'appeared during serialization')
    monkeypatch.setattr(ad.AnnData, 'write_h5ad', writer)
    with pytest.raises(FileExistsError): write_h5ad_atomic(data, path)
    assert path.read_bytes() == b'appeared during serialization'
    assert list(tmp_path.iterdir()) == [path]


def test_atomic_writer_does_not_create_directories_or_accept_views(fitted, tmp_path):
    data = fitted[0]
    with pytest.raises(FileNotFoundError): write_h5ad_atomic(data, tmp_path/'missing'/'x.h5ad')
    assert not (tmp_path/'missing').exists()
    with pytest.raises(ScissorError): write_h5ad_atomic(data[:], tmp_path/'view.h5ad')
    with pytest.raises(ScissorError): write_h5ad_atomic(data, tmp_path/'bad.h5ad', overwrite='yes')


def test_legacy_writer_is_rejected_before_creating_any_file(fitted, tmp_path, monkeypatch):
    def legacy_writer(self, filename=None, compression=None):
        pytest.fail('unsupported writer must never be invoked')
    monkeypatch.setattr(ad.AnnData, 'write_h5ad', legacy_writer)
    with pytest.raises(ScissorError, match='AnnData version supporting'):
        write_h5ad_atomic(fitted[0], tmp_path/'unsupported.h5ad')
    assert not list(tmp_path.iterdir())
