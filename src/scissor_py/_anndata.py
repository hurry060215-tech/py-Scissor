"""Explicit, optional result annotation and JSON-safe audit summaries.

AnnData is imported only when annotation or H5AD writing is requested. No
expression matrix is copied into ``uns`` and no pickle is used.
"""
import hashlib
import inspect
from importlib.metadata import version
import json
import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from ._types import ScissorError, ScissorResult

SUMMARY_SCHEMA = "scissor-py-result-summary/1"
_NONFINITE = "__scissor_nonfinite_float__"


def _json_value(value, location="metadata"):
    """Preserve ordinary JSON types; tag nonfinite diagnostics explicitly."""
    if isinstance(value, np.ndarray):
        return _json_value(value.tolist(), location)
    if isinstance(value, np.bool_):
        value = bool(value)
    elif isinstance(value, np.integer):
        value = int(value)
    elif isinstance(value, np.floating):
        value = float(value)
    elif isinstance(value, np.str_):
        value = str(value)
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if np.isfinite(value):
            return value
        return {_NONFINITE: "nan" if np.isnan(value) else ("+inf" if value > 0 else "-inf")}
    if isinstance(value, dict):
        if _NONFINITE in value:
            raise ScissorError(f"reserved nonfinite marker in {location}")
        if not all(isinstance(key, str) for key in value):
            raise ScissorError(f"JSON audit metadata requires string keys at {location}")
        return {key: _json_value(item, f"{location}.{key}") for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item, f"{location}[{i}]") for i, item in enumerate(value)]
    raise ScissorError(f"unsupported JSON audit value at {location}: {type(value).__name__}")


def _ids(values, location):
    out = []
    for value in values:
        item = _json_value(value, location)
        if item is None or not isinstance(item, (str, int, float, bool)):
            raise ScissorError(f"{location} must use finite string/number/bool scalar IDs for JSON export")
        out.append(item)
    index = pd.Index(values)
    if not index.is_unique or index.hasnans:
        raise ScissorError(f"{location} must be unique and nonmissing")
    return out


def _array(values, shape, location):
    a = np.asarray(values)
    if a.shape != shape or a.dtype.kind not in "biuf":
        raise ScissorError(f"invalid numeric array at {location}; expected shape {shape}")
    # Convert once, at the top-level boundary; re-encoding tagged diagnostics
    # would incorrectly interpret our own nonfinite markers as caller metadata.
    return {"dtype": str(a.dtype), "shape": list(a.shape), "values": a}


def _graph_summary(network, p):
    if not sparse.issparse(network) or network.shape != (p, p):
        raise ScissorError("result graph must be sparse and match its cell IDs")
    graph = network.tocsr(copy=True)
    graph.sum_duplicates(); graph.eliminate_zeros(); graph.sort_indices()
    if np.iscomplexobj(graph) or not np.isfinite(graph.data).all():
        raise ScissorError("result graph must have finite real values")
    digest = hashlib.sha256()
    for value, dtype in ((graph.shape, "<i8"), (graph.indptr, "<i8"), (graph.indices, "<i8"), (graph.data, "<f8")):
        digest.update(np.asarray(value, dtype=dtype).tobytes())
    return {"format": "csr", "shape": list(graph.shape), "nnz": graph.nnz,
            "dtype": str(graph.dtype), "index_dtype": str(graph.indices.dtype),
            "normalized_network_sha256": digest.hexdigest(),
            "hash_encoding": "canonical CSR shape/indptr/indices little-endian int64 then data float64",
            "storage": "edges are retained in caller/result graph, not duplicated in this audit summary"}


def result_summary(result):
    """Return a schema-versioned, standard-JSON-safe result audit dictionary.

    Includes final coefficients/labels, typed IDs, CV curves/folds, model and
    preprocessing provenance, and a sparse-graph fingerprint. Excludes dense
    expression/similarity and full-path coefficient matrices; it is an audit
    summary, not a complete model reload format. Nonfinite diagnostics have an
    explicit tagged representation, while selected coefficients must be finite.
    Unsupported custom ID/metadata objects raise rather than being stringified.
    """
    if not isinstance(result, ScissorResult):
        raise ScissorError("result must be a ScissorResult")
    from . import __version__
    prepared = result.prepared
    cell_ids = _ids(prepared.cell_ids, "cell_ids")
    sample_ids = _ids(prepared.sample_ids, "sample_ids")
    gene_ids = _ids(prepared.gene_ids, "gene_ids")
    p, n = len(cell_ids), len(sample_ids)
    for name, series in (("coefficients", result.coefficients), ("labels", result.labels)):
        if not isinstance(series, pd.Series) or not series.index.is_unique or set(series.index) != set(prepared.cell_ids):
            raise ScissorError(f"result {name} IDs must equal its prepared cell IDs")
    coefficients = result.coefficients.reindex(prepared.cell_ids).to_numpy()
    labels = result.labels.reindex(prepared.cell_ids).tolist()
    if coefficients.dtype.kind not in "fiu" or not np.isfinite(coefficients).all():
        raise ScissorError("selected coefficients must be finite and real")
    expected = np.where(coefficients > 0, "Scissor+", np.where(coefficients < 0, "Scissor-", "Background"))
    if labels != expected.tolist():
        raise ScissorError("result labels do not match strict coefficient signs")
    if result.family not in ("gaussian", "binomial", "cox") or not np.isfinite([result.alpha, result.lambda_min, result.intercept]).all():
        raise ScissorError("result model parameters must be finite and valid")
    if not 0 <= result.alpha <= 1 or result.lambda_min < 0:
        raise ScissorError("result alpha/lambda are outside their valid range")
    if prepared.similarity.shape != (n, p):
        raise ScissorError("result similarity shape does not match sample/cell IDs")
    cv_records = []
    for i, cv in enumerate(result.cv_results):
        count = len(cv.path.lambdas)
        if cv.path.family != result.family or cv.path.coefficients.shape != (p, count):
            raise ScissorError("CV family/coefficients do not match the result")
        folds = np.asarray(cv.foldid)
        if folds.shape != (n,) or folds.dtype.kind not in "iu":
            raise ScissorError("CV fold IDs do not match bulk samples")
        cv_records.append({"alpha": cv.path.alpha,
                           "lambdas": _array(cv.path.lambdas, (count,), f"cv[{i}].lambdas"),
                           "foldid": _array(folds, (n,), f"cv[{i}].foldid"),
                           "mean_loss": _array(cv.mean_loss, (count,), f"cv[{i}].mean_loss"),
                           "standard_error": _array(cv.standard_error, (count,), f"cv[{i}].standard_error"),
                           "fold_loss": _array(cv.fold_loss, (len(np.unique(folds)), count), f"cv[{i}].fold_loss"),
                           "valid": _array(cv.valid, (count,), f"cv[{i}].valid"),
                           "full_converged": _array(cv.path.converged, (count,), f"cv[{i}].full_converged"),
                           "full_iterations": _array(cv.path.iterations, (count,), f"cv[{i}].full_iterations"),
                           "lambda_min": cv.lambda_min, "lambda_1se": cv.lambda_1se,
                           "path_metadata": cv.path.metadata, "cv_metadata": cv.metadata})
    out = {"schema_version": SUMMARY_SCHEMA, "package_version": __version__, "family": result.family,
           "alpha": result.alpha, "lambda_min": result.lambda_min, "intercept": result.intercept,
           "cutoff_reached": bool(result.cutoff_reached),
           "ids": {"cell_ids": cell_ids, "sample_ids": sample_ids, "gene_ids": gene_ids},
           "coefficients": _array(coefficients, (p,), "coefficients"), "labels": labels,
           "network": _graph_summary(prepared.network, p), "cv_results": cv_records,
           "preprocessing": prepared.metadata, "metadata": result.metadata,
           "nonfinite_diagnostic_encoding": {_NONFINITE: "nan, +inf or -inf; never used for selected coefficients"}}
    # The descriptor above is schema documentation, not an encoded diagnostic.
    descriptor = out.pop("nonfinite_diagnostic_encoding")
    out = _json_value(out, "result")
    out["nonfinite_diagnostic_encoding"] = descriptor
    return out


def _check_writeback(adata, key_added, overwrite):
    try:
        from anndata import AnnData
    except ImportError as error:
        raise ScissorError("annotation requires the optional anndata package") from error
    if not isinstance(adata, AnnData):
        raise ScissorError("annotation requires an actual AnnData object")
    if adata.is_view or adata.isbacked:
        raise ScissorError("annotation requires an in-memory, non-view AnnData; explicitly copy/to_memory first")
    if not isinstance(key_added, str) or not key_added.strip() or key_added in (".", "..") or "/" in key_added or "\x00" in key_added:
        raise ScissorError("key_added must be a nonempty H5AD-safe string without '/' or NUL")
    if not isinstance(overwrite, (bool, np.bool_)):
        raise ScissorError("overwrite must be boolean")
    coefficient_key = key_added+"_coef"
    if not overwrite and (key_added in adata.obs or coefficient_key in adata.obs or key_added in adata.uns):
        raise ScissorError("Scissor output key already exists; choose another key_added or explicitly set overwrite=True")
    if not adata.obs_names.is_unique or adata.obs_names.hasnans:
        raise ScissorError("AnnData cell IDs must be unique and nonmissing")
    return coefficient_key


def annotate_anndata(adata, result, *, key_added="scissor", overwrite=False):
    """Write labels/coefficients to obs and the JSON audit summary to uns.

    Cells are aligned by exact ID, not position. Validation and serialization
    finish before mutation. The two in-memory assignments are rolled back if
    either fails. Views/backed objects are explicitly rejected. No file is
    written, and unrelated obs/uns/layers/graphs remain unchanged. Returns None.
    """
    coefficient_key = _check_writeback(adata, key_added, overwrite)
    summary = result_summary(result)
    if set(adata.obs_names) != set(result.prepared.cell_ids) or len(adata.obs_names) != len(result.prepared.cell_ids):
        raise ScissorError("AnnData cell IDs must exactly equal result cell IDs")
    summary["annotation"] = {"obs_label_key": key_added, "obs_coefficient_key": coefficient_key,
                              "display_label_mapping": {"Scissor+": "Scissor+", "Scissor-": "Scissor-", "Background": "background"},
                              "target_cell_order_matches_model": list(adata.obs_names) == list(result.prepared.cell_ids)}
    text = json.dumps(summary, ensure_ascii=False, allow_nan=False, sort_keys=True, separators=(",", ":"))
    old_obs, old_uns = adata.obs, adata.uns
    new_obs, new_uns = old_obs.copy(deep=True), old_uns.copy()
    display_labels = result.labels.reindex(adata.obs_names).replace({"Background": "background"})
    new_obs[key_added] = pd.Categorical(display_labels, categories=["Scissor+", "Scissor-", "background"])
    new_obs[coefficient_key] = result.coefficients.reindex(adata.obs_names).to_numpy(copy=True)
    def encoded(value):
        return json.dumps(value, ensure_ascii=False, allow_nan=False, sort_keys=True, separators=(",", ":"))
    new_uns[key_added] = {
        "schema_version": SUMMARY_SCHEMA, "result_json": text,
        "params": {"family": result.family, "alpha": float(result.alpha),
                   "lambda_min": float(result.lambda_min), "intercept": float(result.intercept),
                   "cutoff_reached": bool(result.cutoff_reached),
                   "semantics": str(result.metadata.get("semantics", "unspecified"))},
        "encoding": encoded(summary["metadata"].get("positive_class")),
        "sharedgenes": encoded(summary["ids"]["gene_ids"]),
        "source": {"reference_commit": str(result.metadata.get("reference_commit", "unspecified")),
                   "normalized_network_sha256": summary["network"]["normalized_network_sha256"]},
        "CV": encoded(summary["cv_results"]),
        "convergence": encoded({"refit": summary["metadata"].get("refit"),
                                 "failed_alpha_candidates": summary["metadata"].get("failed_alpha_candidates", [])}),
        "versions": {name: version(name) for name in ("numpy", "scipy", "pandas", "anndata")},
        "json_fields": "encoding, sharedgenes, CV, convergence, result_json",
    }
    new_uns[key_added]["versions"]["scissor_py"] = summary["package_version"]
    try:
        adata.obs = new_obs
        adata.uns = new_uns
    except BaseException:
        adata.obs = old_obs
        adata.uns = old_uns
        raise


def write_h5ad_atomic(adata, filename, *, overwrite=False):
    """Explicitly save an in-memory AnnData via a sibling temporary file.

    Successful file replacement is atomic. Serialization errors leave an
    existing destination unchanged and remove the temporary file. Requires an
    AnnData writer supporting convert_strings_to_categoricals=False so saving
    does not alter the input's metadata dtypes. This function never annotates or
    fits, and never creates parent directories automatically.
    """
    try:
        from anndata import AnnData
    except ImportError as error:
        raise ScissorError("H5AD writing requires the optional anndata package") from error
    if not isinstance(adata, AnnData) or adata.is_view or adata.isbacked:
        raise ScissorError("H5AD writing requires an in-memory, non-view AnnData")
    if not isinstance(overwrite, (bool, np.bool_)):
        raise ScissorError("overwrite must be boolean")
    target = Path(filename)
    parameters = inspect.signature(adata.write_h5ad).parameters
    if "convert_strings_to_categoricals" not in parameters and not any(
            p.kind == inspect.Parameter.VAR_KEYWORD for p in parameters.values()):
        raise ScissorError("atomic H5AD writing requires an AnnData version supporting convert_strings_to_categoricals=False")
    if target.exists() and not overwrite:
        raise FileExistsError(target)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=target.parent, prefix="."+target.name+"-", suffix=".h5ad", delete=False) as handle:
            temporary = Path(handle.name)
        adata.write_h5ad(temporary, convert_strings_to_categoricals=False)
        with temporary.open("rb") as handle:
            os.fsync(handle.fileno())
        if overwrite:
            os.replace(temporary, target)
        else:
            # Atomic no-clobber create, including a destination appearing while
            # H5AD serialization was in progress. Both paths share a filesystem.
            os.link(temporary, target)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()
