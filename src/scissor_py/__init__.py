"""R-free Scissor research prototype with auditable numerical semantics."""
from ._types import (ScissorError, ConvergenceError, PreparedData, PathResult,
                     CVResult, ScissorResult, REFERENCE_COMMIT)
from ._preprocessing import prepare, normalize_graph, quantile_normalize
from ._api import fit_path, cross_validate, scissor, scissor_from_anndata, DEFAULT_ALPHAS
from ._anndata import result_summary, annotate_anndata, write_h5ad_atomic, SUMMARY_SCHEMA

__version__ = "0.1.0.dev5"
__all__ = ["scissor", "scissor_from_anndata", "fit_path", "cross_validate", "prepare",
           "quantile_normalize", "normalize_graph", "ScissorError", "ConvergenceError",
           "PreparedData", "PathResult", "CVResult", "ScissorResult", "REFERENCE_COMMIT", "DEFAULT_ALPHAS",
           "result_summary", "annotate_anndata", "write_h5ad_atomic", "SUMMARY_SCHEMA"]
