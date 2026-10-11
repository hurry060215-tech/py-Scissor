# API guide

Import from `scissor_py`. The core has no R or AnnData dependency. These APIs
cover APML1-Net selection, not every function in the original Scissor package.

## High-level expression workflow

```python
from scissor_py import scissor

result = scissor(
    bulk, single_cell, phenotype, graph,
    family="binomial", graph_cell_ids=cell_ids,
    positive_label="mutant", semantics="corrected",
    alpha=[.05, .2, .6], cutoff=.2,
    foldid=foldid, nlambda=100, thresh=1e-10,
)
```

This fragment uses caller-provided data. For a complete runnable example, use
[`examples/quickstart.py`](../examples/quickstart.py).

### Data contract

- `bulk` and `single_cell`: finite real DataFrames, genes in rows and observations
  in columns. Gene/sample/cell IDs must be unique and nonmissing. Common genes
  retain bulk order; zero shared genes or undefined Pearson columns fail.
- `phenotype`: Series indexed by exactly the bulk sample IDs, or Cox DataFrame
  with `time` and `status`. Order is aligned by ID. No sample is silently dropped.
- `graph`: raw nonnegative finite cell adjacency, not an already normalized W.
  It must be square and undirected after binarization. Self-loops are removed;
  isolated cells remain. Sparse/array graphs require `graph_cell_ids` at the
  high level; a square DataFrame supplies its own row/column cell IDs.
- Binary strings require `positive_label`. Numeric 0/1 uses 1 as positive unless
  explicitly configured otherwise. Class mapping is saved in metadata.
- Cox times are nonnegative and status is 0/1 (1 means event). Tied censor/event
  times share the Breslow risk set. At least one event is required, and every
  held-out fold and training set must satisfy the family validation rules.
- No implicit log transform, imputation, HVG selection, batch correction, graph
  creation or graph symmetrization is performed.

### Main controls

| Parameter | Default | Meaning |
|---|---|---|
| `family` | `"gaussian"` | `gaussian`, `binomial`, or `cox` |
| `semantics` | `"official"` | Explicit compatibility behavior; corrected mode is opt-in |
| `alpha` | `[.005,.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9]` | Scalar or ordered candidates in [0,1]; controls L1 versus graph penalty |
| `cutoff` | `.2` | Stop at first successful alpha with selected fraction **strictly less** than cutoff |
| `foldid` | `None` | Actual folds; high-level Series aligns by sample ID |
| `nfolds` | up to 10 | Used only when generating folds |
| `random_state` | `123` | NumPy fold RNG; equal R/Python seeds do not imply equal folds |
| `lambdas` | `None` | Optional explicit nonnegative descending path |
| `nlambda` | `100` | Automatic-path length |
| `rlambda` | shape-dependent | Final/initial automatic lambda ratio; 1e-4 if samples > cells, else 1e-2 |
| `thresh` | `1e-7` | Solver stopping setting; 1e-10 was explicitly used in tightened full-data validation |
| `maxit` | `100000` | Iteration cap, with family-specific diagnostics |
| `engine` | `"numpy"` | `"numba"` is optional and only supports corrected Gaussian |
| `max_dense_bytes` | `1000000000` | Dense-work estimate, not an OS/process memory limit |

Changing alpha, the graph, folds, lambda path or stopping settings changes the
experiment. Report them with the results. Review [SEMANTICS.md](SEMANTICS.md)
before treating the two modes as interchangeable.

### Selection and failure behavior

Every winning lambda must converge on the full data and every fold. Its final
single-lambda cold refit must also converge. CV uses held-out MSE for Gaussian,
binomial deviance for Binomial, and event-weighted full-minus-training Breslow
partial likelihood for Cox, not C-index. The high-level workflow chooses
`lambda_min`, not `lambda_1se`, for refitting.

Unconverged alpha candidates are recorded and skipped. If all fail,
`ConvergenceError` is raised. If no successful alpha reaches the cutoff, the last
successful result is returned with `cutoff_reached=False` and a warning.
`ScissorError` covers invalid inputs and unsupported official-compatibility
conditions. No tolerance-based relabeling hides near-zero coefficient differences.

## Results and intermediate arrays

| Record | Important fields |
|---|---|
| `PreparedData` | `similarity` (samples × cells), `network` (sparse cells × cells), IDs, normalized bulk/sc matrices, metadata |
| `PathResult` | `coefficients` (cells × lambdas), intercepts, reported lambdas, convergence flags, iterations, diagnostics |
| `CVResult` | full-data path, actual folds, mean/SE and per-fold losses, valid-candidate mask, lambda.min and lambda.1se |
| `ScissorResult` | ID-indexed coefficients/labels, selected alpha/lambda/intercept, cutoff flag, visited CV results, prepared data and metadata |

`positive_cells` and `negative_cells` are strict-sign label views. For Cox,
positive coefficients indicate higher fitted log hazard; for Binomial they
point toward the encoded positive class. None of these imply causation.

`result_summary(result)` creates a JSON-safe audit summary, not a complete
model reload format. It excludes expression and full-path coefficient matrices;
diagnostic nonfinite values are explicitly tagged. See [ANNDATA.md](ANNDATA.md).

## Low-level and preprocessing APIs

```python
import pandas as pd
from scissor_py import prepare, fit_path, cross_validate

# graph is a SciPy sparse adjacency in cell_ids order.
prepared = prepare(bulk, single_cell, graph, graph_cell_ids=cell_ids)
order = pd.Index(cell_ids).get_indexer(prepared.cell_ids)
aligned_graph = graph.tocsr()[order, :][:, order]
# numeric_y and fold_vector must follow prepared.sample_ids.
path = fit_path(prepared.similarity, numeric_y, aligned_graph, family="gaussian",
                semantics="corrected", alpha=.2, lambdas=[.1, .04])
cv = cross_validate(prepared.similarity, numeric_y, aligned_graph, family="gaussian",
                    semantics="corrected", alpha=.2, foldid=fold_vector)
```

Low-level X is samples × cells and y is positional (n-vector, or n×2 time/status
for Cox). Their graphs are raw adjacency in the same cell order as X. Reorder the caller's
raw graph as above; do not substitute the already-normalized `prepared.network`.
`quantile_normalize` exposes joint column normalization, and `normalize_graph`
exposes loop removal, binarization and symmetric degree normalization.

## AnnData and OmicVerse

`scissor_from_anndata(adata, bulk, phenotype, layer="lognorm",
graph_key="connectivities", ...)` requires both explicit keys. In this core
adapter `layer="X"` selects `adata.X`. It is nonmutating unless
`writeback=True`; writeback requires an actual in-memory, non-view AnnData.

`annotate_anndata` aligns cells by ID and defaults to collision refusal.
`write_h5ad_atomic` separately saves an H5AD with no-clobber by default.
Validation/serialization finish before annotation; failed assignments roll back.
The [full contract](ANNDATA.md) documents JSON fields and durability limits.

The proposed OmicVerse class uses `layer=None` for X and reserves the string
`"X"`. It is a separate thin adapter, not automatically installed by this
package. See its [usage and contribution status](../integration/omicverse/README.md).
