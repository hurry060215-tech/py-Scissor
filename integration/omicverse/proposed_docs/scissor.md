# Scissor phenotype-associated cell selection

`ov.single.Scissor` is a thin optional adapter around the independent,
GPL-3.0-or-later `scissor_py` backend. It adds no solver or graph construction to
OmicVerse. Source and current research-validation status:
https://github.com/hurry060215-tech/py-Scissor

The backend is not yet published on PyPI. Build/install its supplied source or
wheel in the same Python environment. The unrelated `pyscissor` package is not a
substitute. OmicVerse can import without the backend; `run()` checks it lazily.
This contribution is a draft while remaining high-precision real-data oracle
checks finish; it is not a claim of complete official Scissor equivalence.

```python
model = ov.single.Scissor(
    adata, bulk=bulk_expression, phenotype=phenotype,
    family="binomial", layer="lognorm", graph_key="connectivities",
    positive_label="mutant", semantics="corrected",
)
result = model.run(alpha=None, cutoff=.2, random_state=123, foldid=foldid)
```

- Bulk: genes × samples DataFrame, unique gene/sample IDs.
- Gaussian: numeric sample-ID-indexed Series. Binomial: exactly two classes;
  strings require positive_label, numeric 0/1 defaults to positive 1.
- Cox: sample-ID-indexed DataFrame with exactly time/status; 1 event, 0 censored,
  Breslow ties. Positive coefficients imply higher log hazard.
- `layer=None` explicitly selects X; otherwise provide a chosen existing layer.
  The literal string X is reserved by the backend and is rejected to avoid
  confusing X with an actual layer named X. No fallback or automatic transform.
- `graph_key` must already exist in obsp. Preprocess/build neighbors explicitly.
  A Scanpy graph need not match the original Seurat SNN graph.
- The backend aligns sample/gene/cell IDs, jointly quantile-normalizes shared
  genes and calculates sample × cell correlations. It does not drop samples,
  choose HVGs, reconcile duplicate genes or perform batch correction silently.

`semantics='official'` is the default and preserves defined published behavior,
including documented numerical quirks. It refuses potentially out-of-bounds
Gaussian/Cox graph indexing in the fixed upstream reference.
`semantics='corrected'` is an explicit alternative with aligned graph indexing,
true Gaussian curvature, convergence checks and consistent Cox lambda paths.
Read the backend's semantics/validation report before choosing between them.

Each successful alpha uses fixed-graph CV lambda.min followed by a cold refit.
Selection stops at the first alpha with fraction strictly below cutoff. An unmet
cutoff is reported; no coefficient truncation conceals it. Gaussian CV uses MSE,
Binomial deviance, and Cox event-weighted full-minus-training likelihood, not
C-index. Pass actual foldid for cross-language reproducibility; equal seeds do
not guarantee equal folds. Preserve all convergence diagnostics and label
near-zero coefficient instability separately from stable selection.

After complete successful validation/serialization, default outputs are:

- obs scissor: categorical Scissor+/Scissor-/background.
- obs scissor_coef: signed coefficients aligned by exact cell ID.
- uns scissor: params, encoding, sharedgenes, source, CV, convergence and versions.
  Encoding/sharedgenes/CV/convergence are JSON strings to preserve mixed primitive
  types through H5AD; use json.loads. A complete result_json audit is retained.

Use constructor key_added for multiple result namespaces. Existing output keys
fail unless overwrite=True is explicit. Views/backed objects need an explicit
in-memory copy before annotation. A failed run preserves old annotations and
model.result_. The returned core ScissorResult retains Background as its zero
label; AnnData uses lowercase background, with the mapping recorded.

The optional Numba engine is explicit and currently supports corrected Gaussian
only. Default runtime requires no Numba, R, rpy2 or subprocess. Use the backend's
separate write_h5ad_atomic helper for an explicit no-half-file save operation.
Association labels are not causal, clinical or diagnostic conclusions.

Run the optional tests with the backend installed:
`python -m pytest tests/single/test_scissor.py`. The scientific kernel tests and full
LUAD oracle comparisons live in the backend project, not this thin wrapper.
