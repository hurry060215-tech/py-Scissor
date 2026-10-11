# AnnData annotation and portable audit records

`scissor_from_anndata` remains nonmutating by default. Set `writeback=True`, or
call `annotate_anndata(adata, result)` after examining a result. Annotation
requires a real in-memory, non-view AnnData. Backed inputs remain supported for
read-only fitting; explicitly call `to_memory()` before requesting annotation.

- `obs['scissor']`: categorical `Scissor+`, `Scissor-`, `background`.
- `obs['scissor_coef']`: signed final coefficients, aligned by exact cell ID.
- `uns['scissor']`: params, encoding, sharedgenes, source, CV, convergence,
  versions and a complete schema-versioned audit JSON string.
- `key_added='TP53'` changes all three names to TP53/TP53_coef/TP53.
- Existing output keys fail unless `overwrite=True` is explicit. All fitting,
  ID checks and JSON serialization finish before mutation. Assignment failure
  rolls back obs and uns. No expression, layer, graph or other metadata changes.

The returned dataclass retains the original backend label `Background`.
AnnData's lowercase display mapping is explicit in the audit record and does
not alter coefficient signs or introduce a numerical threshold.

`encoding`, `sharedgenes`, `CV`, `convergence`, and `result_json` in `uns` are
standard JSON strings: use `json.loads(adata.uns['scissor']['CV'])`. This avoids
loss of numeric versus string ID/class types through H5AD heterogeneous arrays.
`params`, `source`, and `versions` are directly inspectable dictionaries.

`result_summary(result)` returns the full JSON-safe audit dictionary. It records
final coefficients, typed IDs, all visited CV curves and actual fold IDs,
convergence/iteration metadata, reference SHA and the normalized sparse graph's
canonical SHA-256. It excludes dense normalized expression/similarity and full
path coefficients; it is not a full model-reload format. Graph edges remain in
the caller's `obsp` and in `result.prepared.network`. Nonfinite diagnostic values
are tagged explicitly as `{"__scissor_nonfinite_float__":"nan"}`, `+inf`, or
`-inf`; selected coefficients/model parameters must be finite. Unsupported
custom ID/metadata objects and non-string metadata dictionary keys raise rather
than being stringified or pickled. JSON metadata preserves primitive types;
array records include dtype and shape.

`write_h5ad_atomic(adata, filename, overwrite=False)` is an explicit separate
file operation. It writes a sibling temporary H5AD, fsyncs it, and atomically
creates the destination without clobbering a racing writer; `overwrite=True`
uses atomic replacement. Failed serialization leaves an existing destination
unchanged and removes the temporary file. Parent directories are not created.
This does not promise power-loss durability of the directory entry or atomicity
across multiple files. Annotation and disk writing are separate operations.

The writer disables AnnData's implicit string-to-category conversion to preserve
input metadata dtypes. It requires an AnnData writer with that public option;
actual integration/H5AD tests use AnnData 0.12.2 and h5py 3.16.0. Core package
import, numerical fitting and JSON audit export do not require AnnData.
