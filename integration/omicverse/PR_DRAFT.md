Title: Draft: add an optional ov.single.Scissor adapter

## Summary

- Add one thin `ov.single.Scissor` class and two registration/export lines.
- Load the independent `scissor_py` backend only when `run()` is called. Missing
  backend support does not change normal OmicVerse import behavior.
- Require explicit bulk phenotypes, cell graph and expression-layer choice;
  write ID-aligned annotations only after fitting and serialization succeed.
- Add a focused usage guide and ten self-contained adapter tests. No solver,
  preprocessing, graph creation, dependency auto-install, or data download is
  added to OmicVerse.

The separate GPL-3.0-or-later backend is at
https://github.com/hurry060215-tech/py-Scissor and is not published on PyPI. The
unrelated `pyscissor` package does not provide this backend. Install its supplied
source/wheel into the same environment to run the optional numerical tests.

## API

```python
model = ov.single.Scissor(
    adata, bulk=bulk_expression, phenotype=phenotype,
    family="binomial", layer="lognorm", graph_key="connectivities",
    positive_label="mutant", semantics="corrected",
)
result = model.run(alpha=None, cutoff=0.2, random_state=123, foldid=foldid)
```

`layer=None` explicitly selects X. Graph construction remains the caller's
choice. Default obs keys are scissor/scissor_coef and uns key is scissor;
constructor `key_added` supports multiple results. Collisions fail unless
`overwrite=True` is explicit. A failed run preserves prior annotations and
`model.result_`. H5AD-safe audit JSON records encoding, shared genes, actual folds,
CV/convergence, graph fingerprint and package versions.

## Scientific scope and draft status

The fixed official Scissor source has reproduced graph-index and numerical
quirks. The backend keeps explicit `official` and `corrected` semantics; it
rejects potentially unsafe official Gaussian/Cox graph indexing rather than
silently changing the model. It does not claim complete original Seurat/Scissor
workflow equivalence. See the backend's current validation and semantics docs.

Full public LUAD preprocessing, all-cell AnnData pipelines, unchanged official
Binomial paths, and both-sided high-precision Gaussian/Binomial comparisons have
passed. High-precision Cox comparisons have also passed all three alphas,
including every fold and cold refit after independently verified checkpoint
recovery. Some
baseline-threshold coefficient comparisons failed and remain recorded. This
draft is reviewable interface code. The declared numerical grid is complete;
it does not establish biological effectiveness or universal official equivalence. No raw data, individual phenotype/cell outputs or clinical
effectiveness claims are included in this contribution.

## Local validation

- Actual installed OmicVerse at f9e1db8875d750d7f679514243a86fd2aa37d770, with
  only the proposed adapter/registration change; no OmicVerse module stubs.
- Eighteen integration tests and ten self-contained upstream-candidate tests
  pass with the actual installed backend, private HOME/cache and HTTP/DNS/socket
  traffic blocked before imports.
- Missing backend import/run boundaries, all three families, X/named layers,
  reversed IDs and string class encoding, collisions/failure preservation,
  registry lookup and actual H5AD roundtrip are covered.
- The separate backend has 387 passing local tests, including packaging boundaries; wheel tests outside source
  pass (382 plus five optional Numba skips), and all three runtime families work
  with child processes forbidden. R is development-only.

This change deliberately leaves existing OmicVerse package-import design and
analytics code untouched. The local test network guard isolates those existing
side effects; it does not replace any OmicVerse numerical implementation.
