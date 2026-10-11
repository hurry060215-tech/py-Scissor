# Local OmicVerse Scissor contribution

Target: official `omicverse/omicverse` commit
`f9e1db8875d750d7f679514243a86fd2aa37d770` (version 2.3.2).
Archive URL/SHA and unchanged single/__init__.py SHA are in `source.json`.
This is a local contribution candidate until a linked upstream PR/release says
otherwise. It is not silently installed by scissor-py.

## Scope

- One new `omicverse/single/_scissor.py`, with a registered `Scissor` class.
- Two registration lines in existing `single/__init__.py`, in
  `registration.patch`. No unrelated import refactor or numerical implementation.
- `scissor_py` is imported only on `run()`. The unrelated package `pyscissor`
  does not satisfy this dependency. The missing-backend error points to the
  supplied wheel; no install/download occurs automatically.

```python
model = ov.single.Scissor(
    adata, bulk=bulk_expression, phenotype=phenotype,
    family="binomial", layer="lognorm", graph_key="connectivities",
    positive_label="mutant", semantics="corrected",
)
result = model.run(alpha=None, cutoff=.2, random_state=123, foldid=foldid)
```

`layer=None` explicitly selects X. The literal string `"X"` is rejected by this
wrapper because the underlying adapter reserves it as an X sentinel; rename an
actual layer named X explicitly. No layer/graph is guessed. `official` remains
the default semantics and its unsafe graph guards apply; use `corrected` only as
an explicit scientific choice. See the core semantics documentation.

Default outputs are obs scissor/scissor_coef and uns scissor. Background display
labels are lowercase; the returned core dataclass retains Background, with the
mapping recorded. `key_added` and `overwrite` are constructor parameters;
existing keys fail unless overwriting is explicit. `result_` changes only after
success. See `docs/ANNDATA.md` for complete serialization and file boundaries.

## Apply and verify

Run only against a private copy of the pinned source, then install that source
in an isolated environment with OmicVerse's declared dependencies:

```sh
python integration/omicverse/apply_candidate.py /path/to/omicverse-source
python -m pip install /path/to/omicverse-source
# Before installing the backend, verify actual OV imports and helpful error:
python integration/omicverse/test_missing_backend.py
python -m pip install /path/to/scissor_py-0.1.0.dev4-py3-none-any.whl
PYTHONPATH=tests:integration/omicverse python -m pytest \
  -c integration/omicverse/pytest.ini integration/omicverse/test_scissor_adapter.py
```

The application helper refuses an unexpected source hash/existing module and
is idempotent for the exact candidate. It never fetches, pushes or creates a PR.
The tests use actual installed OmicVerse and installed backend wheels, not
package stubs. They compare every family/layer to the direct backend, verify
ID/class mapping, registry entry, no-half-write errors, output collisions,
H5AD roundtrip and missing optional backend. Core tests are separate.

## Test isolation and existing upstream import side effects

Full integration tests run with a private HOME/cache and a process-local offline
guard. The guard blocks requests, DNS and INET socket traffic; it does not mock
OmicVerse modules or solver functions. The notebook permits only loopback kernel
transport. No analytics modification is part of this contribution.

At the pinned commit, `omicverse/_settings.py`, `omicverseConfig.__init__`, imports
`utils._analytics_sender` and then calls `send_analytics_full_silent`. That
function in `utils/_analytics_sender.py` reads IP geolocation and attempts to
send environment/location fields. The initial unguarded baseline attempt
failed during the line-12 import, in `utils/harness/runtime_state.py` creating a
read-only HOME directory, before that constructor's line-15 sender call. This
establishes only that call-chain boundary, not an audit of all earlier imports.
Subsequent guarded actual single-module imports succeeded; two attempted HTTP
calls were blocked. The Scissor module itself adds no network operation.

## Proposed upstream tests and guide

`proposed_tests/test_scissor.py` is self-contained for inclusion as
`tests/single/test_scissor.py` in the upstream contribution. Its ten tests passed
against the actual installed candidate and backend, with HTTP/DNS/socket traffic
blocked before imports and a private HOME/cache. The numerical tests skip when
the optional backend is absent; import/error-boundary tests remain runnable.
`proposed_docs/scissor.md` is the corresponding focused usage guide. These are
local contribution files until a draft PR is actually created.

The same ten candidate tests also pass from the patched upstream source tree
with its real shared tests/conftest.py and pytest configuration (19.16 seconds).
Only the focused Scissor file was selected; this is not a claim that the entire
OmicVerse suite was run.
