# Reproducing the validation

There are three different levels of evidence. Do not present ordinary CI as a
new full-LUAD run or a corrected oracle as the unchanged official implementation.
Measured results and retained failures are in [VALIDATION_STATUS.md](VALIDATION_STATUS.md).

## 1. Fast R-free regression tests

```sh
python -m pip install -e '.[test,integration-test,acceleration]'
python -m pytest -q
python examples/quickstart.py
python examples/synthetic.py
python -m build
python scripts/check_distribution.py dist/*
```

Twenty-two committed synthetic R fixtures and the fixed upstream source have
SHA-256 records in `reference/manifest.json`. Tests require these files and
verify their hashes. Missing required fixtures are failures, not skips.
Optional AnnData/Numba tests can skip without their extras; the optional CI job
installs both. No R or data download is used by this test suite.

## 2. Regenerate the small R references deliberately

The original runs used R 4.5.0 with Rcpp, RcppEigen, Matrix, jsonlite and pinned
preprocessCore 1.75.1. Use an isolated R development library and C++ compiler.
Follow [reference/README.md](../reference/README.md) for pinned dependency
retrieval and these scripts:

```sh
Rscript reference/generate_fixtures.R
Rscript reference/audit_gaussian_cv.R
Rscript reference/generate_preprocessing.R
Rscript reference/generate_end_to_end.R
```

Regeneration writes synthetic development outputs. Keep a copy or clean checkout
of the original fixtures before comparing hashes. Unsafe upstream graph indexing
is diagnosed from indices, never deliberately executed as an out-of-bounds call.

## 3. Full public LUAD experiment

These are explicit developer actions, not downloads triggered by an API call.
Original files are about 254 MB; preparing all frozen expression references
wrote about 4.45 GB and peaked at 4.56 GB process RSS. Reserve ample RAM/disk and
run large preprocessing tasks serially. Kernel fits retain all cells and
samples but only need the prepared bulk × cell designs.

```sh
python -m pip install -e '.[validation]'
python scripts/download_luad_validation.py validation-data/luad/original
python scripts/prepare_luad_validation.py validation-data/luad/original validation-data/luad/prepared
```

The downloader reads only named official public files, checks expected sizes
and records hashes. Inspect source changes rather than substituting data if
an endpoint fails. Preparation freezes all IDs, raw-input hashes, the graph,
actual folds, all three outcomes and alpha/lambda paths in a local manifest.
Its explicit first-row duplicate-gene policy reproduces R indexing; the normal
library still rejects duplicate IDs. See [INPUT_AUDIT.md](luad/INPUT_AUDIT.md).

With the isolated R dependencies ready, first test the independent oracle and
the unchanged official Binomial source:

```sh
Rscript reference/luad/selftest.R
Rscript reference/luad/verify_official_golden.R
Rscript reference/luad/run_preprocessing.R validation-data/luad/prepared
```

### Keep baseline and tighter runs separate

The prepared manifest uses a corrected **baseline 1e-8** threshold. The ordinary
API default remains **1e-7**. The completed tighter comparison explicitly used
**1e-10 on both sides**, with the same design, graph, folds and paths.

On a filesystem supporting symbolic links, create a new metadata-only contract:

```sh
python scripts/derive_luad_validation_contract.py validation-data/luad/prepared validation-data/luad/tight-inputs --thresh 1e-10
```

The helper leaves original files unchanged and links them read-only by
convention; it does not enforce filesystem permissions. Keep the prepared
directory present and immutable. Only stopping settings and explicit derivation
metadata change. Do not reuse baseline outputs as tighter results.

Run each family/alpha with the original input manifest and an explicit Python
override. For example:

```sh
python scripts/run_luad_python.py validation-data/luad/prepared validation-data/luad/python-tight --family gaussian --semantics corrected --alpha .05 --thresh 1e-10
Rscript reference/luad/run_corrected.R validation-data/luad/tight-inputs validation-data/luad/r-tight gaussian
python scripts/compare_luad_oracle.py validation-data/luad/prepared validation-data/luad/python-tight/gaussian_corrected_alpha_0.05 validation-data/luad/r-tight/gaussian/alpha-0.05 validation-data/luad/comparisons/gaussian-tight-0.05.json
```

The R family command evaluates all three alphas. Repeat the Python command and
comparison for `.2` and `.6`, then for `binomial` and `cox`. This runs all twenty
lambdas, five fixed folds and the cold refit, even where high-level alpha search
would stop earlier. Corrected Gaussian may explicitly use `--engine numba` when
its optional extra is installed; it never falls back silently.

For the **baseline**, use the original prepared directory on both sides and omit
the Python threshold override. For **unchanged official Binomial**, use
`--semantics official` and `reference/luad/run_official_binomial.R`; official
1e-7 stopping remains distinct from exact KKT optimization. Do not run official
Gaussian/Cox against this unsafe connected graph.

`compare_luad_oracle.py` checks full/fold/refit coefficients, intercepts,
predictions, objectives, CV losses/SE, convergence, chosen lambdas and signs at
the original tolerances. It exits nonzero on a failure and preserves the report.
`summarize_luad_validation.py` exports aggregate metrics without individual IDs.
The [R oracle guide](../reference/luad/README.md) describes certificates and
interrupted-run validation. Never resume over a directory whose writer may
still be running. The historical [isolated recovery](../reference/luad/ISOLATED_RECOVERY.md)
is evidence about that run, not a generic restart command.

No original LUAD matrices, individual phenotypes, per-cell real outputs or private
checkpoint history are distributed here. Reproduction requires deliberate
access to public sources and sufficient local resources. No long R job is rerun
automatically by installation, CI or ordinary examples.
