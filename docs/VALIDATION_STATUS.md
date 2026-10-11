# Validation evidence and current scope

Historical numerical validation used 0.1.0.dev4; current packaging is 0.1.0.dev5. **Declared full-data numerical grid completed.**
No claim of complete equivalence with all official Scissor inputs or the original
Seurat workflow. GPL-3.0-or-later; R is a development reference only.

## Historical local checks completed (2026-10-10)

- 387 core/preprocessing/AnnData/packaging tests passed under Python 3.12.14, NumPy 2.3.5,
  SciPy 1.17.0 and pandas 2.2.3. Ordinary tests use saved synthetic R fixtures;
  no R process or data download is required. Installed-wheel testing outside
  source passes 382 tests with five explicit optional-Numba skips.
- A direct build from the development checkout passes the explicit distribution
  member and notebook-content audit. Saved notebook warning paths are normalized
  with explicit metadata, and private execution prefixes fail the gate. Source archives allowlist public documents and exclude real
  cell-level diagnostics, including when stale build manifests exist. Wheels
  contain only the Python package and distribution metadata. This gate also
  runs after builds in the proposed CI workflow.
- 18 tests against actual installed OmicVerse 2.3.2 pinned source and the installed
  backend pass; a further self-contained upstream test file has ten passing
  tests under the same offline guard. Missing scissor_py does not break OV import; Scissor.run raises
  an explicit optional-dependency error. No OV package stubs were used.
- H5AD roundtrip (AnnData 0.12.2 / h5py 3.16.0), ID/class types, sparse graph,
  layers, annotations, collisions, failure rollback and atomic file replacement
  tested. JSON summaries retain explicit nonfinite-diagnostic tags.
- Seven exact notebook code cells and all assertions pass in a shared Python
  namespace. This executor rejects both TCP and Unix IPC binding, so Jupyter
  kernel execution is **not verified**. The notebook records the actual method.
- Complete NumPy versus explicit Numba Gaussian 20-lambda/five-fold/refit
  comparison passes, max floating difference 5.06e-15 with identical iterations,
  flags and chosen lambdas. Observed total elapsed 2677.0 versus 321.2 seconds,
  including 4.19 seconds cold JIT; concurrent CPU loads mean this is not a general
  speed guarantee. Default import/runtime requires neither Numba nor R.

## Full public LUAD experiment

All original 4,102 single cells, 471 survival samples and 498 TP53 samples are
retained. Both cohorts have 22,323 common genes. Explicit validation-only
first-row duplicate-gene handling matches R character indexing and is disclosed;
the normal API still rejects duplicates. The graph is an all-gene PCA20/kNN20
union graph, **not** historical Seurat SNN. Same frozen graph, actual five folds,
alpha [.05,.2,.6], and twenty-lambda paths are shared between Python and R.
Twenty lambdas/five folds are a declared bounded audit, not the 100-lambda default.

All 204,768,879 joint/normalized entries were independently checked against
actual R/preprocessCore; normalized max error <=1.01e-11, Pearson <=4.78e-15.
The current blocked implementation and full public AnnData pipelines for all
three families also pass, with reversed gene/sample/fold-ID orders. Full
pipeline peak RSS was about 2.2–2.8 GB; preparing all large reference intermediates
peaked at 4.56 GB and wrote 4.45 GB. These are measured dataset-specific values.

- **Unchanged official Binomial:** all three alpha full paths, every fold and cold
  refit pass. Coefficients agree around 1e-14, CV curves/choices and stable signs
  agree. This is the official stopping/clipping behavior, not an exact-objective
  optimality guarantee.
- **Official Gaussian/Cox on the full graph:** correctly rejected because the
  published wrapper's graph indexing is unsafe. No undefined R access was run.
  Safe synthetic fixtures separately validate defined official behavior.
- **Corrected baseline validation threshold 1e-8 (the API default is 1e-7):** alpha .05 passes full path/CV/refit in all
  three families. Higher-alpha Gaussian, Binomial and Cox have coefficient-tolerance
  failures even though stable signs and chosen lambdas agree. Gaussian max
  absolute discrepancy reaches 1.42e-4; Binomial reaches 7.82e-5, and Cox
  reaches 1.45e-4. All baseline R families have completed; Cox exits cleanly. Failures remain
  failures; no tolerance has been relaxed.
- **Additional both-sided threshold 1e-10:** Gaussian all three alphas pass full/
  fold/refit coefficient tests (maximum differences 1.39e-7, 5.65e-7 and 1.40e-6);
  the independent R driver exits cleanly. Binomial all three alphas also pass, with
  maximum coefficient differences 9.57e-8, 2.13e-7 and 9.93e-7; its independent R
  driver exits cleanly. Cox alpha .05 also passes the entire path/five-fold/refit
  comparison (maximum coefficient difference 9.37e-8). Cox alpha .2 passes too
  (maximum difference 3.00e-7), including every fold and refit.
  Cox alpha .6 also passes (maximum difference 1.4131e-6). All nine tightened
  family/alpha full paths, every CV fold and cold refit meet the original
  coefficient, prediction, objective, CV and sign checks. Earlier Python
  tightening alone was insufficient at high alpha, hence both sides were rerun.

Initial acceptance limits: normalization/Pearson atol1e-8, rtol1e-6;
coefficients/intercepts/predictions atol1e-5, rtol1e-4; objective atol1e-8,
rtol1e-6; CV loss/SE atol1e-6, rtol1e-4. Stable nonzero signs and CV choices are
checked separately; near-zero differences are reported, not relabelled away.
Independent corrected R optimization and certificates do not modify the fixed
original Scissor source.

Full real input matrices, phenotype/sample listings, per-cell reference outputs
and private checkpoint history are intentionally not distributed here. Explicit
download/preparation/comparison scripts and aggregate findings support independent
reproduction. Synthetic fixtures, licensed attribution notices and original source hashes are included; full upstream code is in the pinned development archive.

Machine-readable [aggregate numerical checks](luad/numerical_summary.json)
retain each failed comparison and show absent reports as not_reported, never a
pass. They omit individual sample/cell results.

## Recovery and remaining boundaries

The original tightened Cox execution channel lost contact after an executor-key
change; its lifecycle remains unknown, and no successful exit is attributed to
that original channel. Ten complete nodes were copied, hash-checked and
independently recertified into a new output root. The other eleven complete
nodes were computed there under unchanged inputs and thresholds. An old-log
hash/mtime guard prevented new scheduling if the original writer advanced. The
new isolated driver completed with exit 0 in 5,116.55 seconds, peak RSS 462 MB.
The final independent Cox recertification also exits 0 (25.00 seconds). Across
all three tightened families, 63 nodes contain 1,089 fitted lambda vectors and
4,467,078 coefficients, each independently recomputed R KKT residual <=1e-10.
No unfinished node was reused. See the [recovery protocol](https://github.com/hurry060215-tech/py-Scissor/tree/2cb65ba5af3e1ffac05226f66d0e0ddf1eb7dc50/reference/luad/ISOLATED_RECOVERY.md).

The numerical grid is complete: fifteen same-threshold comparisons pass and
six lower-precision baseline comparisons remain failures. This does not validate
every possible input or promise coefficient equivalence at the API default
threshold. The stricter 1e-10 setting above is explicit, not a changed default.
The original unsafe official Gaussian/Cox full graph, biological effectiveness,
and actual Jupyter-kernel execution are not claimed as passed.

Exact-version source, installed-wheel and archive-boundary checks accompany the
delivered artifacts, with file hashes and a local source-snapshot identifier.
Current hosted CI results are shown in [GitHub Actions](https://github.com/hurry060215-tech/py-Scissor/actions/workflows/tests.yml).
This historical local validation report does not establish hosted CI or upstream
acceptance. The OmicVerse adapter remains a contribution candidate until accepted upstream. No biological effectiveness,
causal, clinical or diagnostic validation is asserted.

## Repository refresh (2026-10-11)

The restored source passed all 387 existing tests in a fresh Python 3.12 environment.
Six repository-contract tests were then added for required fixture integrity,
version metadata, the runnable README example and explicit threshold derivation.
The resulting suite passed 393 tests. This is a packaging/documentation and
Python-regression refresh; the historical full-data R experiments above were not
rerun. Source, wheel and hosted-CI results are checked separately.

The source-outside-tree clean wheel suite passed 388 tests with five explicit
optional-Numba skips. All three model smoke fits and the exact README example
also passed with network and child-process creation blocked. Hosted CI for
[source commit 2cb65ba5](https://github.com/hurry060215-tech/py-Scissor/commit/2cb65ba5af3e1ffac05226f66d0e0ddf1eb7dc50)
completed successfully in all three jobs: Python 3.10/3.12 core tests, package
build/boundary checks and isolated-wheel examples, plus Python 3.12 tests with
AnnData and Numba. [Recorded run](https://github.com/hurry060215-tech/py-Scissor/actions/runs/38105094573).

## Core-only repository split (2026-10-11)

The original R/C++ references and external OmicVerse/tutorial/skill contributions
have moved out of the main tree and core packages. All Python numerical and
AnnData tests remain, with the original JSON fixtures unchanged. Historical
external-interface results above describe the archived candidate, not files
bundled with version 0.1.0.dev5. See [architecture](ARCHITECTURE.md).
