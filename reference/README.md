# Fixed-source R reference and numerical audit

This directory contains development-only reference material for a Python Scissor prototype. It does not require installing Scissor or Seurat and does not establish end-to-end equivalence with the high-level `Scissor()` workflow.

## Provenance

- Scissor: Duanchen Sun and contributors, repository `sunduanchen/Scissor`, commit `311560a1f160f665917e7da0fcf289914ee99773`.
- Exact computational sources are under `upstream/`. Text line endings were normalized to LF when saving; computational source statements were not patched.
- The original LICENSE, DESCRIPTION and README are retained. DESCRIPTION says `GPL (>= 2)` while the LICENSE text and README identify GPLv3. Source translation and redistribution must preserve applicable license obligations and attribution; academic citation alone is insufficient.
- preprocessCore: Ben Bolstad and contributors, repository `bmbolstad/preprocessCore`, commit `5f4522fcbce0820b4983f453e3d56341da62da85`, version 1.75.1. Its source was used in the isolated development environment but is deliberately omitted from the source distribution; the pinned retrieval step is below. This is a development reference version, not a claim about the user's original R installation.
- Actual execution used R 4.5.0, Rcpp, RcppEigen, Matrix and jsonlite. See `sessionInfo.txt` and `manifest.json`.

## Regeneration

Use an R installation with Rcpp, RcppEigen, Matrix and jsonlite, and a working C++ compiler. All generated compilation products remain below this directory. `Makevars` uses `-O0 -g0` to keep the development build small.

1. `Rscript reference/generate_fixtures.R`
2. `Rscript reference/audit_gaussian_cv.R`
3. To reproduce preprocessing, explicitly obtain the official pinned development dependency (no library call downloads it):
   ```sh
   git clone --no-checkout https://github.com/bmbolstad/preprocessCore.git reference/preprocessCore-source
   git -C reference/preprocessCore-source checkout 5f4522fcbce0820b4983f453e3d56341da62da85
   mkdir -p reference/library
   R CMD INSTALL --library=reference/library --configure-args=--disable-threading reference/preprocessCore-source
   Rscript reference/generate_preprocessing.R
   ```
4. `Rscript reference/generate_end_to_end.R` creates automatic-path and expression-to-selection harness fixtures after preprocessing is available.

The scripts derive the reference root from their own path. No shared R library is modified. The package source's `configure` and `cleanup` scripts may need their executable bits restored after transport.

All four regeneration scripts were rerun from `/tmp`, outside the project directory. The 22 fixture files reproduced byte-for-byte against the recorded SHA-256 manifest; the nine saved upstream files also retained their recorded hashes.

## Fixture scope

JSON matrices are row-major nested arrays. Each regression fixture includes numeric design matrix, outcome, graph, alpha, explicit lambda values and fold IDs, plus reference coefficients, CV scores and flags.

- `*_default_explicit_cv.json`: unchanged APML1 Net computation with explicit lambda/fold IDs and the normal `thresh=1e-7`. The unused L0-trimming step is explicitly disabled (`inzero=FALSE`), and the full coefficient path is retained.
- `*_tight_explicit_cv.json`: same inputs with `thresh=1e-12`. This only changes the public APML1 threshold; it does not alter hard-coded internal tolerances.
- `*_default_with_trimming.json`: separately runs the original default `inzero=TRUE`; all three supplied cases completed successfully.
- `gaussian_default_automatic_cv.json` and `binomial_default_automatic_cv.json`: automatic 20-point lambda paths with the same fixed design, graph and folds.
- `*_expression_harness.json`: raw bulk/sc expression and IDs, deliberately reordered single-cell genes and phenotype samples, one platform-specific gene in each matrix, aligned joint expression, actual normalization/similarity, all alpha CV/refit history and final labels. The harness uses explicit graph/folds with alpha `[0.05, 0.2, 0.6]`, cutoff 0.4 and 20 lambda values. It preserves high-level selection semantics but bypasses Seurat. Cox reaches the cutoff with 2/6 cells; Gaussian/Binomial finish their grids with 6/6 cells and report the cutoff unmet. All full-model and selected-refit flags are zero in these three saved cases; individual CV-fold flags are not exposed by APML1's returned object.
- `*_graph_alignment_audit.json`: unchanged direct C++ kernels are given either the official R wrapper's padded graph or a correctly aligned graph. The aligned call is a diagnostic comparison, not an unchanged high-level official oracle.
- `graph_index_audit.json`: exports graph preparation and neighbor indices, including an unsafe graph. The unsafe C++ call is never executed.
- `gaussian_curvature_audit.json` and `gaussian_cv_curvature_audit.json`: one-predictor problems with an independently known closed-form coefficient.
- `cox_lambda_audit.json`: compares internally generated, explicitly scaled and unscaled lambda paths.
- `preprocessing_ties.json`: joint expression, rank targets, actual preprocessCore output and Pearson matrix. It distinguishes midpoint-rank interpolation from averaging an entire tied block.
- `preprocessing_constant_column.json`: documents that normalization can leave a column constant and Pearson correlation undefined. The Python API should reject this situation explicitly.

The main six-cell graph deliberately isolates the last two cells so the official padded graph and Gaussian inactive-set offset remain within allocated bounds. It must not be generalized as evidence that arbitrary graphs are safe.

## Findings established by execution

1. **Graph alignment matters.** On the same safe irregular graph, changing only direct-kernel graph alignment changed coefficients by up to 0.1392131873 (Gaussian) and 0.0761088051 (Cox). Official R wrappers add a zero intercept graph node although these two kernels retain only p coefficients. For a graph with the last cell connected, R preparation produced a neighbor location 7 when p was 6: `beta[loc-1]` would index 6 outside a six-entry zero-based vector. Gaussian also has an inactive-set access using `beta[loc]`. These out-of-bounds calls were diagnosed from actual prepared indices, never executed or emulated.

2. **Gaussian default convergence does not guarantee the requested coefficient tolerance.** A known optimum is 1.328072177342894; the unchanged default kernel returned 1.3197105661060884 with a success flag. Tightening both direct-kernel tolerances produced 1.3280716360654257. Its reported RSS was 3.8802974133 while recomputing RSS from the returned coefficient gave 0.8782219018, confirming inconsistent RSS bookkeeping for non-unit predictor variance.

3. **Gaussian CV has a hidden stopping threshold.** A training-fold optimum is 1.281498588245183. The default CV kernel returned 1.2729607317733347. Tightening only the exposed APML1 threshold to 1e-12 did not change it, because `thresh2=1e-5` remains fixed in the R wrapper. Directly tightening both thresholds gave 1.2814980646191785. A compatibility implementation should not silently replace this behavior while claiming exact unchanged-default parity.

4. **Cox automatic path scaling differs from wrapper labels/CV input.** In the supplied case, lambda_max is 0.2368093421544794. The automatic C++ path and explicitly scaled path were exactly equal; comparing against unscaled ratios changed coefficients by up to 0.3393124548. The R wrapper still reported and used original ratios `[1, 0.2, 0.04, 0.008]` for its CV. Explicit numeric lambda fixtures avoid ambiguity and should precede automatic-path compatibility claims.

5. **Quantile tie handling is midpoint-rank interpolation.** Actual preprocessCore 1.75.1 output matched interpolation at the average tied rank within 1.78e-15. Averaging target values over the whole tied block differed by 0.2. This distinction is material for zero-heavy expression matrices.

## Interpretation and limits

These fixtures establish selected finite, reproducible upstream behaviors. They do not prove a mathematically corrected solver is equivalent to the original implementation, and undefined memory access must never be treated as a compatibility requirement. Corrected behavior should carry an explicit label and separate mathematical validation. A success flag alone is not a KKT or objective-optimality certificate.

No performance claim, LUAD biological replication, full Seurat workflow, or validation of the user's separate R 4.5.2 environment is implied. The reference installation and compilation are development dependencies only; the Python runtime must not invoke R.
