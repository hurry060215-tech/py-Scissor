# Full LUAD, independent R oracle

This development harness validates the full official tutorial matrices against
the R-free Python kernel. Every run uses all 4,102 cells and all 471 survival or
498 TP53 samples. The Gaussian outcome is synthetic on the full real design;
it is not a clinical endpoint or evidence of biological replication.

## Two distinct comparisons

1. `run_official_binomial.R` compiles the unchanged computational sources at
   `sunduanchen/Scissor@311560a1f160f665917e7da0fcf289914ee99773`. It runs `APML1`
   with explicit lambdas/folds, `thresh=1e-7`, `maxit=100000`, `inzero=FALSE` and
   `keep.beta=TRUE`. The original wrapper discards fold objects, so its exact
   `cvNetLogC` calls are run again to export coefficients and flags. Their
   weighted CV mean/SE must equal the unchanged wrapper's values to 1e-12.
   This comparison preserves the original clipping and inactive-entry behavior.
   Full-LUAD execution uses `-O2 -g0` without fast-math after rechecking the
   saved binomial golden fixtures; the historical small fixtures used `-O0`.
   Generated Rcpp entrypoints have external logging-only wrappers that pass
   every argument and returned object through unchanged. The first slow `-O0`
   full-data attempt was interrupted before returning any alpha result and is
   explicitly excluded from passed validation.
2. `run_corrected.R` uses a new independent mathematical optimizer in
   `corrected_oracle.cpp`, invoked from R. It is **not an upstream source patch**
   and must not be called unchanged-source compatibility. It minimizes the
   exact objectives in `docs/SEMANTICS.md` with a cell-aligned graph, true
   Gaussian curvature, exact logistic likelihood, and Breslow Cox likelihood
   divided by the original sample count. Proximal Newton subproblems use cyclic
   coordinate descent and an exact-objective Armijo line search. Cox uses a
   positive diagonal majorizer of the current sample-space Hessian. An
   independent R certificate recomputes objective and absolute KKT residual;
   its Cox implementation enumerates risk sets rather than using the C++
   cumulative-risk evaluator. Convergence requires KKT <=1e-8 unless the frozen
   manifest explicitly specifies a tighter tolerance. This differs from the
   Python accelerated proximal-gradient implementation for binary/Cox models.

The original Gaussian/Cox wrappers are never called on this fully connected
graph: their published graph-index mismatch can cause undefined memory reads.
No node is disconnected or silently removed to force an unsafe comparison.

## Frozen exchange

The validation-only prepared manifest supplies complete sample/cell/gene IDs,
outcomes, five explicit folds, three alphas `[0.05, 0.2, 0.6]`, and 20 explicit
geometric lambda values per family/alpha, from lambda_max to 0.01 lambda_max.
The numerical threshold is 1e-8 for corrected fits and 1e-7 for official fits.
Full-data, all fold paths, and cold lambda.min refits are exported separately.
Each family/alpha is evaluated even when a lower alpha already meets the
high-level selection cutoff, so the audit covers the entire frozen grid.

`exp1/X.bin`, `exp2/X.bin`, `joint.bin`, `normalized.bin`, and all exported
coefficient binaries are little-endian float64 in column-major order, with
dimensions stored in JSON. No R/Python seed equivalence is presumed.

The graph is supplied once as Matrix Market: an explicit union-symmetrized
20-nearest-neighbor graph on 20 PCs of all-gene log1p(CP10K) counts. It is not
claimed to reproduce Seurat's SNN graph. Both language paths use this same
all-cell graph; binary adjacency and symmetric degree normalization are
recomputed independently in R.

## Full preprocessing verification

`run_preprocessing.R` checks every original input against the frozen joint
matrix, using the chosen log1p(CP10K) single-cell transform and original bulk
values. Both original bulk matrices have duplicate gene names. The explicit
validation adapter uses R's `intersect` (unique names in bulk order) and
character matrix indexing (first matching source row). The strict Python
library continues to reject duplicate IDs rather than silently doing this.
The duplicate mapping is preserved in the frozen exchange.

All values in both 22,323-gene joint matrices are normalized by pinned
`preprocessCore@5f4522fcbce0820b4983f453e3d56341da62da85` (1.75.1, built with
`--disable-threading`) and compared against the Python output in bounded-size
blocks. All Pearson correlations are independently computed with R `cor()`;
blocking changes memory use, never which entries are checked. This is a full
cross-language preprocessing check, not a random spot check. It is still an
explicit-preprocessing kernel validation, not the complete Seurat workflow.

## Reproduction

Use R with Rcpp, RcppEigen, Matrix, jsonlite and the pinned preprocessCore.
Keep development dependencies under `reference/library`, or point read-only
library paths at compatible dependencies. Do not alter a user's R installation.

```sh
Rscript reference/luad/selftest.R
Rscript reference/luad/verify_official_golden.R
Rscript reference/luad/run_preprocessing.R validation-data/luad/prepared
Rscript reference/luad/run_corrected.R validation-data/luad/prepared
Rscript reference/luad/run_official_binomial.R validation-data/luad/prepared
```

`run_corrected.R` accepts an optional output directory and a third argument
`gaussian`, `binomial`, or `cox` to run one family. `run_measured.py` can wrap
each command to record actual exit code, elapsed time and child-process peak
RSS. Actual R/package versions are captured in `sessionInfo.txt`, not inferred
from earlier fixtures. Large raw-expression and quantile jobs should be
serialized with other memory-intensive validation work. Numerical design fits
need only the small bulk-by-cell matrices.

`selftest.R` uses analytical Gaussian cases and independent base-R L-BFGS-B
split-variable optimization for all three families. These small mathematical
checks test the oracle itself and are never presented as full-LUAD evidence.

## Interrupted-run recovery

`resume_corrected.R` is separate from the frozen running driver. Its default
`validate` mode checks already completed full/fold/refit files without fitting
or changing their model outputs. It verifies the input manifest and design/graph
SHA-256 values, coefficient shapes and byte counts, all convergence flags, and
every independently recomputed objective and absolute KKT residual. Audit
sidecars bind exact float64 alpha/lambda/threshold values, family, training-row
indices, input hashes and frozen optimizer-source hashes. Older checkpoints
without parameter sidecars are enrolled only after this full mathematical
reverification; this is explicitly recorded, not inferred from filenames.

```sh
Rscript reference/luad/resume_corrected.R cox \
  validation-data/luad/prepared validation-data/luad/r-oracle/corrected validate
```

After confirming that the original writer has stopped, replace `validate` with
`resume`, then pass the exact expected manifest SHA-256 and the literal argument
`confirmed-stopped`. Completed verified paths are reused; absent paths are
recomputed in full. A partial binary without its completed checkpoint JSON is
unfinished work. Contradictory or corrupt completed checkpoints cause an error,
never silent acceptance. New binaries, checkpoint JSON and final result JSON
are published by atomic renames. No failed-convergence path is accepted.
Do not run recovery concurrently with the original writer. A tighter tolerance
requires its separate explicitly derived input manifest and separate outputs;
default-threshold results cannot be relabeled as a tighter run.

`selftest_resume.R` simulates interruption after a full path and one fold,
leaves a partial second-fold binary, and checks recovery against an uninterrupted
small run. It also tests refusal of missing convergence flags, wrong shapes,
wrong input hashes, altered threshold signatures and corrupted coefficients.
This recovery test is distinct from full-LUAD numerical comparisons.

## License and provenance

The Scissor authors' unmodified source and notices remain in `../upstream`,
with provenance and license details in `../README.md` and `../../NOTICE`.
The new harness/independent optimizer are GPL-3.0-or-later, like this prototype;
see `../../LICENSE`. preprocessCore is an external LGPL development dependency,
not part of the Python runtime. Its checkout and compiler caches are ignored
and must not be accidentally vendored as new project sources.
