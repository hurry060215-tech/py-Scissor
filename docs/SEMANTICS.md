# Numerical contract and compatibility boundaries

## Scope

This is the phenotype-guided APML1 **Net** selection kernel, not a full Seurat,
Omicverse, APML1-penalty, or Scissor UI reimplementation. The caller supplies an
explicit undirected graph and a chosen expression layer. Default alpha order is
0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9;
default cutoff is 0.2 and high-level nlambda is 100. Folds are reused across
alphas. NumPy seed 123 is deterministic within this implementation; it does not
reproduce R's random fold allocation. Explicit foldid is the cross-language
contract. The cutoff comparison is strictly less than, as in Scissor.R.

## Graph and memory

Nonzero nonnegative edge weights become one, loops are removed, and undirected
adjacency A becomes W = D^(-1/2) A D^(-1/2). Isolated rows stay zero. The penalty
uses I-W, so isolated cells still receive the diagonal quadratic penalty.
Directed graphs are rejected rather than silently symmetrized. Cell ordering is
explicit for sparse/array graphs and label aligned for DataFrames.

Graph storage and operations stay sparse. Centered X, expression normalization,
and bulk-by-cell similarity are dense. `max_dense_bytes` estimates preprocessing
workspace and rejects before materialization; it is not a hard process-memory
cap. Path coefficients use O(cells × lambdas) storage. Cox uses O(samples) risk
arrays, without an N×N risk-set matrix.

## Published Gaussian/Cox graph-index boundary

The fixed R wrappers prepend a zero graph node. Binomial also prepends an
intercept coefficient and is aligned. Gaussian and Cox keep p coefficients,
causing the published network to refer to shifted cells and sometimes beyond
the vector. Gaussian additionally uses loc without subtracting one in its
inactive-set test.

`reference/fixtures/graph_index_audit.json` constructs an actual p=6 graph in R:
the unsafe case produces loc=7, which becomes beta[6] in zero-based C++, beyond
the last valid beta[5]. The dangerous kernel access was **not executed**. This is
an R-produced graph/index bounds proof, not a claimed observed crash. Safe
irregular graph inputs were executed: padded vs aligned graph inputs changed
Gaussian coefficients by up to 0.139213 and Cox by up to 0.076109 in the saved
audit fixtures.

For 0 ≤ alpha < 1, official Gaussian conservatively requires the final two
original graph nodes to be isolated; official Cox requires the final node to
be isolated. This ensures every potential published offset stays in bounds.
It then preserves the padded/shifted semantics on the admitted subset. The
bound can reject cases where a particular active-set trajectory would happen
not to touch an unsafe node. It never guesses the contents of out-of-bounds
memory or silently removes offending edges. Alpha=1 follows the unnetworked
Lasso branch; Gaussian p=1 is also a safe upstream special case.
The official Cox wrapper fails for p=1 because its R row subset drops the
matrix dimension. Official mode rejects that input; corrected mode supports it.

Corrected mode uses the graph in true cell order for all families. This is an
explicit model change relative to the padded published Gaussian/Cox wrapper.
It is not presented as full high-level official equivalence.

## Objectives and stops

Let beta be cell coefficients, a the unpenalized intercept (Gaussian/binomial),
and n the original number of bulk samples. Corrected mode minimizes:

- Gaussian: ||y-a-X beta||²/(2n) + lambda alpha ||beta||₁ + lambda(1-alpha) betaᵀ(I-W)beta/2
- Binomial: mean(log(1+exp(a+X beta))-y(a+X beta)) plus the same penalty
- Cox: -Breslow partial log likelihood/n plus the same penalty; no intercept

Gaussian uses actual centered coordinate second moments. Binomial and Cox use
stable likelihood calculations, sparse graph products, monotone accelerated
proximal-gradient steps, backtracking, and absolute KKT tolerance `thresh`.
Binomial fitting in corrected mode uses the exact logistic objective without
probability clipping. Corrected Cox retains the original n divisor after
removing uninformative records censored before the first event; tie/event order
matches the Breslow convention. No baseline hazard is estimated.

Official mode translates the published coordinate order, warm starts,
active-set rules, probability clipping, and stops. Gaussian retains the unit
coordinate denominator and RSS bookkeeping even though X is only centered.
The recurrence is not the true RSS; its stopping condition can end before the
desired coefficient tolerance. Tightening public APML1 thresh alone leaves its
hardcoded secondary threshold 1e-5 unchanged. Binomial retains the original
signed network term outside abs in the inactive-set admission test. Cox retains
ifast behavior and returned-path truncation; unavailable columns are excluded.
Official convergence means these published stops, and metadata includes KKT
diagnostics to show the distinction from mathematical convergence.

The Python maximum-iteration bound is enforced even when newly admitted
coordinates would cause the original local-restart loop to exceed it. Such a
candidate is marked failed. Candidate acceptance also rejects any failed fold
or full-data fit, whereas the R wrappers do not consistently use every flag.
These are explicit safety deviations from upstream failure handling.

## Lambda and CV

Automatic ratios use rlambda=0.0001 when n>p, else 0.01; geometric path length
is explicit. Lambda maximum uses the null score divided by n and alpha (0.001
when alpha=0), following the reference. A zero score yields a zero path.

Official Cox computes the full-data fit with scaled lambdas but exposes the
unscaled ratio vector through the R wrapper, which CV and lambda.min refit
then use. This quirk is reproduced deliberately; PathResult.lambdas contains
the wrapper labels and metadata.solver_lambdas contains actual full-fit values.
For explicit lambdas these agree. Corrected mode always reports and reuses
actual solver lambdas. There is no silent rescaling after selection.

CV statistics are fold-size-weighted held-out MSE for Gaussian and clipped
binomial deviance (probabilities limited to [1e-5, 1-1e-5]) for binomial.
Cox uses -(LL_full(beta_train)-LL_train(beta_train))/held_out_event_count,
weighted by held-out event counts. This is the published full-minus-training
partial-likelihood criterion, not C-index or standalone held-out likelihood.
The standard error is the square root of weighted squared fold deviations
divided by number_of_folds-1. The first minimum wins ties, and lambda.1se is
the first valid path entry within one standard error.

All folds must have at least two training observations; binary training folds
must contain both classes, and Cox training/validation folds both need events.
Invalid folds are rejected, not omitted. Constant binomial X columns are
rejected in official mode because the R wrapper drops columns without updating
the graph. Corrected mode can retain constant predictors; expression profiles
constant across genes are rejected before undefined Pearson correlation.

## IDs, direction, and stability

Expression axes and phenotype IDs must be unique/nonmissing. The shared-gene
intersection follows bulk order, with counts of nonshared genes recorded.
Phenotype order is aligned exactly to bulk sample IDs. A nonnumeric binary
label requires explicit positive_label; no arbitrary factor ordering is used.

Scissor+/- use strict beta>0/beta<0. Values 0<abs(beta)<1e-5 are reported as
near-zero diagnostics without silently relabeling them. Nonunique solutions
cannot guarantee individual coefficient or sign stability; prediction and
objective agreement must supplement coefficient comparisons. The prototype
does not automatically prove uniqueness or certify biological reliability.

## Deliberate exclusions

No missing-value imputation, automatic graph creation, count/log transformation,
automatic expression downloads, performance claims for atlas scale, adaptive
wbeta/sgn, other penalties, L0 trimming, R runtime, or real-dataset biological
validation. Inzero=TRUE golden outputs establish that untrimmed Scissor inputs
are unchanged on the audited cases, but the unused trimming stage is not ported.
