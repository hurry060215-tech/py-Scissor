# Full LUAD input audit — 2026-10-10

## Verified inputs

All five files were fetched in full directly from the original endpoints in the
[official Scissor tutorial](https://sunduanchen.github.io/Scissor/vignettes/Scissor_Tutorial.html).
No mirrors, subsampling, replacement cohort or automatic library download was
used. Exact URLs, byte counts, timestamps and SHA-256 values are in
`download_manifest.json`.

- Single-cell source: 33,694 genes × **4,102 cells**, all cell IDs unique
- Survival bulk: 56,716 source rows × **471 samples**, all samples retained;
  121 events and 350 censored, 351 distinct observed times, range 1–6,812
- TP53 bulk: 56,716 source rows × **498 samples**, all samples retained;
  255 mutant (positive class 1), 243 wild-type (0)
- Both phenotype files already match their own bulk sample ID order. The two
  cohorts were not intersected into a smaller common sample set
- Both shared-gene intersections contain 22,323 unique genes, in bulk order
- Expression and phenotypes are finite, every cell has positive library size,
  and every normalized profile has nonzero variance

## Duplicate genes require an explicit policy

Each bulk file has 1,485 duplicate gene-name rows. Twenty-two of these rows have
names found in the single-cell data, and **all 22 have different expression from
the first row with that name**. This is a material data property.

The validation script explicitly selects the first source row per gene name,
then takes the common-gene intersection. This matches the fixed R Scissor
`intersect` followed by character matrix indexing, and is not a biological
claim that averaging, summing or another reconciliation would be equivalent.
The locally generated `duplicate_gene_map.tsv` records every source/retained row and equality check.
The package's normal API still rejects duplicate IDs and never performs this
reconciliation silently. No sample or cell is removed by this policy.

## Frozen preprocessing and explicit graph

The validation script applies log1p(CP10K) to the full single-cell count matrix,
using each cell's total across all source genes. The bulk matrices are used as
provided, without a further transform. Per cohort, all common genes undergo
joint bulk-plus-cell quantile normalization and Pearson similarity.

The cell graph is an explicit all-cell graph: centered, unscaled PCA of all
33,694 log-normalized genes, 20 components, followed by Euclidean 20-nearest
neighbors and union symmetrization. There is no HVG selection, random cell
sampling, batch correction or dropped cell. The graph has 117,170 directed
nonzeros and no isolated node; neighbor-cutoff distance ties were absent.
This graph is **not the historical Seurat SNN graph**. Consequently this kernel
validation does not claim to reproduce the old tutorial's selected-cell counts.
All R/Python families receive this same exact graph.

The frozen manifest includes every gene/sample/cell ID, actual 1-based five-fold
assignments, outcomes, three alpha values (0.05, 0.2, 0.6), and each family's
20-point explicit lambda path from its full-data lambda maximum to 0.01 times
that maximum. Cox folds are event-stratified; Gaussian uses the same 471-sample
folds; TP53 has its own 498-sample class-stratified folds. Equal seeds are not
being used as a substitute for actual shared fold IDs.

Gaussian uses a declared synthetic continuous outcome on the full exp1 design:
`2 + 2*X[,1] - 1.5*X[,138] + X[,778] + 0.02*sin(1:n)` (R's 1-based indexing).
It has no biological phenotype interpretation.

完整 input_manifest.json 含实际样本/phenotype清单，仅由准备脚本在本地生成，不随公共源码发布。其 SHA-256:
`001f3cd90dd02953ed7379e01827699bc4da7a621ef055500289d4e908991eee`.
The locally generated `source_provenance.json` freezes preparation-code hashes. Full intermediate
binaries are retained locally with manifest checksums and can be regenerated
with the explicit development scripts; they are not bundled into the wheel.

Full preparation took 40.1 seconds in this environment, with measured process
peak RSS 4,564,004,864 bytes. The 4.45 GB intermediate storage includes both
complete joint/normalized matrices and both complete similarity matrices. These
are measurements for this dataset, not general speed or memory guarantees.

## Not yet established

At this checkpoint, model/CV/refit R comparisons are running. Connected-graph
official Gaussian/Cox execution remains disallowed because of verified upstream
unsafe indexing; rejecting those calls is not a successful numerical comparison.
Unchanged-source binomial and separately labelled corrected mathematical
oracles are distinct evidence. No biological reliability test or clinical
interpretation is claimed.
