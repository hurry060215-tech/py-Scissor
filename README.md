# py-Scissor

[![Tests](https://github.com/hurry060215-tech/py-Scissor/actions/workflows/tests.yml/badge.svg)](https://github.com/hurry060215-tech/py-Scissor/actions/workflows/tests.yml)
[English](README.md) | [中文](README.zh-CN.md)

**Phenotype-guided single-cell selection with an independent, R-free Python
APML1-Net kernel.** Supports continuous, binary and right-censored survival
phenotypes through Gaussian, Binomial and Cox models.

Core runtime dependencies are NumPy, SciPy and pandas. AnnData integration and
Gaussian Numba acceleration are optional. No R, data downloads, graph creation,
HVG selection, normalization-layer guessing or dependency installation occurs
on import or during a library call.

This is a **research pre-alpha**, not an official release by the Scissor authors.
The fixed reference is [Scissor 311560a1](https://github.com/sunduanchen/Scissor/tree/311560a1f160f665917e7da0fcf289914ee99773).
R is used only as a development oracle. See the [numerical semantics](docs/SEMANTICS.md)
and [measured validation limits](docs/VALIDATION_STATUS.md) before interpreting results.

## Install

Python 3.10 or newer. Install this repository, not the unrelated PyPI project
`pyscissor`. The repository is **py-Scissor**, distribution **scissor-py**, import
**scissor_py**. This project has **not been published on PyPI**.

```bash
python -m venv .venv
# Activate .venv using your shell, then:
python -m pip install "git+https://github.com/hurry060215-tech/py-Scissor.git"
```

For reproducible work, replace the branch with a reviewed commit SHA:
`git+https://github.com/hurry060215-tech/py-Scissor.git@<commit-sha>`.
For a local checkout: `python -m pip install .`; for a supplied wheel:
`python -m pip install /path/to/scissor_py-0.1.0.dev4-py3-none-any.whl`.

Optional extras from a checkout:

```bash
python -m pip install '.[integration-test]'  # AnnData/H5AD support and tests
python -m pip install '.[acceleration]'      # explicit corrected Gaussian Numba
python -m pip install '.[notebook]'          # tutorial dependencies, excluding OmicVerse
```

## Runnable quickstart

This deterministic **synthetic** example is [examples/quickstart.py](examples/quickstart.py).
It needs no network, R, AnnData or downloaded dataset. Its three lambdas and
three folds are chosen to keep the example small; they are not the defaults.

<!-- quickstart:start -->
```python
import numpy as np
import pandas as pd
from scipy import sparse
from scissor_py import scissor

rng = np.random.default_rng(123)
genes = [f"gene_{i}" for i in range(60)]
samples = [f"sample_{i}" for i in range(24)]
cell_ids = [f"cell_{i}" for i in range(8)]
bulk = pd.DataFrame(rng.gamma(2, 2, (60, 24)), index=genes, columns=samples)
cells = pd.DataFrame(rng.gamma(2, 2, (60, 8)), index=genes, columns=cell_ids)
phenotype = pd.Series(["wildtype", "mutant"] * 12, index=samples)
graph = sparse.diags([np.ones(7), np.ones(7)], [-1, 1], format="csr")
foldid = pd.Series(np.arange(24) % 3, index=samples)

result = scissor(
    bulk, cells, phenotype, graph, graph_cell_ids=cell_ids,
    family="binomial", positive_label="mutant", semantics="corrected",
    alpha=[.2, .6], lambdas=[.10, .04, .015], foldid=foldid,
    cutoff=.5, thresh=1e-10,
)
print(result.coefficients.to_frame().join(result.labels))
print("alpha:", result.alpha, "lambda:", result.lambda_min)
print("cutoff reached:", result.cutoff_reached)
```
<!-- quickstart:end -->

`Scissor+` means a positive coefficient toward the encoded `mutant` class;
`Scissor-` means a negative coefficient. `Background` means exactly zero.
An unmet cutoff is reported as `False` with a warning, rather than changing
coefficient signs or silently truncating the selected cells. Synthetic output
is an API demonstration, not biological evidence.

## Choose numerical semantics explicitly

- **`official` is the default.** Reproduces defined behavior of the pinned R
  implementation, including documented stopping and Cox lambda-label quirks.
  Potentially out-of-bounds Gaussian/Cox graph inputs are rejected.
- **`corrected` is opt-in.** Aligns graph nodes, uses the actual Gaussian
  curvature and exact-objective KKT stopping, and labels Cox lambdas consistently.
  This is not a claim of equivalence to unchanged official code.

The Binomial quickstart explicitly chooses the corrected logistic objective.
The same connected graph would be rejected in official Gaussian/Cox mode because
of the original wrapper indexing. The implementation never switches modes automatically. [Affected inputs and equations](docs/SEMANTICS.md)
are documented alongside [actual R reproductions](reference/README.md).

## Inputs, models and outputs

| Input | Required layout / meaning |
|---|---|
| `bulk` | DataFrame: genes × bulk samples, unique IDs |
| `single_cell` | DataFrame: genes × cells, unique IDs |
| `phenotype` | Sample-indexed Series; for Cox, DataFrame with `time`, `status` |
| `graph` | Explicit cell × cell raw adjacency; sparse input needs `graph_cell_ids` |
| `foldid` | Optional explicit fold vector, or sample-indexed Series at the high-level API |

Gaussian expects numeric outcomes. Binary strings require `positive_label`;
numeric 0/1 defaults to positive class 1. Cox uses `status=1` for an event and
0 for right-censoring, with Breslow ties. Positive Cox coefficients increase
log hazard within the fitted model; this is not a clinical conclusion.

The pipeline aligns IDs, jointly quantile-normalizes shared genes, computes
bulk-sample × cell Pearson similarities, binarizes and degree-normalizes the
supplied graph, cross-validates lambda, cold-refits `lambda.min`, then searches
alpha in order. Labels retain strict coefficient signs.

`ScissorResult` exposes `coefficients`, `labels`, `positive_cells`, `negative_cells`,
`alpha`, `lambda_min`, `intercept`, `cutoff_reached`, visited `cv_results`,
`prepared` intermediates and audit `metadata`.

## API map and examples

| Task | Python API |
|---|---|
| Expression-to-selection workflow | `scissor` |
| Explicit preprocessing / similarity | `prepare` |
| APML1-Net path, positional similarity matrix | `fit_path` |
| Fixed-graph, family-specific CV | `cross_validate` |
| Optional expression-layer adapter | `scissor_from_anndata` |
| JSON audit, annotation, atomic H5AD save | `result_summary`, `annotate_anndata`, `write_h5ad_atomic` |

The first four cover the numerical roles of upstream `Scissor`, preprocessing,
`APML1` and `APML1.cv`; they do not reproduce Seurat preprocessing or unrelated
APML1 penalties. Read the [API guide](docs/API.md), [AnnData contract](docs/ANNDATA.md),
[three-family example](examples/synthetic.py) and [tutorial notebook](tutorials/scissor_walkthrough.ipynb).
A proposed optional [`ov.single.Scissor` adapter](integration/omicverse/README.md)
is included for review; upstream OmicVerse availability must not be assumed.

## What has been validated?

- Saved synthetic fixtures generated by the pinned, unchanged R source run in
  ordinary CI without R; missing required fixtures fail rather than count as parity.
- Full public LUAD preprocessing retained all 4,102 cells, 471 survival samples
  and 498 TP53 samples, with 22,323 shared genes per cohort. The explicit
  PCA20/kNN20 union graph is not the historical Seurat SNN graph.
- All nine corrected family/alpha comparisons passed the predeclared tolerances
  when **both Python and the independent R oracle used `thresh=1e-10`**:
  three alphas, twenty lambdas, five fixed folds and cold refits. Unchanged
  official Binomial also passed all three alphas.
- Six coefficient comparisons at the looser **1e-8 validation threshold failed**
  and remain recorded. The API default is still **1e-7**. These are distinct
  settings; no default-precision or universal-equivalence guarantee is implied.
- Gaussian's full-data phenotype is a known synthetic signal on the complete
  real design. No biological effectiveness, causal or diagnostic claim is made.

See [exact counts, tolerances and limitations](docs/VALIDATION_STATUS.md), the
[aggregate comparison results](docs/luad/numerical_summary.json) and
[reproduction guide](docs/REPRODUCING.md). Ordinary CI is not the full LUAD audit.
The tutorial's exact code cells passed plain-Python execution; Jupyter-kernel
execution remains unverified in the original restricted environment.

Dense expression and similarity still consume memory. `max_dense_bytes` is an
estimate, not a process-RSS limit. Graph storage stays sparse; no dense cell²
network or sample² Cox risk matrix is required.

## Development, citation and license

```bash
python -m pip install -e '.[test,integration-test,acceleration]'
python -m pytest -q
python -m build
python scripts/check_distribution.py dist/*
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for numerical-change requirements and
[SECURITY.md](SECURITY.md) for safe reporting. GitHub Actions tests core and
optional-dependency paths, packages the project, and checks a clean installed
wheel outside the source tree. It does not publish to PyPI or merge changes.

Please cite the original method: Sun D et al., *Identifying phenotype-associated
subpopulations by integrating bulk and single-cell sequencing data*, Nature
Biotechnology, [doi:10.1038/s41587-021-01091-3](https://doi.org/10.1038/s41587-021-01091-3).
[CITATION.cff](CITATION.cff) separates this independent software from that paper.

**GPL-3.0-or-later.** [LICENSE](LICENSE) and [NOTICE](NOTICE) retain upstream
copyright, the fixed source provenance and modification history. Licensed
upstream sources are included in the source distribution. No upstream
endorsement, affiliation or ownership of the original method is claimed.
