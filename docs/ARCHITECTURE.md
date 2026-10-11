# Python core and development boundaries

## What actually runs

`scissor_py` contains nine Python modules. It is a Python reimplementation of
APML1-Net, not a wrapper around R or an Rcpp shared library.

| Component | Python implementation |
|---|---|
| ID alignment, joint quantile normalization, Pearson, sparse graph | `_preprocessing.py` |
| Gaussian and Binomial coordinate/optimization paths | `_solvers.py` |
| Cox likelihood, score, tied-event risk sets and paths | `_cox.py` |
| Lambda paths, family-specific CV, refit and alpha search | `_api.py` |
| Input contracts and result records | `_validation.py`, `_types.py` |
| Explicit AnnData annotation and atomic H5AD writing | `_anndata.py` |
| Optional explicitly requested Gaussian JIT | `_gaussian_numba.py` |
| Public exports | `__init__.py` |

The default engine uses NumPy/SciPy/pandas. These dependencies themselves use
compiled numerical libraries; that does not make this an R backend. No R,
rpy2, sourceCpp, subprocess or C++ bridge is used by the core. Numba is optional,
imported only for explicit `engine='numba'`; no JIT is needed by the default path.

Both wheel and source distribution exclude R/C/C++ code, reference compiler
settings, full tutorials, skills and OmicVerse contribution patches. The archive
checker enforces this for all paths, including files hidden inside the package.
The installed-wheel smoke test empties PATH, forbids child processes/network,
and runs all three families plus the README example.

## Why the first repository snapshot contained R/C++

The initial public snapshot bundled development references with the core:
26 R files and two C++ files, 4,608 lines in total. Six were unchanged upstream
source files (3,366 lines); the remainder generated fixtures or independently
checked the mathematics. They were never wheel runtime components, but were
included in the earlier source distribution. That layout mixed responsibilities.

Version 0.1.0.dev5 separates them without rewriting Git history or hiding language
statistics. The complete previous snapshot remains at
[2cb65ba5](https://github.com/hurry060215-tech/py-Scissor/tree/2cb65ba5af3e1ffac05226f66d0e0ddf1eb7dc50),
also named `validation-reference`. It is an opt-in development archive, not the
installation branch. Main retains the 22 golden JSON fixtures, two additional
synthetic sweep fixtures, hashes, historical aggregate evidence and GPL notices.
All numerical and AnnData regression tests remain; missing required fixtures fail.

## Repository responsibilities

- **py-Scissor:** independent Python kernel, API, tests and small runnable examples.
- **[omicverse](https://github.com/omicverse/omicverse):** optional `ov.single.Scissor`
  wrapper contribution, separately reviewed; not installed by this core package.
- **[omicverse-tutorials](https://github.com/omicverse/omicverse-tutorials):** full
  workflow notebook contribution.
- **[omicverse-skills](https://github.com/omicverse/omicverse-skills):** workflow
  skill contribution.

The external contributions have been separated for their own repositories.
This does not claim they are already submitted, merged or available upstream.
Historical interface/notebook checks in validation reports remain clearly dated;
they are not core CI results. See [reproduction](REPRODUCING.md) and
[license/provenance](../NOTICE).
