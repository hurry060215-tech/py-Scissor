# Reproducing validation

## Ordinary Python checks: no R or compiler

```sh
python -m pip install -e '.[test,integration-test,acceleration]'
python -m pytest -q
python examples/quickstart.py
python examples/synthetic.py
python -m build
python scripts/check_distribution.py dist/*
```

Twenty-two required JSON golden fixtures and two synthetic sweep fixtures stay
in `reference/`. Their expected numerical results were generated previously by
R, but Python tests only read JSON. No R process, compilation or data download
is involved. CI separately tests a checkout, the extracted source distribution,
and a clean wheel installed outside the source tree. The wheel smoke forbids
network/child processes and removes R/compilers from PATH.

Optional AnnData/Numba tests may skip when their extras are absent; the optional
CI job installs both. Skips are not counted as successful optional validation.

## Deliberately regenerate R development references

R/C++ reference implementation and generators are outside the core repository's
main tree and both package formats. Use the immutable development archive:

```sh
git clone https://github.com/hurry060215-tech/py-Scissor.git scissor-reference
cd scissor-reference
git checkout --detach 2cb65ba5af3e1ffac05226f66d0e0ddf1eb7dc50
```

This is an explicit separate development checkout, not the installation recipe.
Read its [small-fixture audit](https://github.com/hurry060215-tech/py-Scissor/blob/2cb65ba5af3e1ffac05226f66d0e0ddf1eb7dc50/reference/README.md)
and [full-data oracle instructions](https://github.com/hurry060215-tech/py-Scissor/blob/2cb65ba5af3e1ffac05226f66d0e0ddf1eb7dc50/reference/luad/README.md).
They retain all original scripts, sources, GPL licenses and fixed upstream hashes.
The original oracle used R 4.5.0, Rcpp/RcppEigen/Matrix/jsonlite and pinned
preprocessCore 1.75.1. Install those only in an isolated development environment.

Ordinary tests never fetch this archive. Regeneration should write to a separate
workspace; compare saved output hashes before considering any golden update.
Never execute the known unsafe official Gaussian/Cox graph-index cases.

## Full-data numerical experiment

[Current evidence and retained failures](VALIDATION_STATUS.md) cover all 4,102
cells, the complete supplied bulk cohorts, three alphas, twenty lambdas and five
fixed folds with cold refits. Gaussian uses a known synthetic phenotype on the
full real design. Ordinary CI does not rerun this experiment.

The prepared corrected baseline threshold was **1e-8**, whereas the API default
remains **1e-7**. Tightened validation explicitly used **1e-10 on both sides**
with identical data, graph, folds and paths. Six baseline coefficient failures
remain recorded; no tolerance was loosened. The unchanged official Binomial
reference is distinct from the independently corrected oracle.

The core keeps explicit Python acquisition/preparation/comparison helpers for
researchers who choose to reproduce the experiment. They do not run on import
or installation. Full oracle drivers and interruption certificates belong to the
separate archive. Original downloads were about 254 MB; complete preprocessing
wrote about 4.45 GB and peaked at 4.56 GB RSS. Reserve resources and serialize
large tasks. No raw LUAD matrices, individual phenotypes or real per-cell outputs
are bundled in either the main tree or its distributions.

OmicVerse, full tutorial and skill validation are separate contribution tasks,
not hidden dependencies of this kernel. See [architecture](ARCHITECTURE.md).
