# Contributing

This is an independent GPL-3.0-or-later research implementation. Small focused
changes are easier to review than new frameworks or unrelated refactors.

## Local checks

Use an isolated Python 3.10+ environment and install from this checkout:

```sh
python -m pip install -e '.[test,integration-test,acceleration]'
python -m pytest -q
python examples/quickstart.py
python examples/synthetic.py
python -m build
python scripts/check_distribution.py dist/*
```

Ordinary tests use committed synthetic R fixtures and require neither R nor
network access. Optional dependency skips should be reported as skips, not
as evidence that AnnData or Numba passed. CI also checks a freshly installed
wheel outside the source tree. Full LUAD/R runs are deliberate development
experiments described in [REPRODUCING.md](docs/REPRODUCING.md), not routine CI.

## Numerical changes

- Preserve explicit `official` and `corrected` semantics. Do not silently repair
  the compatibility mode or reproduce undefined memory access.
- Include a minimal synthetic regression with explicit graph, IDs, phenotype
  encoding, actual folds, alpha/lambda path and stopping settings.
- Compare coefficients, predictions, objectives, CV choices and stable signs.
  Preserve failed comparisons. Do not loosen acceptance limits to obtain a pass.
- Keep licensed upstream source unchanged. New oracle modifications must be
  identified, separately versioned and distinguishable from official results.
- Do not add R, Scanpy, OmicVerse or Numba to the core runtime dependencies.
  Optional imports should remain lazy, with explicit missing-dependency errors.
- Check no-half-write behavior, H5AD roundtrip and package contents when changing
  integration or serialization. Do not place input expression in `uns`.

## Reports and pull requests

Use the repository's [issue tracker](https://github.com/hurry060215-tech/py-Scissor/issues)
for non-sensitive bugs or proposals. Include package/dependency versions,
family, semantics, a small synthetic reproducer, expected behavior and actual
error or diagnostic. Do not post credentials, individual clinical records,
real patient identifiers or large raw datasets. See [SECURITY.md](SECURITY.md).

Describe the reason for a change and the checks actually run. A green unit-test
job does not establish full scientific parity or biological effectiveness.
Preserve GPL attribution and the fixed-source notices. Contributions should be
compatible with this project's GPL-3.0-or-later licensing.
