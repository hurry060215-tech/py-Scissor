import json
from pathlib import Path
import numpy as np
import pytest
from scissor_py import fit_path, scissor, ScissorError


@pytest.mark.parametrize("semantics",["official","corrected"])
def test_zero_event_time_is_valid_and_time_shift_invariant(semantics):
    d=json.loads((Path(__file__).parents[1]/"reference/fixtures/cox_default_explicit_cv.json").read_text())["inputs"]
    y=np.asarray(d["y"],dtype=float)
    baseline=fit_path(d["x"],y,d["Omega"],family="cox",alpha=.4,lambdas=[.12,.06],semantics=semantics)
    y[:,0]-=y[:,0].min()
    shifted=fit_path(d["x"],y,d["Omega"],family="cox",alpha=.4,lambdas=[.12,.06],semantics=semantics)
    np.testing.assert_allclose(shifted.coefficients,baseline.coefficients,atol=1e-12)
    assert shifted.converged.all()


@pytest.mark.parametrize("family",["gaussian","binomial","cox"])
def test_boolean_maxit_is_not_an_iteration_count(family):
    d=json.loads((Path(__file__).parents[1]/f"reference/fixtures/{family}_default_explicit_cv.json").read_text())["inputs"]
    with pytest.raises(ScissorError,match="maxit"):
        fit_path(d["x"],d["y"],d["Omega"],family=family,alpha=.4,lambdas=[.12],maxit=True)
@pytest.mark.parametrize("part", ["X", "y"])
def test_numeric_complex_solver_input_is_not_silently_cast(part):
    X = np.array([[1., 2.], [2., 3.], [3., 1.]])
    y = np.array([0., 1., 2.])
    if part == "X":
        X = X.astype(complex)+1j
    else:
        y = y.astype(complex)+1j
    with pytest.raises(ScissorError, match="real-valued"):
        fit_path(X, y, np.zeros((2, 2)), semantics="corrected", lambdas=[.1])


@pytest.mark.parametrize("parameters", [
    {"lambdas": np.array([.1+1j])},
    {"alpha": .4+1j, "lambdas": [.1]},
    {"alpha": [.4], "lambdas": [.1]},
    {"rlambda": np.complex128(.1+.2j)},
    {"thresh": np.complex128(1e-8+1j), "lambdas": [.1]},
])
def test_path_parameters_cannot_silently_drop_imaginary_or_vector_parts(parameters):
    X = np.array([[1., 2.], [2., 3.], [3., 1.]])
    with pytest.raises(ScissorError):
        fit_path(X, np.arange(3.), np.zeros((2, 2)), semantics="corrected", **parameters)


@pytest.mark.parametrize("alpha", [np.array([.4+1j]), [[.4, .5]]])
def test_high_level_alpha_shape_and_complex_values_fail_before_preparation(alpha):
    with pytest.raises(ScissorError, match="alpha"):
        scissor(None, None, None, None, alpha=alpha)


def test_complex_cutoff_cannot_be_reinterpreted_as_real():
    with pytest.raises(ScissorError, match="cutoff"):
        scissor(None, None, None, None, cutoff=np.complex128(.2+1j))
