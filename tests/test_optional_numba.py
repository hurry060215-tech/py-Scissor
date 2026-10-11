import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse

from scissor_py import ScissorError, cross_validate, fit_path


@pytest.mark.parametrize("alpha", [0., .05, .6, 1.])
def test_numba_same_full_path_and_cv(alpha):
    pytest.importorskip("numba")
    rng = np.random.default_rng(731)
    X = rng.normal(size=(30, 14))*.4
    graph = sparse.diags([np.ones(13), np.ones(13)], [-1, 1], format="csr")
    y = X[:, 2]*2-X[:, 7]+.1*rng.normal(size=30)
    kw = dict(family="gaussian", semantics="corrected", alpha=alpha,
              lambdas=[.3, .1, .03, .01], thresh=1e-11)
    a = cross_validate(X, y, graph, foldid=np.arange(30) % 3, engine="numpy", **kw)
    b = cross_validate(X, y, graph, foldid=np.arange(30) % 3, engine="numba", **kw)
    np.testing.assert_allclose(a.path.coefficients, b.path.coefficients, atol=1e-9, rtol=1e-8)
    np.testing.assert_allclose(a.path.intercepts, b.path.intercepts, atol=1e-10)
    np.testing.assert_allclose(a.fold_loss, b.fold_loss, atol=1e-10)
    np.testing.assert_array_equal(a.valid, b.valid)
    assert a.lambda_min == b.lambda_min and a.lambda_1se == b.lambda_1se
    assert b.path.metadata["engine"] == "numba"
    assert b.path.metadata["jit_fastmath"] is False and b.path.metadata["jit_parallel"] is False


def test_missing_numba_is_explicit_and_default_has_no_import():
    # Independent process also prevents an installed/previously imported Numba
    # from hiding an unwanted default import or silent fallback.
    root = Path(__file__).resolve().parents[1]
    program = """
import builtins, sys
original = builtins.__import__
def guarded(name, *a, **kw):
    if name == 'numba' or name.startswith('numba.'):
        raise ModuleNotFoundError("No module named 'numba'", name='numba')
    return original(name, *a, **kw)
builtins.__import__ = guarded
import numpy as np
from scipy import sparse
from scissor_py import fit_path, ScissorError
x=np.array([[1.,2.],[2.,-1.],[3.,1.]])
kw=dict(family='gaussian',semantics='corrected',lambdas=[.1])
r=fit_path(x,np.array([1.,2.,3.]),sparse.csr_matrix((2,2)),**kw)
assert r.converged.all() and 'numba' not in sys.modules
try:
    fit_path(x,np.array([1.,2.,3.]),sparse.csr_matrix((2,2)),engine='numba',**kw)
except ScissorError as e:
    assert 'optional' in str(e)
else:
    raise AssertionError('Missing dependency was hidden')
print('NO_NUMBA_DEFAULT_OK')
"""
    env = dict(os.environ, PYTHONPATH=str(root/"src"))
    result = subprocess.run([sys.executable, "-c", program], env=env, capture_output=True, text=True, check=True)
    assert "NO_NUMBA_DEFAULT_OK" in result.stdout


@pytest.mark.parametrize("family,semantics", [("gaussian", "official"), ("binomial", "corrected"), ("cox", "corrected")])
def test_numba_unsupported_semantics_never_silently_fallback(family, semantics):
    X = np.array([[1., 2.], [2., -1.], [3., 1.], [4., 2.]])
    y = np.array([0., 1., 0., 1.]) if family != "cox" else np.array([[1., 1.], [2., 0.], [3., 1.], [4., 1.]])
    with pytest.raises(ScissorError, match="corrected Gaussian only"):
        fit_path(X, y, sparse.csr_matrix((2, 2)), family=family, semantics=semantics, engine="numba", lambdas=[.1])


def test_numba_does_not_hide_nonconvergence():
    pytest.importorskip("numba")
    X = np.random.default_rng(712).normal(size=(20, 8))
    r = fit_path(X, X[:, 0]-X[:, 3], sparse.csr_matrix((8, 8)), family="gaussian",
                 semantics="corrected", engine="numba", lambdas=[.001], maxit=1, thresh=1e-15)
    assert not r.converged[0]
