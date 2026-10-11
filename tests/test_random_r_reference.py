"""Independent extra R oracle sweep, including unequal CV fold sizes."""
import json
from pathlib import Path
import numpy as np
import pytest
from scissor_py import fit_path, cross_validate, ScissorError

ROOT = Path(__file__).resolve().parents[1]/"reference"/"review"
PATH_CASES = json.loads((ROOT/"random_paths.json").read_text())
CV_CASES = json.loads((ROOT/"unequal_fold_cv.json").read_text())


@pytest.mark.parametrize("case", PATH_CASES)
def test_random_official_paths(case):
    i, e = case["inputs"], case["output"]
    args = dict(family=i["family"], alpha=i["alpha"], lambdas=i["lambda"], thresh=i["thresh"], maxit=i["maxit"])
    if "error" in e:
        assert i["family"] == "cox" and len(i["Omega"]) == 1
        with pytest.raises(ScissorError, match="p=1"):
            fit_path(i["x"], i["y"], i["Omega"], **args)
        return
    r = fit_path(i["x"], i["y"], i["Omega"], **args)
    beta = np.vstack([r.intercepts, r.coefficients]) if i["family"] == "binomial" else r.coefficients
    np.testing.assert_allclose(beta, e["Beta"], atol=1e-7, rtol=1e-7)
    np.testing.assert_array_equal(np.sign(beta), np.sign(e["Beta"]))


@pytest.mark.parametrize("case", CV_CASES)
def test_random_unequal_fold_cv(case):
    i, e = case["inputs"], case["output"]
    r = cross_validate(i["x"], i["y"], i["Omega"], family=i["family"], alpha=i["alpha"],
                        lambdas=i["lambda"], foldid=i["foldid"], thresh=i["thresh"], maxit=i["maxit"])
    beta = np.vstack([r.path.intercepts, r.path.coefficients]) if i["family"] == "binomial" else r.path.coefficients
    np.testing.assert_allclose(beta, e["Beta"], atol=1e-7, rtol=1e-7)
    np.testing.assert_allclose(r.mean_loss, [v["cvm"] for v in e["fit"]], atol=1e-7, rtol=1e-7)
    np.testing.assert_allclose(r.standard_error, [v["cvse"] for v in e["fit"]], atol=1e-7, rtol=1e-7)
    assert r.lambda_min == e["lambda.min"]
