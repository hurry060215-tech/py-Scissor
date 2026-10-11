"""Golden outputs generated once by the actual pinned R source; no R in CI."""
import json
from pathlib import Path
import numpy as np
import pytest
import pandas as pd
import warnings
from scissor_py import cross_validate, fit_path, ScissorError, scissor, quantile_normalize

FIXTURES = Path(__file__).resolve().parents[1]/"reference"/"fixtures"


@pytest.mark.parametrize("family", ["gaussian", "binomial", "cox"])
@pytest.mark.parametrize("precision", ["default", "tight"])
def test_official_explicit_paths_and_cv(family, precision):
    f = json.loads((FIXTURES/f"{family}_{precision}_explicit_cv.json").read_text())
    inputs, expected = f["inputs"], f["output"]["result"]
    assert f["output"]["ok"]
    result = cross_validate(inputs["x"], inputs["y"], inputs["Omega"], family=family,
                            alpha=inputs["alpha"], lambdas=inputs["lambda"], foldid=inputs["foldid"], thresh=inputs["thresh"])
    actual = np.vstack([result.path.intercepts, result.path.coefficients]) if family == "binomial" else result.path.coefficients
    np.testing.assert_allclose(actual, expected["Beta"], atol=1e-5, rtol=1e-4)
    np.testing.assert_allclose(result.mean_loss, [v["cvm"] for v in expected["fit"]], atol=1e-6, rtol=1e-5)
    np.testing.assert_allclose(result.standard_error, [v["cvse"] for v in expected["fit"]], atol=1e-6, rtol=1e-5)
    assert result.lambda_min == expected["lambda.min"]
    assert result.valid.all()
    np.testing.assert_array_equal(np.sign(actual), np.sign(expected["Beta"]))


def test_gaussian_default_reproduces_known_nonoptimal_coefficient():
    d = json.loads((FIXTURES/"gaussian_curvature_audit.json").read_text())
    i = d["inputs"]
    official = fit_path(i["x"], i["y"], [[0]], alpha=i["alpha"], lambdas=[i["lambda"]])
    corrected = fit_path(i["x"], i["y"], [[0]], alpha=i["alpha"], lambdas=[i["lambda"]], semantics="corrected", thresh=1e-12)
    np.testing.assert_allclose(official.coefficients, d["official_default"]["Beta"], atol=1e-12)
    assert abs(official.coefficients[0, 0]-d["closed_form_coefficient"]) > .001
    np.testing.assert_allclose(corrected.coefficients[0, 0], d["closed_form_coefficient"], atol=1e-11)


@pytest.mark.parametrize("family", ["gaussian", "cox"])
def test_official_unsafe_graph_rejected(family):
    d = json.loads((FIXTURES/f"{family}_default_explicit_cv.json").read_text())["inputs"]
    g = np.array(d["Omega"]); g[0, -1] = g[-1, 0] = 1
    with pytest.raises(ScissorError, match="out-of-bounds"):
        fit_path(d["x"], d["y"], g, family=family, alpha=.4, lambdas=[.1])


def test_actual_preprocesscore_ties():
    d = json.loads((FIXTURES/"preprocessing_ties.json").read_text())
    np.testing.assert_allclose(quantile_normalize(d["inputs"]["joint_expression"]), d["normalized"], atol=1e-8, rtol=1e-6)


@pytest.mark.parametrize("family", ["gaussian", "binomial"])
def test_automatic_lambda_path_and_cv(family):
    d = json.loads((FIXTURES/f"{family}_default_automatic_cv.json").read_text())
    i, e = d["inputs"], d["output"]["result"]
    r = cross_validate(i["x"], i["y"], i["Omega"], family=family, alpha=i["alpha"], foldid=i["foldid"], nlambda=i["nlambda"])
    beta = np.vstack([r.path.intercepts, r.path.coefficients]) if family == "binomial" else r.path.coefficients
    np.testing.assert_allclose(beta, e["Beta"], atol=1e-5, rtol=1e-4)
    np.testing.assert_allclose(r.path.lambdas, [v["lambda"] for v in e["fit"]], atol=1e-12)
    np.testing.assert_allclose(r.mean_loss, [v["cvm"] for v in e["fit"]], atol=1e-7, rtol=1e-6)
    assert r.lambda_min == pytest.approx(e["lambda.min"], abs=1e-12)


@pytest.mark.parametrize("family", ["gaussian", "binomial", "cox"])
def test_expression_to_selection_against_actual_r_harness(family):
    d = json.loads((FIXTURES/f"{family}_expression_harness.json").read_text())
    i, e = d["inputs"], d["result"]
    bulk = pd.DataFrame(i["bulk"], index=i["bulk_gene_ids"], columns=i["sample_ids"])
    cells = pd.DataFrame(i["sc"], index=i["sc_gene_ids"], columns=i["cell_ids"])
    phenotype = (pd.DataFrame(i["phenotype"], index=i["phenotype_sample_ids"], columns=["time", "status"])
                 if family == "cox" else pd.Series(i["phenotype"], index=i["phenotype_sample_ids"]))
    graph = pd.DataFrame(i["graph"], index=i["graph_cell_ids"], columns=i["graph_cell_ids"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        r = scissor(bulk, cells, phenotype, graph, family=family, alpha=i["alpha"], cutoff=i["cutoff"],
                     nlambda=i["nlambda"], foldid=i["foldid"])
    np.testing.assert_allclose(np.column_stack([r.prepared.normalized_bulk, r.prepared.normalized_single_cell]),
                               d["intermediates"]["normalized"], atol=1e-8, rtol=1e-6)
    np.testing.assert_allclose(r.prepared.similarity, d["intermediates"]["similarity"], atol=1e-8, rtol=1e-6)
    np.testing.assert_allclose(r.coefficients, e["coefficients"], atol=1e-5, rtol=1e-4)
    assert r.labels.tolist() == e["labels"]
    assert r.alpha == e["alpha"]
    assert r.lambda_min == pytest.approx(e["lambda"], abs=1e-10)
    assert r.cutoff_reached == e["cutoff_reached"]


def test_cox_auto_reports_published_ratio_but_records_actual_solver_lambda():
    d = json.loads((FIXTURES/"cox_lambda_audit.json").read_text())
    i = d["inputs"]
    r = fit_path(i["x"], i["y"], i["graph"], family="cox", alpha=i["alpha"],
                 nlambda=4, rlambda=.008, thresh=1e-10, ifast=False)
    np.testing.assert_allclose(r.lambdas, i["ratios"], atol=1e-12)
    np.testing.assert_allclose(r.metadata["solver_lambdas"], d["actual_lambda"], atol=1e-12)
    np.testing.assert_allclose(r.coefficients, d["automatic"]["Beta"], atol=1e-10)
    corrected = fit_path(i["x"], i["y"], i["graph"], family="cox", alpha=i["alpha"], semantics="corrected", nlambda=4, rlambda=.008)
    np.testing.assert_allclose(corrected.lambdas, d["actual_lambda"], atol=1e-12)

