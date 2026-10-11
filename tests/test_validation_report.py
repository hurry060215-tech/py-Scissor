"""Public aggregate reports must preserve failures and omit individual data."""
import importlib.util
import json
from pathlib import Path


def report_builder():
    path = Path(__file__).parents[1]/"scripts"/"summarize_luad_validation.py"
    spec = importlib.util.spec_from_file_location("validation_summary", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.summarize


def test_aggregate_keeps_failure_and_omits_private_result_details(tmp_path):
    metric = {"max_absolute_error": 2e-4, "outside_tolerance": 2}
    check = {key: metric for key in ("coefficients", "intercepts", "linear_predictions_all_samples", "objective")}
    check.update(strict_sign_differences=3, stable_sign_differences=1,
                 near_zero_sign_examples=[{"cell_id": "PRIVATE_CELL", "python": .01}])
    report = {"family": "gaussian", "semantics": "corrected", "alpha": .2,
              "python_engine": "numpy", "python_threshold": 1e-8, "r_threshold": 1e-8,
              "manifest_sha256": "frozen", "tolerances": {"coefficients": {"atol": 1e-5, "rtol": 1e-4}},
              "passed": False, "contract": {"ids_match": True}, "paths": {"full": check},
              "stable_coefficient_threshold": 1e-5, "cv": {"lambda_min_exact": True},
              "sample_ids": ["PRIVATE_SAMPLE"], "phenotype": [1]}
    path = tmp_path/"case.json"
    path.write_text(json.dumps(report))
    result = report_builder()([path])
    text = json.dumps(result, allow_nan=False)
    assert "PRIVATE_" not in text and "near_zero_sign_examples" not in text
    assert result["comparisons"][0]["passed"] is False
    matching = [row for row in result["expected_same_threshold_grid"] if row["family"] == "gaussian" and row["alpha"] == .2 and row["threshold"] == 1e-8]
    assert matching[0]["status"] == "failed"
    assert all(row["status"] != "passed" for row in result["expected_same_threshold_grid"])


def test_empty_aggregate_never_infers_numerical_pass():
    result = report_builder()([])
    assert len(result["expected_same_threshold_grid"]) == 21
    assert all(row["status"] == "not_reported" for row in result["expected_same_threshold_grid"])
