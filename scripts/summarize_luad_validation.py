"""Export aggregate numerical checks without any sample IDs or cell results.

Reads existing comparison reports only. Missing reports remain not_reported;
this tool never fits, downloads, changes tolerances, or infers a pass.
"""
import argparse
import hashlib
import json
from pathlib import Path


def summarize(paths):
    records = []
    for path in sorted(paths):
        payload = path.read_bytes()
        report = json.loads(payload)
        checks = list(report["paths"].values())
        metrics = {}
        for name in ("coefficients", "intercepts", "linear_predictions_all_samples", "objective"):
            errors = [check[name]["max_absolute_error"] for check in checks]
            metrics[name] = {
                "max_absolute_error": None if None in errors else max(errors),
                "outside_tolerance": sum(check[name]["outside_tolerance"] for check in checks),
            }
        records.append({
            "comparison_report": path.name,
            "comparison_report_sha256": hashlib.sha256(payload).hexdigest(),
            "family": report["family"], "semantics": report["semantics"],
            "alpha": report["alpha"], "python_engine": report["python_engine"],
            "python_threshold": report["python_threshold"], "r_threshold": report["r_threshold"],
            "same_requested_tolerance": report["python_threshold"] == report["r_threshold"],
            "manifest_sha256": report["manifest_sha256"],
            "tolerances": report["tolerances"],
            "passed": report["passed"], "contract_passed": all(report["contract"].values()),
            "path_count": len(checks), "path_metrics": metrics,
            "strict_sign_differences": sum(check["strict_sign_differences"] for check in checks),
            "stable_sign_differences": sum(check["stable_sign_differences"] for check in checks),
            "stable_coefficient_threshold": report["stable_coefficient_threshold"],
            "cv_checks": report["cv"],
        })
    grid = []
    for semantics, families, threshold in (
        ("official", ("binomial",), 1e-7),
        ("corrected", ("gaussian", "binomial", "cox"), 1e-8),
        ("corrected", ("gaussian", "binomial", "cox"), 1e-10),
    ):
        for family in families:
            for alpha in (.05, .2, .6):
                matches = [record for record in records if
                           record["semantics"] == semantics and record["family"] == family and
                           record["alpha"] == alpha and record["python_threshold"] == threshold and
                           record["r_threshold"] == threshold]
                grid.append({"family": family, "semantics": semantics, "alpha": alpha,
                             "threshold": threshold,
                             "status": "not_reported" if not matches else
                             ("passed" if all(record["passed"] for record in matches) else "failed"),
                             "reports": [record["comparison_report"] for record in matches]})
    return {"schema_version": "scissor-py-luad-aggregate/1",
            "scope": "full 4102-cell explicit-graph experiment; 471 survival and 498 TP53 samples; synthetic Gaussian phenotype",
            "grid": "three alphas, twenty lambdas, five fixed folds and cold refit; not the 100-lambda API default",
            "privacy": "aggregate checks only; no sample IDs, cell IDs, phenotypes, expression or coefficients",
            "missing_report_policy": "not_reported does not distinguish pending from never run and is never a pass",
            "expected_same_threshold_grid": grid, "comparisons": records}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("comparisons", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = summarize(args.comparisons.glob("*.json"))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"comparison_reports": len(result["comparisons"]),
                      "status_counts": {status: sum(row["status"] == status for row in result["expected_same_threshold_grid"])
                                        for status in ("passed", "failed", "not_reported")}}))


if __name__ == "__main__":
    main()
