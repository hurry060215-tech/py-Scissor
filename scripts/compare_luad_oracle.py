"""Compare all saved full-LUAD paths, folds, CV and cold refits against R.

No fitting, downloads, selection of cells or tolerance adaptation occurs here.
The report preserves failures and all near-zero sign differences separately.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from scissor_py import normalize_graph
from scissor_py._cox import cox_log_likelihood
from scissor_py._solvers import network_quadratic


TOLERANCES = {
    "coefficients": {"atol": 1e-5, "rtol": 1e-4},
    "intercepts": {"atol": 1e-5, "rtol": 1e-4},
    "linear_predictions": {"atol": 1e-5, "rtol": 1e-4},
    "objective": {"atol": 1e-8, "rtol": 1e-6},
    "cv_loss_and_se": {"atol": 1e-6, "rtol": 1e-4},
}
STABLE_COEFFICIENT = 1e-5


def compare(a, b, kind):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape:
        return {"passed": False, "left_shape": list(a.shape), "right_shape": list(b.shape)}
    finite = np.isfinite(a) & np.isfinite(b)
    error = np.abs(a-b)
    tol = TOLERANCES[kind]
    ok = finite & (error <= tol["atol"] + tol["rtol"]*np.abs(b))
    return {"passed": bool(ok.all()), "entries": a.size,
            "max_absolute_error": float(error.max(initial=0)) if finite.all() else None,
            "outside_tolerance": int((~ok).sum())}


def objective(X, y, W, beta, intercept, lambdas, family, alpha):
    result = []
    for k, lam in enumerate(lambdas):
        b = beta[:, k]
        eta = X@b+intercept[k]
        if family == "gaussian":
            loss = np.mean((eta-y)**2)/2
        elif family == "binomial":
            loss = np.mean(np.logaddexp(0, eta)-y*eta)
        else:
            loss = -cox_log_likelihood(X, y, b)/len(y)
        result.append(loss + lam*alpha*np.abs(b).sum()
                      + lam*(1-alpha)*network_quadratic(b, W)/2)
    return np.array(result)


def run(prepared, python_folder, r_folder):
    manifest_bytes = (prepared/"manifest.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    py = json.loads((python_folder/"result.json").read_text())
    r = json.loads((r_folder/"result.json").read_text())
    family, alpha = py["family"], py["alpha"]
    config = manifest["families"][family]
    cohort = manifest["cohorts"][config["cohort"]]
    X = np.load(prepared/config["cohort"]/"X.npy")
    y = np.asarray(config["y"])
    foldid = np.asarray(config["foldid"])
    levels = np.unique(foldid)
    W = normalize_graph(sparse.load_npz(prepared/"graph.npz"))
    expected_lambdas = next(p["lambdas"] for p in config["paths"] if p["alpha"] == alpha)
    metadata = r["metadata"]
    r_semantics = "official" if metadata["semantics"] == "unchanged APML1 binomial" else metadata["semantics"]
    contract = {
        "python_completed": py["status"] == "passed",
        "manifest_hash_matches": py["manifest_sha256"] == hashlib.sha256(manifest_bytes).hexdigest(),
        "family_matches": metadata["family"] == family,
        "alpha_matches": metadata["alpha"] == alpha,
        "semantics_matches": r_semantics == py["semantics"],
        "foldid_matches": np.array_equal(metadata["foldid"], foldid),
        "lambda_path_matches": np.array_equal(metadata["lambdas"], expected_lambdas),
        "sample_count_matches": metadata["n_samples"] == X.shape[0],
        "cell_count_matches": metadata["n_cells"] == X.shape[1],
        "r_no_sampling": metadata["no_sampling"] is True,
        "python_no_sampling": py["sampling"] == "none",
    }
    for ids in ("sample_ids", "cell_ids"):
        if ids in metadata:
            contract[f"{ids}_match"] = metadata[ids] == cohort[ids]
    report = {"family": family, "semantics": py["semantics"], "alpha": alpha,
              "python_engine": py.get("engine", "numpy"), "manifest_sha256": py["manifest_sha256"],
              "tolerances": TOLERANCES, "stable_coefficient_threshold": STABLE_COEFFICIENT,
              "python_threshold": py.get("thresh", 1e-7 if py["semantics"] == "official" else config["thresh"]),
              "r_threshold": metadata["thresh"], "contract": contract,
              "r_source_modified": metadata.get("upstream_source_patched"),
              "r_semantics_description": metadata["semantics"], "paths": {}}
    report["same_requested_tolerance"] = report["python_threshold"] == report["r_threshold"]
    report["comparison_kind"] = ("same requested tolerance" if report["same_requested_tolerance"]
                                 else "additional convergence experiment with explicitly different requested tolerances")
    records = [("full", r["full"], np.ones(len(y), dtype=bool), expected_lambdas)]
    records.extend((f"fold_{level}", rec, foldid != level, expected_lambdas)
                   for level, rec in zip(levels, r["folds"], strict=True))
    records.append(("refit", r["refit"], np.ones(len(y), dtype=bool), [r["lambda_min"]]))
    for label, rec, train, lambdas in records:
        p = np.load(python_folder/f"{label}.npz")
        rb = np.fromfile(r_folder/rec["coefficient_file"], dtype="<f8").reshape(rec["coefficient_shape"], order="F")
        pb = p["coefficients"]
        ri, pi = np.atleast_1d(rec["intercepts"]), p["intercepts"]
        sign_diff = np.sign(pb) != np.sign(rb)
        stable_diff = sign_diff & (np.maximum(np.abs(pb), np.abs(rb)) >= STABLE_COEFFICIENT)
        positions = np.argwhere(sign_diff)
        rp = X@rb+ri
        pp = X@pb+pi
        item = {
            "coefficients": compare(pb, rb, "coefficients"),
            "intercepts": compare(pi, ri, "intercepts"),
            "linear_predictions_all_samples": compare(pp, rp, "linear_predictions"),
            "objective": compare(objective(X[train], y[train], W, pb, pi, lambdas, family, alpha),
                                 objective(X[train], y[train], W, rb, ri, lambdas, family, alpha), "objective"),
            "python_all_converged": bool(p["converged"].all()),
            "r_all_converged": bool(np.all(rec["converged"])) if "converged" in rec else bool(np.all(np.asarray(rec["upstream_flag"]) == 0)),
            "python_lambda_path_matches": bool(np.array_equal(p["lambdas"], lambdas)),
            "strict_sign_differences": int(sign_diff.sum()),
            "stable_sign_differences": int(stable_diff.sum()),
            "near_zero_sign_examples": [{"cell_id": cohort["cell_ids"][i], "lambda_index_0based": int(j),
                                          "lambda": float(lambdas[j]), "python": float(pb[i, j]), "R": float(rb[i, j])}
                                         for i, j in positions[:20]],
        }
        item["passed"] = all(item[k]["passed"] for k in ("coefficients", "intercepts", "linear_predictions_all_samples", "objective")) and all(item[k] for k in ("python_all_converged", "r_all_converged", "python_lambda_path_matches")) and not stable_diff.any()
        report["paths"][label] = item
    cv = np.load(python_folder/"cv.npz")
    report["cv"] = {k: compare(cv[pkey], r[rkey], "cv_loss_and_se") for k, pkey, rkey in
                    [("mean_loss", "mean_loss", "cvm"), ("standard_error", "standard_error", "cvse"),
                     ("fold_losses", "fold_loss", "fold_losses")]}
    report["cv"].update(lambda_min_exact=py["lambda_min"] == r["lambda_min"],
                         valid_matches=bool(np.array_equal(cv["valid"], r.get("valid", np.ones(len(expected_lambdas), dtype=bool)))))
    if "lambda_1se" in r:
        report["cv"]["lambda_1se_exact"] = py["lambda_1se"] == r["lambda_1se"]
    else:
        best = int(np.argmin(r["cvm"]))
        one = int(np.flatnonzero(np.asarray(r["cvm"]) <= r["cvm"][best]+r["cvse"][best])[0])
        report["cv"]["lambda_1se_exact"] = py["lambda_1se"] == expected_lambdas[one]
    report["passed"] = all(contract.values()) and all(x["passed"] for x in report["paths"].values()) and all(
        x["passed"] if isinstance(x, dict) else x for x in report["cv"].values())
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prepared", type=Path)
    parser.add_argument("python_folder", type=Path)
    parser.add_argument("r_folder", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = run(args.prepared, args.python_folder, args.r_folder)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps({k: result[k] for k in ("family", "semantics", "alpha", "passed")}))
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
