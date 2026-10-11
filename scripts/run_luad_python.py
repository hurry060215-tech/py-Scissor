"""Run the public Python CV/refit API on frozen full-LUAD inputs, no downloads.

The development-only spy preserves every actual internal full/fold path for
independent R comparison; it does not change the API's numerical computation.
"""
import argparse
import hashlib
import json
import resource
import time
from pathlib import Path
from unittest.mock import patch

import numpy as np
from scipy import sparse

import scissor_py._api as api
from scissor_py import ScissorError, cross_validate, fit_path, normalize_graph
from scissor_py._cox import cox_log_likelihood
from scissor_py._solvers import network_quadratic


def plain(value):
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    if isinstance(value, np.ndarray):
        return plain(value.tolist())
    if isinstance(value, np.generic):
        return plain(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return str(value)
    return value


def write_json(path, record):
    path.write_text(json.dumps(plain(record), indent=2, allow_nan=False)+"\n")


def objectives(X, y, W, path):
    out = []
    for j, lam in enumerate(path.lambdas):
        beta = path.coefficients[:, j]
        eta = X@beta + path.intercepts[j]
        if path.family == "gaussian":
            loss = np.mean((eta-y)**2)/2
        elif path.family == "binomial":
            loss = np.mean(np.logaddexp(0, eta)-y*eta)
        else:
            loss = -cox_log_likelihood(X, y, beta)/len(y)
        out.append(loss+lam*path.alpha*np.abs(beta).sum()+lam*(1-path.alpha)*network_quadratic(beta, W)/2)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prepared", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--family", choices=["gaussian", "binomial", "cox"], required=True)
    parser.add_argument("--semantics", choices=["official", "corrected"], required=True)
    parser.add_argument("--alpha", type=float, required=True)
    parser.add_argument("--engine", choices=["numpy", "numba"], default="numpy")
    parser.add_argument("--thresh", type=float, help="Explicit additional convergence experiment; never changes frozen inputs")
    args = parser.parse_args()
    manifest_bytes = (args.prepared / "manifest.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    settings = manifest["families"][args.family]
    default_threshold = 1e-7 if args.semantics == "official" else settings["thresh"]
    threshold = default_threshold if args.thresh is None else args.thresh
    if not np.isfinite(threshold) or threshold <= 0:
        parser.error("--thresh must be finite and positive")
    path_config = next(x for x in settings["paths"] if x["alpha"] == args.alpha)
    X = np.load(args.prepared / settings["cohort"] / "X.npy")
    y = np.asarray(settings["y"])
    graph = sparse.load_npz(args.prepared / "graph.npz")
    W = normalize_graph(graph)
    folder = args.output / f"{args.family}_{args.semantics}_alpha_{args.alpha:g}"
    folder.mkdir(parents=True, exist_ok=False)
    report = {"status": "running", "manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
              "family": args.family, "semantics": args.semantics, "alpha": args.alpha, "engine": args.engine,
              "shape": list(X.shape), "full_data": True, "sampling": "none", "paths": [],
              "thresh": threshold, "default_validation_thresh": default_threshold,
              "additional_convergence_experiment": args.thresh is not None}
    start = time.monotonic()
    original = api._fit_path
    folds = np.asarray(settings["foldid"])
    levels = np.unique(folds)
    count = 0

    def capture(*pos, **kw):
        nonlocal count
        before = time.monotonic()
        result = original(*pos, **kw)
        label = "full" if count == 0 else f"fold_{levels[count-1]}"
        count += 1
        np.savez_compressed(folder / f"{label}.npz", coefficients=result.coefficients,
                            intercepts=result.intercepts, lambdas=result.lambdas,
                            converged=result.converged, iterations=result.iterations)
        item = {"label": label, "seconds": time.monotonic()-before,
                "converged": result.converged, "iterations": result.iterations,
                "metadata": result.metadata, "objective": objectives(pos[0], pos[1], W, result)}
        report["paths"].append(item)
        write_json(folder / "result.json", report)
        print(json.dumps(plain({k: item[k] for k in ("label", "seconds", "converged", "iterations")})), flush=True)
        return result

    try:
        if args.engine == "numba":
            before_compile = time.monotonic()
            fit_path(np.zeros((2, 2)), np.array([0., 1.]), sparse.csr_matrix((2, 2)),
                     family="gaussian", semantics="corrected", engine="numba", alpha=.5, lambdas=[1.])
            report["jit_cold_warmup_seconds"] = time.monotonic()-before_compile
            report["jit_warmup_note"] = "2x2 dummy compilation only; every scored validation path still uses complete frozen LUAD inputs"
        if args.semantics == "official" and args.family in ("gaussian", "cox"):
            try:
                fit_path(X, y, graph, family=args.family, semantics="official", alpha=args.alpha,
                         lambdas=path_config["lambdas"])
            except ScissorError as error:
                report.update(status="expected_unsafe_graph_rejection", error=str(error))
                return
            raise AssertionError("Unsafe connected official Gaussian/Cox graph was not rejected")
        with patch.object(api, "_fit_path", side_effect=capture):
            cv = cross_validate(X, y, graph, family=args.family, semantics=args.semantics,
                                alpha=args.alpha, lambdas=path_config["lambdas"], foldid=folds,
                                thresh=threshold, maxit=settings["maxit"], engine=args.engine)
        np.savez_compressed(folder / "cv.npz", mean_loss=cv.mean_loss, standard_error=cv.standard_error,
                            fold_loss=cv.fold_loss, valid=cv.valid, foldid=cv.foldid)
        before = time.monotonic()
        refit = fit_path(X, y, graph, family=args.family, semantics=args.semantics,
                         alpha=args.alpha, lambdas=[cv.lambda_min], thresh=threshold, maxit=settings["maxit"], engine=args.engine)
        np.savez_compressed(folder / "refit.npz", coefficients=refit.coefficients,
                            intercepts=refit.intercepts, converged=refit.converged, lambdas=refit.lambdas)
        beta = refit.coefficients[:, 0]
        report.update(status="passed" if bool(refit.converged.all()) else "refit_not_converged",
                      lambda_min=cv.lambda_min, lambda_1se=cv.lambda_1se, cv_metadata=cv.metadata,
                      refit_metadata=refit.metadata, refit_seconds=time.monotonic()-before,
                      refit_objective=objectives(X, y, W, refit), refit_iterations=refit.iterations,
                      positive_count=int(np.count_nonzero(beta > 0)), negative_count=int(np.count_nonzero(beta < 0)),
                      near_zero_nonzero_count=int(np.count_nonzero((np.abs(beta) < 1e-5) & (beta != 0))),
                      selected_fraction=float(np.count_nonzero(beta)/len(beta)),
                      cutoff_reached=bool(np.count_nonzero(beta)/len(beta) < settings["cutoff"]))
    except Exception as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report.update(seconds=time.monotonic()-start,
                      process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
        write_json(folder / "result.json", report)
        print(json.dumps(plain({k: v for k, v in report.items() if k not in ("paths", "cv_metadata", "refit_metadata")})), flush=True)


if __name__ == "__main__":
    main()
