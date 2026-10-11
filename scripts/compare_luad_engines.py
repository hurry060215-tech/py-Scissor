"""Audit complete saved NumPy/Numba Gaussian paths, CV, and cold refit."""
import argparse
import json
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("numpy_folder", type=Path)
    parser.add_argument("numba_folder", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--legacy-numpy-metadata", action="store_true",
                        help="Explicitly identify the recorded NumPy-only run launched before the engine option existed")
    args = parser.parse_args()
    a = json.loads((args.numpy_folder/"result.json").read_text())
    b = json.loads((args.numba_folder/"result.json").read_text())
    assert a["status"] == b["status"] == "passed"
    if "engine" not in a:
        assert args.legacy_numpy_metadata
        assert all("numba_version" not in path["metadata"] for path in a["paths"])
        assert "numba_version" not in a["refit_metadata"]
    else:
        assert a["engine"] == "numpy"
    assert b["engine"] == "numba"
    for key in ("manifest_sha256", "family", "semantics", "alpha", "shape", "full_data", "sampling"):
        assert a[key] == b[key], key
    report = {"family": a["family"], "alpha": a["alpha"], "manifest_sha256": a["manifest_sha256"],
              "absolute_tolerance": 1e-12, "relative_tolerance": 1e-10,
              "arrays": {}, "numpy_seconds": a["seconds"], "numba_seconds_including_compile": b["seconds"],
              "numba_compile_seconds": b["jit_cold_warmup_seconds"],
              "numpy_engine_metadata_present": "engine" in a,
              "numpy_engine_provenance": "Explicit NumPy-only command started before engine option was introduced" if "engine" not in a else "engine=numpy recorded in result",
              "numpy_peak_rss_bytes": a["process_peak_rss_bytes"], "numba_peak_rss_bytes": b["process_peak_rss_bytes"],
              "timing_scope": "observed single-thread full data, all 20 lambdas, all 5 folds, cold selected-lambda refit; concurrent task load may differ"}
    for file in sorted(args.numpy_folder.glob("*.npz")):
        left, right = np.load(file), np.load(args.numba_folder/file.name)
        assert left.files == right.files
        entries = {}
        for key in left.files:
            x, y = left[key], right[key]
            assert x.shape == y.shape
            exact = np.issubdtype(x.dtype, np.integer) or x.dtype == np.bool_ or key == "lambdas"
            passed = np.array_equal(x, y) if exact else np.allclose(x, y, atol=1e-12, rtol=1e-10)
            entries[key] = {"passed": bool(passed), "entries": x.size, "requires_exact_equality": exact,
                            "max_absolute_difference": float(np.abs(x.astype(float)-y.astype(float)).max(initial=0))}
        report["arrays"][file.name] = entries
    report["selection"] = {key: a[key] == b[key] for key in ("lambda_min", "lambda_1se", "positive_count", "negative_count", "cutoff_reached")}
    report["passed"] = all(v["passed"] for item in report["arrays"].values() for v in item.values()) and all(report["selection"].values())
    report["observed_total_speed_ratio"] = a["seconds"]/b["seconds"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(json.dumps({k: report[k] for k in ("passed", "numpy_seconds", "numba_seconds_including_compile", "numba_compile_seconds", "observed_total_speed_ratio")}))
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
