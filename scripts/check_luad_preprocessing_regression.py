"""Check every frozen LUAD element after a preprocessing implementation edit."""
import argparse
import gc
import hashlib
import json
import resource
import time
from pathlib import Path

import numpy as np

from scissor_py import quantile_normalize
from scissor_py._preprocessing import _pearson_columns


def compare_all(actual, expected, atol=1e-8, rtol=1e-6):
    largest = 0.
    mismatches = 0
    for start in range(0, actual.shape[1], 64):
        a, b = actual[:, start:start+64], expected[:, start:start+64]
        delta = np.abs(a-b)
        largest = max(largest, float(delta.max(initial=0)))
        mismatches += int(np.count_nonzero(delta > atol+rtol*np.abs(b)))
    return {"elements": actual.size, "max_absolute_error": largest, "outside_tolerance": mismatches,
            "atol": atol, "rtol": rtol}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prepared", type=Path)
    parser.add_argument("r_output", type=Path)
    parser.add_argument("report", type=Path)
    args = parser.parse_args()
    manifest = json.loads((args.prepared/"manifest.json").read_text())
    report = {"source_sha256": hashlib.sha256((Path(__file__).resolve().parents[1]/"src/scissor_py/_preprocessing.py").read_bytes()).hexdigest(),
              "sampling": "none", "cohorts": {}}
    start = time.monotonic()
    for name, spec in manifest["cohorts"].items():
        shape = tuple(spec["joint_shape"])
        n, p = spec["X_shape"]
        joint = np.memmap(args.prepared/name/"joint.bin", dtype="<f8", mode="r", shape=shape, order="F")
        frozen = np.memmap(args.prepared/name/"normalized.bin", dtype="<f8", mode="r", shape=shape, order="F")
        actual = quantile_normalize(joint)
        norm = compare_all(actual, frozen)
        centered = _pearson_columns(actual[:, :n])
        x = np.empty((n, p))
        for j in range(0, p, 64):
            x[:, j:j+64] = centered.T@_pearson_columns(actual[:, n+j:n+j+64])
        np.clip(x, -1, 1, out=x)
        r_x = np.memmap(args.r_output/f"{name}-R-X.bin", dtype="<f8", mode="r", shape=(n, p), order="F")
        info = {"normalized_vs_frozen": norm,
                "pearson_vs_frozen": compare_all(x, np.load(args.prepared/name/"X.npy")),
                "pearson_vs_actual_R": compare_all(x, r_x)}
        assert all(v["outside_tolerance"] == 0 for v in info.values())
        report["cohorts"][name] = info
        print(name, json.dumps(info), flush=True)
        del actual, joint, frozen, centered, x, r_x
        gc.collect()
    report.update(status="passed", seconds=time.monotonic()-start,
                  process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    args.report.write_text(json.dumps(report, indent=2)+"\n")


if __name__ == "__main__":
    main()
