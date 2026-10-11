"""Recreate hash-exact R X.bin files from complete checkpoint X.npy files.

This restores model-validation inputs only. It does not reconstruct the large
joint/normalized matrices, download data, or claim a new preprocessing run.
"""
import argparse
import hashlib
import json
import os
import tempfile
from pathlib import Path

import numpy as np


def verified(path, expected):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    if path.stat().st_size != expected["bytes"] or digest.hexdigest() != expected["sha256"]:
        raise ValueError(f"frozen input checksum mismatch: {path.name}")


def restore(prepared):
    manifest = json.loads((prepared/"manifest.json").read_text())
    result = {}
    for graph in ("graph.npz", "graph.mtx"):
        verified(prepared/graph, manifest["files"][graph])
    for cohort in ("exp1", "exp2"):
        source, target = prepared/cohort/"X.npy", prepared/cohort/"X.bin"
        verified(source, manifest["files"][f"{cohort}/X.npy"])
        expected = manifest["files"][f"{cohort}/X.bin"]
        if target.exists():
            verified(target, expected)
            result[cohort] = "existing binary checksum verified"
            continue
        X = np.load(source, allow_pickle=False)
        if list(X.shape) != manifest["cohorts"][cohort]["X_shape"] or not np.isfinite(X).all():
            raise ValueError(f"invalid frozen design: {cohort}")
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(dir=target.parent, prefix=".X.bin-", delete=False) as handle:
                temporary = Path(handle.name)
                np.asarray(X, dtype="<f8").ravel(order="F").tofile(handle)
                handle.flush()
                os.fsync(handle.fileno())
            verified(temporary, expected)
            # Explicit recovery tool; existing destination was already checked
            # above. Do not invoke this concurrently on the same destination.
            os.replace(temporary, target)
            result[cohort] = "recreated binary checksum verified"
        finally:
            if temporary is not None and temporary.exists():
                temporary.unlink()
    return {"scope": "complete frozen model designs and graph only", "status": "passed", "cohorts": result}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prepared", type=Path)
    args = parser.parse_args()
    print(json.dumps(restore(args.prepared)))
