"""Make a separate explicit-threshold development contract without copying data.

The original manifest and matrices are never modified. All non-manifest entries
are symbolic links to the existing prepared directory, which must remain present
and immutable. This is a development helper, never called by the library.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path


def derive(source, destination, threshold):
    source, destination = Path(source).resolve(), Path(destination).resolve()
    if not math.isfinite(threshold) or threshold <= 0:
        raise ValueError("threshold must be finite and positive")
    if destination.exists() or destination.is_relative_to(source):
        raise ValueError("destination must be a new directory outside the source")
    payload = (source / "manifest.json").read_bytes()
    manifest = json.loads(payload)
    if set(manifest["families"]) != {"gaussian", "binomial", "cox"}:
        raise ValueError("expected the three-family LUAD development manifest")
    previous = {name: config["thresh"] for name, config in manifest["families"].items()}
    for config in manifest["families"].values():
        config["thresh"] = threshold
    manifest["explicit_threshold_derivation"] = {
        "source_manifest_sha256": hashlib.sha256(payload).hexdigest(),
        "previous_thresholds": previous, "requested_threshold": threshold,
        "numerical_inputs_changed": False,
        "storage": "symbolic links to the immutable prepared input; not independent copies",
    }
    destination.mkdir()
    try:
        for entry in source.iterdir():
            if entry.name != "manifest.json":
                (destination / entry.name).symlink_to(entry, target_is_directory=entry.is_dir())
        output = json.dumps(manifest, indent=2, allow_nan=False) + "\n"
        (destination / "manifest.json").write_text(output)
    except Exception:
        # Only remove entries created by this invocation in its fresh directory.
        for entry in destination.iterdir():
            entry.unlink()
        destination.rmdir()
        raise
    return manifest["explicit_threshold_derivation"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--thresh", type=float, required=True)
    args = parser.parse_args()
    print(json.dumps(derive(args.source, args.destination, args.thresh), indent=2))


if __name__ == "__main__":
    main()
