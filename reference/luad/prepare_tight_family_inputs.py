"""GPL-3.0-or-later. Derive a separate explicitly tightened family contract."""
import argparse
import hashlib
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("family", choices=("gaussian", "binomial", "cox"))
parser.add_argument("destination", type=Path)
args = parser.parse_args()
root = Path(__file__).resolve().parents[2]
original = root / "validation-data/luad/prepared"
original_bytes = (original / "manifest.json").read_bytes()
original_sha = hashlib.sha256(original_bytes).hexdigest()
assert original_sha == "001f3cd90dd02953ed7379e01827699bc4da7a621ef055500289d4e908991eee"
manifest = json.loads(original_bytes)
assert manifest["families"][args.family]["thresh"] == 1e-8
manifest["families"][args.family]["thresh"] = 1e-10
manifest["derived_validation"] = {
    "original_manifest_sha256": original_sha,
    "only_parameter_change": f"families.{args.family}.thresh: 1e-8 -> 1e-10",
    "purpose": "additional tighter cross-language convergence audit; retain original default results",
    "all_design_graph_outcome_fold_lambda_files_unchanged": True,
}
destination = args.destination.resolve()
destination.mkdir(parents=True, exist_ok=True)
for name in manifest["files"]:
    source = original / name
    target = destination / name
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.is_symlink():
        assert target.resolve() == source.resolve()
    else:
        assert not target.exists(), target
        target.symlink_to(source)
derived_bytes = (json.dumps(manifest, indent=2) + "\n").encode()
target = destination / "manifest.json"
if target.exists():
    assert target.read_bytes() == derived_bytes, "Refusing to overwrite a different frozen contract"
else:
    target.write_bytes(derived_bytes)
derived_sha = hashlib.sha256(derived_bytes).hexdigest()
derivation = (json.dumps({**manifest["derived_validation"], "derived_manifest_sha256": derived_sha}, indent=2) + "\n").encode()
target = destination / "derivation.json"
if target.exists():
    assert target.read_bytes() == derivation
else:
    target.write_bytes(derivation)
print(json.dumps({"family": args.family, "destination": str(destination), "derived_manifest_sha256": derived_sha}))
