"""Derive an explicit tighter Gaussian audit without changing frozen inputs."""
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parents[2]
original = root / "validation-data/luad/prepared"
destination = root / "validation-data/luad/r-oracle/tight-inputs"
original_bytes = (original / "manifest.json").read_bytes()
original_sha = hashlib.sha256(original_bytes).hexdigest()
assert original_sha == "001f3cd90dd02953ed7379e01827699bc4da7a621ef055500289d4e908991eee"
manifest = json.loads(original_bytes)
assert manifest["families"]["gaussian"]["thresh"] == 1e-8
manifest["families"]["gaussian"]["thresh"] = 1e-10
manifest["derived_validation"] = {
    "original_manifest_sha256": original_sha,
    "only_parameter_change": "families.gaussian.thresh: 1e-8 -> 1e-10",
    "purpose": "additional tighter cross-language convergence audit; retain original default results",
    "all_design_graph_outcome_fold_lambda_files_unchanged": True,
}
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
(destination / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
derived_sha = hashlib.sha256((destination / "manifest.json").read_bytes()).hexdigest()
(destination / "derivation.json").write_text(json.dumps({
    **manifest["derived_validation"], "derived_manifest_sha256": derived_sha,
}, indent=2) + "\n")
print(json.dumps({"destination": str(destination), "derived_manifest_sha256": derived_sha}))
