"""Hash every frozen exchange file; never sample file contents."""
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parents[2]
prepared = root / "validation-data/luad/prepared"
manifest_path = prepared / "manifest.json"
manifest_sha = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
expected_manifest = "001f3cd90dd02953ed7379e01827699bc4da7a621ef055500289d4e908991eee"
assert manifest_sha == expected_manifest, (manifest_sha, expected_manifest)
manifest = json.loads(manifest_path.read_text())
results = {}
for name, expected in manifest["files"].items():
    path = prepared / name
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(8 * 1024 * 1024):
            digest.update(block)
    found = {"sha256": digest.hexdigest(), "bytes": path.stat().st_size}
    assert found == expected, (name, found, expected)
    results[name] = found
output = root / "validation-data/luad/r-oracle/input-integrity.json"
output.write_text(json.dumps({"manifest_sha256": manifest_sha, "all_files_verified": True,
                              "files": results}, indent=2) + "\n")
print(f"Verified manifest and {len(results)} complete exchange files.")
