"""Keep public examples, version metadata and required reference files honest."""
import hashlib
import json
from pathlib import Path
import re
import runpy
import importlib.util
import os

import pytest

import numpy as np
import scissor_py

ROOT = Path(__file__).resolve().parents[1]


def test_public_exports_and_version_metadata_agree():
    assert len(scissor_py.__all__) == len(set(scissor_py.__all__))
    assert all(hasattr(scissor_py, name) for name in scissor_py.__all__)
    match = re.search(r'^version = "([^"]+)"$', (ROOT / "pyproject.toml").read_text(), re.M)
    assert match.group(1) == scissor_py.__version__
    assert f'version: "{scissor_py.__version__}"' in (ROOT / "CITATION.cff").read_text()


def test_required_reference_manifest_is_complete_and_unchanged():
    manifest = json.loads((ROOT / "reference/manifest.json").read_text())
    assert manifest["fixture_count"] == 22
    for name, expected in manifest["files"].items():
        payload = (ROOT / "reference" / name).read_bytes()
        assert hashlib.sha256(payload).hexdigest() == expected["sha256"], name
        assert len(payload) == expected["bytes"], name


def test_readme_quickstart_is_exact_runnable_script():
    readme = (ROOT / "README.md").read_text()
    block = readme.split("<!-- quickstart:start -->", 1)[1].split("<!-- quickstart:end -->", 1)[0]
    code = block.split("```python\n", 1)[1].rsplit("```", 1)[0]
    assert code.strip() == (ROOT / "examples/quickstart.py").read_text().strip()


def test_quickstart_uses_public_api_and_retains_all_observations():
    namespace = runpy.run_path(str(ROOT / "examples/quickstart.py"))
    result = namespace["result"]
    assert isinstance(result, scissor_py.ScissorResult)
    assert result.prepared.similarity.shape == (24, 8)
    assert result.metadata["semantics"] == "corrected"
    assert result.metadata["positive_class"]["positive"] == "mutant"
    assert np.isfinite(result.coefficients).all()
    assert all(cv.valid.any() for cv in result.cv_results)


def contract_deriver():
    path = ROOT / "scripts/derive_luad_validation_contract.py"
    spec = importlib.util.spec_from_file_location("derive_contract", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.derive


def test_threshold_derivation_only_changes_explicit_metadata(tmp_path):
    probe = tmp_path / "symlink-probe"
    try:
        probe.symlink_to(tmp_path / "target")
        probe.unlink()
    except OSError as exc:
        if os.name == "nt" and getattr(exc, "winerror", None) == 1314:
            pytest.skip("Windows symlink privilege is unavailable for this development helper")
        raise
    source = tmp_path / "prepared"
    source.mkdir()
    manifest = {"families": {name: {"thresh": 1e-8, "y": [1, 2], "foldid": [0, 1],
                                   "paths": [{"alpha": .2, "lambdas": [.1, .01]}]}
                              for name in ("gaussian", "binomial", "cox")}}
    payload = json.dumps(manifest).encode()
    (source / "manifest.json").write_bytes(payload)
    (source / "design.bin").write_bytes(b"unchanged input")
    destination = tmp_path / "tight"
    result = contract_deriver()(source, destination, 1e-10)
    assert (source / "manifest.json").read_bytes() == payload
    assert (destination / "design.bin").is_symlink()
    assert (destination / "design.bin").read_bytes() == b"unchanged input"
    derived = json.loads((destination / "manifest.json").read_text())
    assert result["source_manifest_sha256"] == hashlib.sha256(payload).hexdigest()
    derived.pop("explicit_threshold_derivation")
    for config in derived["families"].values():
        assert config["thresh"] == 1e-10
        config["thresh"] = 1e-8
    assert derived == manifest


def test_threshold_derivation_refuses_bad_threshold_or_existing_destination(tmp_path):
    derive = contract_deriver()
    with pytest.raises(ValueError, match="finite and positive"):
        derive(tmp_path / "missing", tmp_path / "output", float("nan"))
    with pytest.raises(ValueError, match="new directory"):
        derive(tmp_path / "missing", tmp_path, 1e-10)
