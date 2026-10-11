"""The distribution checker rejects private and newly unreviewed report paths."""
import importlib.util
import io
from pathlib import Path
import tarfile
import zipfile

import pytest


def checker():
    path = Path(__file__).parents[1] / "scripts/check_distribution.py"
    spec = importlib.util.spec_from_file_location("distribution_boundary", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.check_distribution


def sdist(tmp_path, name, payload=b"{}"):
    path = tmp_path / "test.tar.gz"
    with tarfile.open(path, "w:gz") as archive:
        member = tarfile.TarInfo("scissor_py-test/" + name)
        member.size = len(payload)
        archive.addfile(member, io.BytesIO(payload))
    return path


@pytest.mark.parametrize("name", ["docs/luad/anndata-cox.json", "docs/luad/input_manifest.json",
                                  "docs/luad/comparisons/real.json", "docs/new-unreviewed-report.json",
                                  "validation-data/output.json", "scripts/export_public_source.py"])
def test_private_or_unreviewed_sdist_member_rejected(tmp_path, name):
    with pytest.raises(ValueError):
        checker()(sdist(tmp_path, name))


def test_public_aggregate_allowed(tmp_path):
    assert checker()(sdist(tmp_path, "docs/luad/numerical_summary.json"))["public_file_boundary"] == "passed"


def test_wheel_rejects_bundled_reference_runtime(tmp_path):
    path = tmp_path / "test.whl"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("reference/oracle.R", "")
    with pytest.raises(ValueError, match="Unexpected runtime"):
        checker()(path)


def test_wheel_package_and_license_allowed(tmp_path):
    path = tmp_path / "test.whl"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("scissor_py/__init__.py", "")
        archive.writestr("scissor_py-test.dist-info/licenses/LICENSE", "GPL")
    assert checker()(path)["public_file_boundary"] == "passed"


def test_private_notebook_warning_output_rejected(tmp_path):
    payload = b'{"cells":[{"outputs":[{"text":"/workspace/private/environment/warning"}]}]}'
    with pytest.raises(ValueError, match="Private execution metadata"):
        checker()(sdist(tmp_path, "tutorials/example.ipynb", payload))


def test_normalized_notebook_warning_output_allowed(tmp_path):
    payload = b'{"cells":[{"outputs":[{"text":"$PROJECT/environment/warning"}]}]}'
    result = checker()(sdist(tmp_path, "tutorials/example.ipynb", payload))
    assert result["notebooks_content_checked"] == 1
    assert result["public_file_boundary"] == "passed"
