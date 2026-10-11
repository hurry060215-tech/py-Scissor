"""Check built artifacts for the project's explicit public-file boundary.

Run after a build, including builds directly from a development checkout. This
checks archive members without extracting them, and notebook content for private
execution paths. It is not a general secret scanner: the source exporter and
review also check source text before publishing.
"""
import argparse
import json
from pathlib import Path, PurePosixPath
import tarfile
import zipfile

PUBLIC_DOCS = {
    "docs/API.md", "docs/REPRODUCING.md",
    "docs/ANNDATA.md", "docs/SEMANTICS.md", "docs/VALIDATION.md",
    "docs/VALIDATION_STATUS.md", "docs/validation_metrics.json",
    "docs/notebook_execution.json", "docs/luad/INPUT_AUDIT.md",
    "docs/luad/download_manifest.json", "docs/luad/current_preprocessing_regression.json",
    "docs/luad/numba_full_cv_comparison.json", "docs/luad/numerical_summary.json",
}


def check_distribution(path):
    path = Path(path)
    wheel = path.suffix == ".whl"
    notebooks = {}
    if wheel:
        with zipfile.ZipFile(path) as archive:
            names = [item.filename for item in archive.infolist() if not item.is_dir()]
            notebooks = {name: archive.read(name) for name in names if name.endswith(".ipynb")}
    else:
        with tarfile.open(path, "r:gz") as archive:
            members = archive.getmembers()
            if any(item.issym() or item.islnk() for item in members):
                raise ValueError("Source distribution must not contain links")
            names = [item.name for item in members if item.isfile()]
            notebooks = {name: archive.extractfile(name).read() for name in names if name.endswith(".ipynb")}
    paths = [PurePosixPath(name) for name in names]
    if not paths or any(p.is_absolute() or ".." in p.parts for p in paths):
        raise ValueError("Invalid archive paths")
    if not wheel:
        if len({p.parts[0] for p in paths}) != 1:
            raise ValueError("Source distribution must have one root")
        paths = [PurePosixPath(*p.parts[1:]) for p in paths]
    for item in paths:
        name = item.as_posix()
        if wheel:
            if item.parts[0] != "scissor_py" and not item.parts[0].endswith(".dist-info"):
                raise ValueError(f"Unexpected runtime wheel member: {name}")
        elif item.parts[0] == "docs" and name not in PUBLIC_DOCS:
            raise ValueError(f"Private or unreviewed documentation in sdist: {name}")
        if any(part in {"validation-data", "deliverables", ".git", ".venv", ".ov-venv",
                        "__pycache__", ".pytest_cache"} for part in item.parts):
            raise ValueError(f"Private/development member: {name}")
        if name == "scripts/export_public_source.py":
            raise ValueError("Private source-export helper must not be distributed")
    for name, payload in notebooks.items():
        # Stored warnings/metadata are distributable content too, not just code.
        text = json.dumps(json.loads(payload), ensure_ascii=False)
        if any(prefix in text for prefix in ("/workspace/", "/home/", "Sentinel" + "_", "libfile" + "_")):
            raise ValueError(f"Private execution metadata in notebook: {name}")
    return {"file": path.name, "kind": "wheel" if wheel else "sdist",
            "files": len(paths), "notebooks_content_checked": len(notebooks),
            "public_file_boundary": "passed"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifacts", nargs="+", type=Path)
    args = parser.parse_args()
    print(json.dumps([check_distribution(path) for path in args.artifacts], indent=2))


if __name__ == "__main__":
    main()
