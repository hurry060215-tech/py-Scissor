"""Reject oracle code, external contributions and private data in core archives."""
import argparse
import json
from pathlib import Path, PurePosixPath
import tarfile
import zipfile

PUBLIC_DOCS = {
    'docs/ARCHITECTURE.md', 'docs/API.md', 'docs/REPRODUCING.md',
    'docs/ANNDATA.md', 'docs/SEMANTICS.md', 'docs/VALIDATION.md',
    'docs/VALIDATION_STATUS.md', 'docs/validation_metrics.json',
    'docs/luad/INPUT_AUDIT.md', 'docs/luad/download_manifest.json',
    'docs/luad/current_preprocessing_regression.json',
    'docs/luad/numba_full_cv_comparison.json', 'docs/luad/numerical_summary.json',
}
NON_PYTHON_CODE = {'.r', '.c', '.cc', '.cpp', '.cxx', '.h', '.hpp', '.so', '.dll', '.pyd', '.dylib'}
EXTERNAL_TREES = {'integration', 'skills', 'tutorials'}
PRIVATE_PARTS = {'validation-data', 'deliverables', '.git', '.venv', '.ov-venv',
                 '__pycache__', '.pytest_cache'}


def check_distribution(path):
    path = Path(path)
    wheel = path.suffix == '.whl'
    if wheel:
        with zipfile.ZipFile(path) as archive:
            names = [item.filename for item in archive.infolist() if not item.is_dir()]
    else:
        with tarfile.open(path, 'r:gz') as archive:
            members = archive.getmembers()
            if any(item.issym() or item.islnk() for item in members):
                raise ValueError('Source distribution must not contain links')
            names = [item.name for item in members if item.isfile()]
    paths = [PurePosixPath(name) for name in names]
    if not paths or any(p.is_absolute() or '..' in p.parts for p in paths):
        raise ValueError('Invalid archive paths')
    if not wheel:
        if len({p.parts[0] for p in paths}) != 1:
            raise ValueError('Source distribution must have one root')
        paths = [PurePosixPath(*p.parts[1:]) for p in paths]
    for item in paths:
        name = item.as_posix()
        if item.suffix.lower() in NON_PYTHON_CODE or item.name == 'Makevars':
            raise ValueError(f'Non-Python/oracle code in core distribution: {name}')
        if item.suffix == '.ipynb' or item.parts[0] in EXTERNAL_TREES:
            raise ValueError(f'External contribution/tutorial in core distribution: {name}')
        if name.startswith('reference/luad/'):
            raise ValueError(f'Development oracle belongs in the separate archive: {name}')
        if wheel:
            if item.parts[0] != 'scissor_py' and not item.parts[0].endswith('.dist-info'):
                raise ValueError(f'Unexpected runtime wheel member: {name}')
            if item.parts[0] == 'scissor_py' and (len(item.parts) != 2 or item.suffix != '.py'):
                raise ValueError(f'Runtime wheel must contain Python modules only: {name}')
        elif item.parts[0] == 'docs' and name not in PUBLIC_DOCS:
            raise ValueError(f'Private or unreviewed documentation in sdist: {name}')
        if any(part in PRIVATE_PARTS for part in item.parts):
            raise ValueError(f'Private/development member: {name}')
        if name == 'scripts/export_public_source.py':
            raise ValueError('Private source-export helper must not be distributed')
    return {'file': path.name, 'kind': 'wheel' if wheel else 'sdist',
            'files': len(paths), 'reference_code_files': 0, 'external_contribution_files': 0,
            'public_file_boundary': 'passed'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('artifacts', nargs='+', type=Path)
    args = parser.parse_args()
    print(json.dumps([check_distribution(path) for path in args.artifacts], indent=2))


if __name__ == '__main__':
    main()
