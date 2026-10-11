"""Execute a notebook with the current Python and a private Unix-IPC kernel.

Requires nbclient, nbformat and ipykernel in this environment. Does not register
an enduring kernel, download data, or modify other Python environments. Run from
the project root so tutorial-relative resources resolve. Notebook cells retain
all assertions; cell errors fail the command and never produce a success report.
"""
import argparse
import asyncio
import json
import os
from pathlib import Path
import sys
import tempfile
import time
import nbformat
from nbclient import NotebookClient
from jupyter_client import AsyncKernelManager
from jupyter_client.kernelspec import KernelSpecManager

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('notebook', type=Path)
parser.add_argument('--output', type=Path)
parser.add_argument('--report', type=Path)
args = parser.parse_args()
notebook = nbformat.read(args.notebook, as_version=4)
start = time.monotonic()
with tempfile.TemporaryDirectory(prefix='scissor-notebook-') as directory:
    directory = Path(directory)
    spec = directory/'kernels'/'scissor-validation'; spec.mkdir(parents=True)
    (spec/'kernel.json').write_text(json.dumps({'argv':[sys.executable, '-m', 'ipykernel_launcher', '-f', '{connection_file}'],
        'display_name':'Scissor validation', 'language':'python'}))
    manager = KernelSpecManager(kernel_dirs=[str(directory/'kernels')])
    kernel = AsyncKernelManager(kernel_name='scissor-validation', kernel_spec_manager=manager,
                                transport='ipc', ip=str(directory/'ipc'))
    client = NotebookClient(notebook, km=kernel, timeout=300,
                            resources={'metadata':{'path':str(Path.cwd())}})
    try:
        client.execute()
    finally:
        if kernel.has_kernel:
            asyncio.run(kernel.shutdown_kernel(now=True))
nbformat.write(notebook, args.output or args.notebook)
report = {'passed':True, 'notebook':str(args.notebook), 'seconds':time.monotonic()-start,
          'code_cells':sum(cell.cell_type == 'code' for cell in notebook.cells),
          'execution_counts':[cell.execution_count for cell in notebook.cells if cell.cell_type == 'code'],
          'transport':'local Unix IPC', 'scope':'all notebook code cells, actual kernel execution'}
if args.report:
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
