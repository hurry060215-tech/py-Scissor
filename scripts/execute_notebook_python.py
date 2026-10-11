"""Run plain-Python notebook cells without a socket/Jupyter kernel.

All supplied tutorial cells are plain Python (no magics, display protocol or
shell cells). Execute the exact saved source, including all assertions, in one
shared namespace, capture stdout/stderr, and mark the execution method explicitly.
This is not evidence that a Jupyter kernel or rich display was tested. It needs
nbformat and the tutorial's actual dependencies in the current environment.
"""
import argparse
import contextlib
import io
import json
from pathlib import Path
import time
import nbformat

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('notebook',type=Path)
parser.add_argument('--output',type=Path)
parser.add_argument('--report',type=Path)
args=parser.parse_args()
project_root = Path(__file__).resolve().parents[1]
notebook=nbformat.read(args.notebook,as_version=4)
namespace={'__name__':'__main__'}
start=time.monotonic();count=0
for cell in notebook.cells:
    if cell.cell_type != 'code': continue
    count+=1
    stdout,stderr=io.StringIO(),io.StringIO()
    # compile rejects unsupported magics rather than silently skipping them.
    with contextlib.redirect_stdout(stdout),contextlib.redirect_stderr(stderr):
        exec(compile(cell.source,f'{args.notebook.name}:cell{count}','exec'),namespace)
    cell.execution_count=count
    cell.outputs=[]
    for name,stream in [('stdout',stdout),('stderr',stderr)]:
        if stream.getvalue():
            text=stream.getvalue().replace(str(project_root), '$PROJECT')
            cell.outputs.append(nbformat.v4.new_output('stream',name=name,text=text))
try:
    notebook_name=str(args.notebook.resolve().relative_to(project_root))
except ValueError:
    notebook_name=args.notebook.name
report={'passed':True,'notebook':notebook_name,'code_cells':count,'seconds':time.monotonic()-start,
        'execution_method':'plain Python shared namespace, exact notebook source and assertions; no Jupyter kernel',
        'jupyter_kernel_verified':False,'network':'tutorial installs process-local external HTTP/DNS/socket guard',
        'scope':'synthetic three-family API, actual installed OmicVerse candidate, JSON/H5AD; no full LUAD refitting',
        'output_normalization':'Local project prefix in stdout/stderr replaced by $PROJECT; code and numerical outputs unchanged'}
notebook.metadata['scissor_execution']=report
nbformat.validate(notebook)
nbformat.write(notebook,args.output or args.notebook)
if args.report:
    args.report.parent.mkdir(parents=True,exist_ok=True)
    args.report.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
