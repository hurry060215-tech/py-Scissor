"""Regenerate bounded numerical/memory/timing evidence; no R or downloads."""
import json
import platform
import time
import tracemalloc
from pathlib import Path
import numpy as np
import scipy
import pandas as pd
from scipy import sparse
from scissor_py import fit_path, cross_validate, normalize_graph
from scissor_py._cox import cox_score
ROOT = Path(__file__).resolve().parents[1]
report = {"date":"2026-10-09", "python":platform.python_version(), "platform":platform.system(),
          "numpy":np.__version__,"scipy":scipy.__version__,"pandas":pd.__version__,"explicit_r_fixtures":[]}
for family in ["gaussian","binomial","cox"]:
 for precision in ["default","tight"]:
  d=json.loads((ROOT/"reference"/"fixtures"/f"{family}_{precision}_explicit_cv.json").read_text()); i=d["inputs"]; e=d["output"]["result"]
  r=cross_validate(i["x"],i["y"],i["Omega"],family=family,alpha=i["alpha"],lambdas=i["lambda"],foldid=i["foldid"],thresh=i["thresh"])
  beta=np.vstack([r.path.intercepts,r.path.coefficients]) if family=="binomial" else r.path.coefficients
  report["explicit_r_fixtures"].append({"family":family,"precision":precision,"max_abs_coefficient_error":float(abs(beta-np.array(e["Beta"])).max()),"max_abs_cv_mean_error":float(abs(r.mean_loss-np.array([v['cvm'] for v in e['fit']])).max()),"max_abs_cv_se_error":float(abs(r.standard_error-np.array([v['cvse'] for v in e['fit']])).max()),"lambda_min_matches":bool(r.lambda_min==e['lambda.min'])})
n=100_000; graph=sparse.diags([np.ones(n-1),np.ones(n-1)],[-1,1],format='csr'); w=normalize_graph(graph)
report['sparse_graph']={"nodes":n,"nnz":w.nnz,"csr_bytes":w.data.nbytes+w.indices.nbytes+w.indptr.nbytes,"dense_equivalent_bytes":n*n*8}
rng=np.random.default_rng(712); n=6000; X=rng.normal(size=(n,3)); y=np.column_stack([np.arange(n),np.ones(n)])
tracemalloc.start(); cox_score(X,y,[.1,-.3,.2]); _,peak=tracemalloc.get_traced_memory();tracemalloc.stop()
report['cox_risk_memory']={"samples":n,"predictors":3,"python_numpy_tracemalloc_peak_bytes":peak,"dense_risk_matrix_bytes":n*n*8,"note":"tracemalloc is not total process RSS"}
rng=np.random.default_rng(521); n=60; p=80; X=rng.normal(size=(n,p))*.3; graph=sparse.diags([np.ones(p-1),np.ones(p-1)],[-1,1],format='csr')
yg=1.3*X[:,0]-.8*X[:,1]+rng.normal(size=n)*.2
yb=(yg>np.median(yg)).astype(float); yc=np.column_stack([np.repeat(np.arange(1,31),2),(np.arange(n)%5!=0).astype(float)])
report['bounded_path_timings']=[]
for family,y in [('gaussian',yg),('binomial',yb),('cox',yc)]:
 start=time.perf_counter(); r=fit_path(X,y,graph,family=family,alpha=.4,lambdas=np.geomspace(.1,.002,6),semantics='corrected',thresh=1e-7); seconds=time.perf_counter()-start
 report['bounded_path_timings'].append({"family":family,"samples":n,"cells":p,"lambdas":6,"seconds":seconds,"converged":bool(r.converged.all()),"max_kkt":float(np.max(r.metadata['kkt_residual']))})
(ROOT/'docs'/'validation_metrics.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
