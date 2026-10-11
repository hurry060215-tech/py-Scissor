"""Run with an installed wheel from outside the source tree and isolated Python."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import socket

def no_child_process(*args, **kwargs):
    raise AssertionError("runtime smoke must not spawn R or any child process")

subprocess.Popen = no_child_process

def no_network(*args, **kwargs):
    raise AssertionError("runtime smoke and quickstart must not access the network")

socket.getaddrinfo = no_network
socket.create_connection = no_network
socket.socket.connect = no_network
socket.socket.connect_ex = no_network
socket.socket.sendto = no_network
import numpy as np
from scipy import sparse
import scissor_py

location = Path(scissor_py.__file__).resolve()
assert "site-packages" in str(location), location
assert importlib.util.find_spec("rpy2") is None
assert not any(name in sys.modules for name in ("anndata", "scanpy", "omicverse", "torch"))
rng = np.random.default_rng(72)
X = rng.normal(size=(30,5))*.3
graph = sparse.diags([np.ones(4),np.ones(4)],[-1,1],format="csr")
outcomes = {
    "gaussian": X[:,0]-X[:,1]+rng.normal(size=30)*.1,
    "binomial": (np.arange(30)%2).astype(float),
    "cox": np.column_stack([np.repeat(np.arange(1,16),2),(np.arange(30)%5!=0).astype(float)]),
}
families = []
for family,y in outcomes.items():
    fit = scissor_py.fit_path(X,y,graph,family=family,alpha=.4,semantics="corrected",lambdas=[.1,.03])
    assert fit.converged.all(), (family,fit.metadata)
    assert np.isfinite(fit.coefficients).all()
    families.append(family)
print(json.dumps({"version":scissor_py.__version__,"python":sys.version.split()[0],
                  "import_from_site_packages":True,"rpy2_present":False,
                  "optional_heavy_modules_imported":False,
                  "child_processes_forbidden":True,"network_forbidden":True,
                  "families_passed":families}))
