"""Build the self-contained synthetic integration tutorial (no downloads)."""
from pathlib import Path
import nbformat as nb
root=Path(__file__).resolve().parents[1]
cells=[]
def md(text): cells.append(nb.v4.new_markdown_cell(text))
def code(text): cells.append(nb.v4.new_code_cell(text))
md('''# Python Scissor: explicit inputs, three phenotypes, auditable outputs

This executable tutorial uses a small **synthetic** example to demonstrate the
API. It does not replace the separate full LUAD validation (4,102 cells,
471 survival samples / 498 TP53 samples, 20 lambdas, five fixed folds and three
alphas). See `docs/luad/` and the full-data reproduction scripts.

Install the supplied `scissor-py` wheel. For the optional OmicVerse cells, install
the local candidate against the pinned upstream source as described in
`integration/omicverse/README.md`. The backend import is `scissor_py`, not the
unrelated `pyscissor`. R is used only to develop validation references.

The core's default semantics are `official`. This tutorial explicitly requests
`corrected` because its connected cell graph is not safe under the original
Gaussian/Cox wrapper's indexing. Corrected results must not be described as
unchanged official Scissor results. No data, model weights or graph is downloaded.''')
code('''from pathlib import Path
import sys, json, tempfile
ROOT = Path.cwd()
if not (ROOT / "pyproject.toml").exists():
    ROOT = ROOT.parent
sys.path.insert(0, str(ROOT / "integration" / "omicverse"))
from offline_guard import install
network_attempts = install(allow_loopback=True)  # keep local kernel transport only

import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad
from scissor_py import prepare, scissor_from_anndata, result_summary, write_h5ad_atomic
print("Python", sys.version.split()[0])
''')
md('''## 1. Choose expression and graph explicitly

Bulk has genes × samples; AnnData has cells × genes. All identifiers must be
unique. Here counts are generated locally and the chosen `lognorm` layer is
explicit log1p expression for the example. Real preprocessing choices belong
upstream. The demonstration chain graph is not a Seurat or Scanpy graph.
The kernel jointly quantile-normalizes shared genes and computes sample × cell
Pearson similarities; it does not silently run HVG selection or build neighbors.''')
code('''rng = np.random.default_rng(20261010)
genes = pd.Index([f"gene_{i}" for i in range(80)])
cells = pd.Index([f"cell_{i}" for i in range(8)])
samples = pd.Index([f"bulk_{i}" for i in range(30)])
bulk = pd.DataFrame(np.log1p(rng.poisson(8, (80, 30))), index=genes, columns=samples)
counts = rng.poisson(8, (8, 80)).astype(float)
adata = ad.AnnData(sparse.csr_matrix(counts), obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))
adata.layers["lognorm"] = sparse.csr_matrix(np.log1p(counts))
adata.obsp["explicit_graph"] = sparse.diags([np.ones(7), np.ones(7)], [-1, 1], format="csr")
adata.obs["original_group"] = pd.Categorical(["A"]*4 + ["B"]*4)
expression = pd.DataFrame(adata.layers["lognorm"].toarray().T, index=genes, columns=cells)
prepared = prepare(bulk, expression, adata.obsp["explicit_graph"], graph_cell_ids=cells)
assert prepared.similarity.shape == (30, 8)
foldid = pd.Series(np.arange(30) % 3, index=samples).iloc[::-1]  # deliberately reordered IDs
signal = 2 * prepared.similarity[:, 0] - 1.5 * prepared.similarity[:, 3]
print("Shared genes / samples / cells:", len(prepared.gene_ids), *prepared.similarity.shape)
''')
md('''## 2. Continuous, binary and survival outcomes

Gaussian uses a numeric Series. For binary strings, `positive_label` is required;
Scissor+ means association toward that class. Cox uses time/status (1 event,
0 right-censoring), Breslow ties, and event-weighted full-minus-training CV
partial likelihood, not C-index. These examples are generated signals, not
clinical outcomes or causal evidence. Explicit fold IDs permit cross-language
reproduction; the same random seed alone does not.''')
code('''gaussian = pd.Series(2 + signal + rng.normal(0, .005, 30), index=samples).iloc[::-1]
binary = pd.Series(np.where(signal > np.median(signal), "mutant", "wildtype"), index=samples).iloc[::-1]
status = (np.arange(30) % 4 != 0).astype(int)
time = np.maximum(.1, np.round(np.exp(-3*signal) * rng.uniform(1, 8, 30), 1))
survival = pd.DataFrame({"time":time, "status":status}, index=samples).iloc[::-1]
assert set(binary) == {"mutant", "wildtype"}
assert all(survival.loc[foldid.index[foldid == k], "status"].sum() > 0 for k in range(3))
print("Binary counts:", binary.value_counts().to_dict(), "Survival events:", int(status.sum()))
''')
md('''## 3. Run the independent backend

This demonstration explicitly uses eight lambdas and three fixed folds to keep
it small. The defaults are 100 lambdas and up to ten folds; the separate full
LUAD audit uses 20 and five. No observations are sampled away. Both convergence
and the final cold lambda.min refit must succeed. The first alpha meeting the
strict selected_fraction < cutoff stops the search. An unmet cutoff is reported,
not silently fixed by truncating coefficients.''')
code('''outcomes = {"gaussian": gaussian, "binomial": binary, "cox": survival}
results = {}
for family, phenotype in outcomes.items():
    kwargs = {"positive_label":"mutant"} if family == "binomial" else {}
    result = scissor_from_anndata(
        adata, bulk, phenotype, layer="lognorm", graph_key="explicit_graph",
        family=family, semantics="corrected", alpha=[.2, .6], cutoff=.75,
        foldid=foldid, nlambda=8, rlambda=.02, thresh=1e-10,
        writeback=True, key_added=family, **kwargs,
    )
    results[family] = result
    assert all(cv.valid.any() for cv in result.cv_results)
    assert result.metadata["fold_source"] == "provided"
print(pd.DataFrame({family:{"alpha":r.alpha, "lambda_min":r.lambda_min,
      "selected":len(r.positive_cells)+len(r.negative_cells), "cutoff_reached":r.cutoff_reached}
      for family,r in results.items()}).T.to_string())
''')
md('''## 4. Optional OmicVerse class

The candidate contributes `ov.single.Scissor` and only two registration lines.
It does not rework OmicVerse's existing import system. Its backend dependency is
checked when `run()` is called. `layer=None` selects X explicitly; a missing graph
asks for upstream preprocessing. It delegates the same numerical work and saves
results only after successful validation and serialization.''')
code('''import omicverse as ov
model = ov.single.Scissor(
    adata, bulk=bulk, phenotype=binary, family="binomial", layer="lognorm",
    graph_key="explicit_graph", positive_label="mutant", semantics="corrected",
    key_added="TP53_demo",
)
ov_result = model.run(alpha=[.2,.6], cutoff=.75, random_state=123,
                      foldid=foldid, nlambda=8, rlambda=.02, thresh=1e-10)
pd.testing.assert_series_equal(ov_result.coefficients, results["binomial"].coefficients)
assert model.result_ is ov_result
print("OmicVerse:", ov.__version__, "coefficient match: exact")
''')
md('''## 5. Inspect and save without losing types

AnnData uses labels Scissor+/Scissor-/background and the coefficient column
`<key>_coef`. The backend dataclass retains the original `Background` label;
the display mapping is recorded. `uns[key]` contains params/source/versions and
JSON-encoded encoding/sharedgenes/CV/convergence. Coefficients use strict signs;
near-zero diagnostics do not change labels.

The summary is an audit record, not a full model-reload artifact: dense expression,
similarity and all-path coefficient matrices remain outside it. Explicit
nonfinite-diagnostic tags keep JSON standards compliant. The separate atomic
H5AD helper does not implicitly convert object-string metadata to categoricals.''')
code('''audit = result_summary(ov_result)
json_text = json.dumps(audit, allow_nan=False)
assert json.loads(adata.uns["TP53_demo"]["encoding"]) == {"positive":"mutant", "negative":"wildtype"}
assert "normalized_network_sha256" in audit["network"]
with tempfile.TemporaryDirectory() as directory:
    path = Path(directory) / "scissor-demo.h5ad"
    write_h5ad_atomic(adata, path)
    restored = ad.read_h5ad(path)
    pd.testing.assert_frame_equal(restored.obs, adata.obs)
    assert restored.uns["TP53_demo"]["result_json"] == adata.uns["TP53_demo"]["result_json"]
    try:
        write_h5ad_atomic(adata, path)
    except FileExistsError:
        pass
    else:
        raise AssertionError("Existing file must not be overwritten implicitly")
print("Standard JSON and H5AD round-trip passed")
''')
md('''## 6. Failure must leave prior results intact

Repeating a run with the same output key fails unless overwrite=True is explicit.
Invalid phenotype IDs, classes, layers, graphs, folds, budgets or convergence
settings also fail without half-written annotations. For view/backed objects,
make an explicit in-memory copy before requesting writeback. Read-only fitting
is available by omitting writeback.''')
code('''old_obs = adata.obs.copy(deep=True)
old_json = adata.uns["TP53_demo"]["result_json"]
try:
    model.run(alpha=.2, lambdas=[1], foldid=foldid)
except ValueError as error:
    assert "already exists" in str(error)
else:
    raise AssertionError("Expected output collision")
pd.testing.assert_frame_equal(adata.obs, old_obs)
assert adata.uns["TP53_demo"]["result_json"] == old_json
print("No partial overwrite; blocked external-network attempts:", len(network_attempts))
''')
md('''## Full-data reproducibility and limits

Use `scripts/download_luad_validation.py` / the source manifest and
`scripts/prepare_luad_validation.py` as documented in `docs/luad/INPUT_AUDIT.md`;
then preserve the resulting graph, common-gene ordering, outcomes, actual fold
IDs and lambda paths for both languages. See `reference/luad/README.md` for the
independent R oracle and its checkpoint verification. Check script filenames
and command help in the supplied release before running; downloads are explicit,
never part of `Scissor.run()`.

The original paper's graph and default 100-lambda workflow are not claimed to
have been reproduced by the bounded full-data audit. Default stopping tolerances
can match predictions/CV labels while failing tight coefficient comparisons on
correlated cells; retain those failures separately from tighter-threshold runs.
See the release's current validation report for passed, failed and pending
checks. This tutorial's synthetic success is not a substitute for that report.
''')
notebook=nb.v4.new_notebook(cells=cells, metadata={'kernelspec':{'display_name':'Scissor validation','language':'python','name':'scissor-validation'},'language_info':{'name':'python'}})
nb.write(notebook,root/'tutorials/scissor_walkthrough.ipynb')
print('Wrote',root/'tutorials/scissor_walkthrough.ipynb')
