import numpy as np
import pandas as pd
from scipy import sparse
from scissor_py import scissor

rng = np.random.default_rng(123)
genes = [f"gene_{i}" for i in range(60)]
samples = [f"sample_{i}" for i in range(24)]
cell_ids = [f"cell_{i}" for i in range(8)]
bulk = pd.DataFrame(rng.gamma(2, 2, (60, 24)), index=genes, columns=samples)
cells = pd.DataFrame(rng.gamma(2, 2, (60, 8)), index=genes, columns=cell_ids)
phenotype = pd.Series(["wildtype", "mutant"] * 12, index=samples)
graph = sparse.diags([np.ones(7), np.ones(7)], [-1, 1], format="csr")
foldid = pd.Series(np.arange(24) % 3, index=samples)

result = scissor(
    bulk, cells, phenotype, graph, graph_cell_ids=cell_ids,
    family="binomial", positive_label="mutant", semantics="corrected",
    alpha=[.2, .6], lambdas=[.10, .04, .015], foldid=foldid,
    cutoff=.5, thresh=1e-10,
)
print(result.coefficients.to_frame().join(result.labels))
print("alpha:", result.alpha, "lambda:", result.lambda_min)
print("cutoff reached:", result.cutoff_reached)
