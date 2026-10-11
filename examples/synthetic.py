"""Small reproducible three-family example; no downloads, R, or AnnData."""
import numpy as np
import pandas as pd
from scipy import sparse
from scissor_py import scissor

rng = np.random.default_rng(123)
bulk = pd.DataFrame(rng.gamma(2, 2, (80, 30)), index=[f"gene_{i}" for i in range(80)], columns=[f"sample_{i}" for i in range(30)])
cells = pd.DataFrame(rng.gamma(2, 2, (80, 8)), index=bulk.index, columns=[f"cell_{i}" for i in range(8)])
graph = sparse.diags([np.ones(7), np.ones(7)], [-1, 1], format="csr")
foldid = np.arange(30) % 3
phenotypes = {
    "gaussian": pd.Series(rng.normal(size=30), index=bulk.columns),
    "binomial": pd.Series(np.arange(30) % 2, index=bulk.columns),
    "cox": pd.DataFrame({"time": np.repeat(np.arange(1, 16), 2), "status": (np.arange(30) % 5 != 0).astype(int)}, index=bulk.columns),
}

for family, phenotype in phenotypes.items():
    result = scissor(bulk, cells, phenotype, graph, graph_cell_ids=cells.columns,
                     family=family, semantics="corrected", alpha=[.2, .5],
                     lambdas=[.1, .04, .015], foldid=foldid, cutoff=.5)
    print(family, "alpha=", result.alpha, "lambda=", result.lambda_min,
          "positive=", result.positive_cells, "negative=", result.negative_cells,
          "cutoff_reached=", result.cutoff_reached)
