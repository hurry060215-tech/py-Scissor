"""Full-cell/full-cohort acceptance of the existing public AnnData entry point.

Uses frozen transformed, common-gene inputs whose raw-to-joint construction was
separately checked against R. This is not a fresh raw-file preprocessing run.
No sampling, downloads, graph reconstruction, or AnnData writeback occurs.
"""
import argparse
import gc
import hashlib
import json
import resource
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from scissor_py import scissor_from_anndata
from run_luad_python import plain


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prepared", type=Path)
    parser.add_argument("saved_paths", type=Path, help="Complete Python tightened 1e-10 path/CV/refit records")
    parser.add_argument("output", type=Path)
    parser.add_argument("--family", choices=["all", "binomial", "cox", "gaussian"], default="all")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    manifest_bytes = (args.prepared/"manifest.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    graph = sparse.load_npz(args.prepared/"graph.npz")
    for family in ("binomial", "cox", "gaussian"):
        if args.family != "all" and family != args.family:
            continue
        before = time.monotonic()
        config = manifest["families"][family]
        cohort = manifest["cohorts"][config["cohort"]]
        n, p = cohort["X_shape"]
        joint = np.memmap(args.prepared/config["cohort"]/"joint.bin", dtype="<f8", mode="r",
                          shape=tuple(cohort["joint_shape"]), order="F")
        cells, genes, samples = cohort["cell_ids"], cohort["common_gene_ids"], cohort["sample_ids"]
        bulk = pd.DataFrame(joint[:, :n], index=genes, columns=samples)
        # All genes are reversed in AnnData to exercise ID alignment. The
        # selected layer is explicit; zero sparse X is intentionally not used.
        data = ad.AnnData(X=sparse.csr_matrix((p, len(genes))),
                          obs=pd.DataFrame({"preserved": pd.Categorical(["source"]*p)}, index=cells),
                          var=pd.DataFrame(index=genes[::-1]))
        data.layers["log1p_cp10k"] = joint[::-1, n:].T
        data.obsp["frozen_network"] = graph.copy()
        data.uns["preserved"] = {"cohort": config["cohort"]}
        if family == "cox":
            phenotype = pd.DataFrame(config["y"], index=samples, columns=["time", "status"])
        elif family == "binomial":
            phenotype = pd.Series(np.where(np.asarray(config["y"]) == 1, "mutant", "wildtype"), index=samples)
        else:
            phenotype = pd.Series(config["y"], index=samples)
        folds = pd.Series(config["foldid"], index=samples).iloc[::-1]
        original_obs, original_var = data.obs.copy(deep=True), data.var.copy(deep=True)
        print(json.dumps({"status": "starting", "family": family, "shape": [n, p], "common_genes": len(genes)}), flush=True)
        result = scissor_from_anndata(data, bulk, phenotype.iloc[::-1], layer="log1p_cp10k", graph_key="frozen_network",
                                     family=family, positive_label="mutant" if family == "binomial" else None,
                                     alpha=[x["alpha"] for x in config["paths"]], cutoff=config["cutoff"],
                                     foldid=folds, nlambda=20, semantics="corrected", thresh=1e-10,
                                     maxit=config["maxit"], max_dense_bytes=3_000_000_000,
                                     engine="numba" if family == "gaussian" else "numpy")
        pd.testing.assert_frame_equal(data.obs, original_obs)
        pd.testing.assert_frame_equal(data.var, original_var)
        assert data.uns == {"preserved": {"cohort": config["cohort"]}}
        assert data.X.nnz == 0 and (data.obsp["frozen_network"] != graph).nnz == 0
        assert np.shares_memory(data.layers["log1p_cp10k"], joint)
        assert result.prepared.gene_ids == tuple(genes)
        assert result.prepared.sample_ids == tuple(samples) and result.prepared.cell_ids == tuple(cells)
        expected_X = np.load(args.prepared/config["cohort"]/"X.npy")
        np.testing.assert_allclose(result.prepared.similarity, expected_X, atol=1e-8, rtol=1e-6)
        comparison = args.saved_paths/f"{family}_corrected_alpha_{result.alpha:g}"
        expected = np.load(comparison/"refit.npz")
        np.testing.assert_allclose(result.coefficients, expected["coefficients"][:, 0], atol=1e-5, rtol=1e-4)
        assert np.array_equal(np.sign(result.coefficients), np.sign(expected["coefficients"][:, 0]))
        expected_alpha = .6 if family == "binomial" else .05
        assert result.alpha == expected_alpha and result.cutoff_reached
        assert len(result.cv_results) == (3 if family == "binomial" else 1)
        cv_comparisons = []
        for cv in result.cv_results:
            saved = args.saved_paths/f"{family}_corrected_alpha_{cv.path.alpha:g}"
            old_cv, old_path = np.load(saved/"cv.npz"), np.load(saved/"full.npz")
            np.testing.assert_array_equal(cv.foldid, config["foldid"])
            np.testing.assert_allclose(cv.path.lambdas, old_path["lambdas"], atol=1e-12, rtol=1e-10)
            np.testing.assert_allclose(cv.mean_loss, old_cv["mean_loss"], atol=1e-6, rtol=1e-4)
            cv_comparisons.append({"alpha": cv.path.alpha, "actual_lambdas": cv.path.lambdas,
                                   "lambda_max_difference": np.max(np.abs(cv.path.lambdas-old_path["lambdas"])),
                                   "cv_mean_max_difference": np.max(np.abs(cv.mean_loss-old_cv["mean_loss"]))})
        record = {"status": "passed", "family": family, "shape": [n, p], "common_genes": len(genes),
                  "scope": "full transformed common-gene input through public AnnData adapter, automatic path, explicit folds, ordered alpha search and cold refit",
                  "raw_to_joint_validation": "separate actual full R/preprocessCore run",
                  "manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
                  "sampling": "none", "threshold": 1e-10, "alpha": result.alpha, "lambda_min": result.lambda_min,
                  "positive_count": len(result.positive_cells), "negative_count": len(result.negative_cells),
                  "cutoff_reached": result.cutoff_reached, "input_not_mutated": True,
                  "reversed_gene_phenotype_and_fold_ID_alignment_passed": True,
                  "similarity_max_difference": np.max(np.abs(result.prepared.similarity-expected_X)),
                  "refit_coefficient_max_difference": np.max(np.abs(result.coefficients.to_numpy()-expected["coefficients"][:, 0])),
                  "metadata": result.metadata, "cv_comparisons": cv_comparisons,
                  "seconds": time.monotonic()-before,
                  "process_peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024}
        (args.output/f"{family}.json").write_text(json.dumps(plain(record), indent=2, allow_nan=False)+"\n")
        result.coefficients.to_csv(args.output/f"{family}-coefficients.csv")
        print(json.dumps(plain({k: v for k, v in record.items() if k not in ("metadata", "cv_comparisons")})), flush=True)
        del result, data, joint, bulk, expected_X, expected
        gc.collect()


if __name__ == "__main__":
    main()
