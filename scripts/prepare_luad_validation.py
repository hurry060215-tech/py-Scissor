"""Prepare every official LUAD cell/sample for an explicit-graph kernel audit.

Development command only; requires rdata. Does not reproduce Seurat clustering
or the published tutorial's selected-cell counts. Input reconciliation is
explicit here, never hidden inside the strict library API.
"""
import argparse
import gc
import hashlib
import json
import platform
import resource
import time
from pathlib import Path

import numpy as np
import pandas as pd
import rdata
import scipy
from scipy import io, sparse
from scipy.sparse.linalg import LinearOperator, svds
from scipy.spatial import cKDTree

from scissor_py import prepare
from scissor_py._api import _lambda_path


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def binary(path, *arrays):
    with path.open("wb") as stream:
        for array in arrays:
            # Column writes bound temporary storage even for strided views.
            for j in range(array.shape[1]):
                stream.write(np.asarray(array[:, j], dtype="<f8").tobytes())


def stratified_folds(labels, seed=123):
    rng = np.random.default_rng(seed)
    folds = np.empty(len(labels), dtype=np.int64)
    for label in np.unique(labels):
        ids = np.flatnonzero(labels == label)
        ids = rng.permutation(ids)
        folds[ids] = np.arange(len(ids)) % 5 + 1
    return folds


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    record = args.output / "manifest.json"
    if record.exists():
        raise FileExistsError(f"Refusing to replace frozen manifest {record}")
    start = time.monotonic()
    metadata = {
        "status": "preparing", "sampling": "none",
        "source_download_manifest": json.loads((args.source / "download-manifest.json").read_text()),
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "scipy": scipy.__version__, "pandas": pd.__version__, "rdata": rdata.__version__},
        "binary_format": "little-endian float64, column-major, no header",
        "single_cell_transform": "log1p(counts / per-cell total of all 33694 source genes * 10000)",
        "bulk_transform": "none; original preprocessed tutorial matrices",
        "duplicate_gene_policy": "validation-only explicit first source row per gene, matching R character matrix indexing after intersect; library itself rejects duplicates",
        "cohorts": {}, "families": {},
    }
    sc = rdata.read_rda(args.source / "scRNA-seq.RData")["sc_dataset"]
    genes = sc.coords[sc.dims[0]].values.astype(str)
    cells = sc.coords[sc.dims[1]].values.astype(str)
    assert pd.Index(genes).is_unique and pd.Index(cells).is_unique
    matrix = np.asarray(sc.values)
    assert np.isfinite(matrix).all() and np.min(matrix) >= 0
    counts = sparse.csc_matrix(matrix)
    del sc, matrix
    gc.collect()
    totals = np.asarray(counts.sum(axis=0)).ravel()
    assert np.all(totals > 0)
    counts.data *= np.repeat(10000 / totals, np.diff(counts.indptr))
    np.log1p(counts.data, out=counts.data)
    metadata["single_cell"] = {"shape": list(counts.shape), "cell_ids": cells.tolist(),
                               "gene_ids": genes.tolist(), "nnz": counts.nnz,
                               "raw_library_size_min": float(totals.min()),
                               "raw_library_size_max": float(totals.max())}
    np.save(args.output / "single_cell_log1p.npy", counts.toarray(order="F"))
    # Explicit all-gene, all-cell PCA; no HVG selection or downsampling.
    a = counts.T.tocsr()
    mean = np.asarray(a.mean(axis=0)).ravel()
    centered = LinearOperator(
        a.shape, dtype=np.float64,
        matvec=lambda v: np.asarray(a @ v).ravel() - np.dot(mean, v),
        rmatvec=lambda u: np.asarray(a.T @ u).ravel() - mean * np.sum(u),
        matmat=lambda v: a @ v - (mean @ v)[None, :],
        rmatmat=lambda u: a.T @ u - mean[:, None] * np.sum(u, axis=0)[None, :],
    )
    print("PCA: all cells and genes", flush=True)
    u, s, _ = svds(centered, k=20, tol=1e-8, solver="arpack", rng=np.random.default_rng(123))
    order = np.argsort(s)[::-1]
    scores = u[:, order] * s[order]
    np.save(args.output / "graph_pca_scores.npy", scores)
    distances, neighbors = cKDTree(scores).query(scores, k=21, workers=1)
    assert np.all(neighbors[:, 0] == np.arange(len(cells)))
    # Reject ambiguous neighbor cutoffs rather than hiding an arbitrary tie.
    extra_distances, _ = cKDTree(scores).query(scores, k=22, workers=1)
    assert np.all(extra_distances[:, 20] != extra_distances[:, 21])
    graph = sparse.csr_matrix((np.ones(len(cells)*20),
                             (np.repeat(np.arange(len(cells)), 20), neighbors[:, 1:].ravel())),
                            shape=(len(cells), len(cells)))
    graph = graph.maximum(graph.T)
    graph.sort_indices()
    sparse.save_npz(args.output / "graph.npz", graph)
    io.mmwrite(args.output / "graph.mtx", graph, symmetry="general")
    metadata["graph"] = {"method": "union-symmetrized Euclidean 20-NN of 20 centered, unscaled all-gene log1p(CP10K) principal components",
                         "not_seurat_snn": True, "random_seed": 123, "shape": list(graph.shape),
                         "nnz": graph.nnz, "isolated_cells": 0, "pca_singular_values": s[order].tolist()}
    del counts, a, centered, u
    gc.collect()
    sc_frame = pd.DataFrame(np.load(args.output / "single_cell_log1p.npy", mmap_mode="r"),
                            index=genes, columns=cells, copy=False)
    survival = rdata.read_rda(args.source / "TCGA_LUAD_survival.RData")["bulk_survival"]
    survival = survival.set_index("TCGA_patient_barcode").rename(columns={"OS_time": "time", "Status": "status"})
    mutation = rdata.read_rda(args.source / "TCGA_LUAD_TP53_mutation.RData")["TP53_mutation"].to_series()
    all_duplicates = []
    for cohort, pheno in (("exp1", survival), ("exp2", mutation)):
        raw = rdata.read_rda(args.source / f"TCGA_LUAD_{cohort}.RData")["bulk_dataset"]
        row_ids = pd.Index(raw.coords[raw.dims[0]].values.astype(str))
        samples = pd.Index(raw.coords[raw.dims[1]].values.astype(str))
        assert samples.is_unique and not samples.hasnans
        assert pheno.index.is_unique and set(pheno.index) == set(samples)
        y = pheno.reindex(samples).to_numpy(dtype=float)
        assert np.isfinite(y).all()
        keep = ~row_ids.duplicated(keep="first")
        for i in np.flatnonzero(~keep):
            first = int(np.flatnonzero(row_ids == row_ids[i])[0])
            all_duplicates.append({"cohort": cohort, "gene_id": row_ids[i],
                                   "source_row_1based": int(i+1), "used_source_row_1based": first+1,
                                   "shared_with_single_cell": bool(row_ids[i] in sc_frame.index),
                                   "identical_expression_to_first": bool(np.array_equal(raw.values[i], raw.values[first]))})
        bulk = pd.DataFrame(raw.values[keep], index=row_ids[keep], columns=samples)
        assert np.isfinite(bulk.to_numpy()).all()
        del raw
        gc.collect()
        common = bulk.index[bulk.index.isin(sc_frame.index)]
        folder = args.output / cohort
        folder.mkdir()
        binary(folder / "joint.bin", bulk.loc[common].to_numpy(), sc_frame.loc[common].to_numpy())
        print(f"PREPARE {cohort}: {len(common)} genes x ({len(samples)} samples + {len(cells)} cells)", flush=True)
        prepared = prepare(bulk, sc_frame, graph, graph_cell_ids=cells, max_dense_bytes=4_000_000_000)
        del bulk
        gc.collect()
        binary(folder / "normalized.bin", prepared.normalized_bulk, prepared.normalized_single_cell)
        binary(folder / "X.bin", prepared.similarity)
        np.save(folder / "X.npy", prepared.similarity)
        info = {"source_bulk_shape": [len(row_ids), len(samples)], "unique_bulk_genes": int(keep.sum()),
                "duplicate_rows": int((~keep).sum()), "common_gene_ids": list(prepared.gene_ids),
                "sample_ids": list(prepared.sample_ids), "cell_ids": list(prepared.cell_ids),
                "joint_shape": [len(common), len(samples)+len(cells)], "X_shape": list(prepared.similarity.shape),
                "phenotype_already_source_order": bool(np.array_equal(pheno.index.to_numpy(dtype=str), samples.to_numpy(dtype=str))),
                "preprocessing_metadata": prepared.metadata, "similarity_quantiles": np.quantile(prepared.similarity, [0,.25,.5,.75,1]).tolist()}
        metadata["cohorts"][cohort] = info
        if cohort == "exp1":
            assert np.all(y[:, 0] >= 0) and set(np.unique(y[:, 1])) == {0, 1}
            folds = stratified_folds(y[:, 1])
            outcomes = {"cox": y,
                        "gaussian": 2+2*prepared.similarity[:, 0]-1.5*prepared.similarity[:, 137]+prepared.similarity[:, 777]+.02*np.sin(np.arange(1, len(samples)+1))}
            info["survival"] = {"events": int(y[:, 1].sum()), "censored": int((y[:, 1] == 0).sum()),
                                "zero_times": int((y[:, 0] == 0).sum()), "unique_times": len(np.unique(y[:, 0])),
                                "range_time": [float(y[:, 0].min()), float(y[:, 0].max())]}
        else:
            assert set(np.unique(y)) == {0, 1}
            folds = stratified_folds(y)
            outcomes = {"binomial": y}
            info["positive_class"] = "TP53 mutation = 1; wild-type = 0"
            info["class_counts"] = {str(i): int((y == i).sum()) for i in (0, 1)}
        for family, outcome in outcomes.items():
            paths = []
            for alpha in (.05, .2, .6):
                lam, _, _ = _lambda_path(prepared.similarity, outcome, family, alpha, None, 20, .01)
                paths.append({"alpha": alpha, "lambdas": lam.tolist()})
            metadata["families"][family] = {"cohort": cohort, "y": outcome.tolist(), "foldid": folds.tolist(),
                                             "paths": paths, "thresh": 1e-8, "maxit": 100000, "cutoff": .2,
                                             "gaussian_signal": "2+2*X[,1]-1.5*X[,138]+X[,778]+0.02*sin(1:n); synthetic continuous response" if family == "gaussian" else None}
        del prepared
        gc.collect()
    pd.DataFrame(all_duplicates).to_csv(args.output / "duplicate-gene-map.tsv", sep="\t", index=False)
    metadata["files"] = {str(p.relative_to(args.output)): {"bytes": p.stat().st_size, "sha256": digest(p)}
                         for p in sorted(args.output.rglob("*")) if p.is_file()}
    metadata.update(status="frozen", elapsed_seconds=time.monotonic()-start,
                    process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    record.write_text(json.dumps(metadata, indent=2)+"\n")
    print(json.dumps({"manifest": str(record), "elapsed_seconds": metadata["elapsed_seconds"],
                      "process_peak_rss_bytes": metadata["process_peak_rss_bytes"]}), flush=True)


if __name__ == "__main__":
    main()
