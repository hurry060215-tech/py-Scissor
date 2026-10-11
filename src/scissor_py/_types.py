"""Public, auditable records; no R object or AnnData runtime dependency."""
from dataclasses import dataclass, field
from typing import Any
import numpy as np
import pandas as pd
from scipy import sparse


REFERENCE_COMMIT = "311560a1f160f665917e7da0fcf289914ee99773"


class ScissorError(ValueError):
    """Invalid scientific input or unsupported compatibility condition."""


class ConvergenceError(RuntimeError):
    """No converged fit satisfies the requested operation."""


@dataclass
class PreparedData:
    similarity: np.ndarray
    network: sparse.csr_matrix
    gene_ids: tuple
    sample_ids: tuple
    cell_ids: tuple
    normalized_bulk: np.ndarray
    normalized_single_cell: np.ndarray
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class PathResult:
    family: str
    alpha: float
    lambdas: np.ndarray
    coefficients: np.ndarray  # cells × lambdas
    intercepts: np.ndarray
    converged: np.ndarray
    iterations: np.ndarray
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class CVResult:
    path: PathResult
    foldid: np.ndarray
    mean_loss: np.ndarray
    standard_error: np.ndarray
    fold_loss: np.ndarray
    lambda_min: float
    lambda_1se: float
    valid: np.ndarray
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class ScissorResult:
    coefficients: pd.Series
    labels: pd.Series
    alpha: float
    lambda_min: float
    intercept: float
    family: str
    cutoff_reached: bool
    cv_results: list[CVResult]
    prepared: PreparedData
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def positive_cells(self) -> list:
        return self.labels.index[self.labels == "Scissor+"].tolist()

    @property
    def negative_cells(self) -> list:
        return self.labels.index[self.labels == "Scissor-"].tolist()

