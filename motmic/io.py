"""Input helpers for public single-cell metastasis datasets."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from scipy import io, sparse


def read_10x_mtx(matrix_path: str | Path, features_path: str | Path, barcodes_path: str | Path) -> pd.DataFrame:
    """Read a 10x matrix into a cell-by-gene DataFrame.

    GEO supplementary files for GSE173958 and GSE249057 are provided as
    `matrix.mtx.gz`, `features.tsv.gz`, and `barcodes.tsv.gz`. The matrix is
    stored as genes by cells, so this helper transposes it to cells by genes.
    """

    matrix = io.mmread(str(matrix_path)).tocsr().T
    features = pd.read_csv(features_path, sep="\t", header=None)
    barcodes = pd.read_csv(barcodes_path, sep="\t", header=None)[0].astype(str)

    if features.shape[1] >= 2:
        genes = features.iloc[:, 1].astype(str)
    else:
        genes = features.iloc[:, 0].astype(str)

    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    return pd.DataFrame(matrix, index=barcodes.values, columns=genes.values)


def downsample_cells(expr: pd.DataFrame, n_cells: int = 2000, random_state: int = 13) -> pd.DataFrame:
    """Downsample cells for fast tutorial runs."""

    if len(expr) <= n_cells:
        return expr
    return expr.sample(n=n_cells, random_state=random_state)


def find_10x_triplet(raw_dir: str | Path, sample_token: str) -> tuple[Path, Path, Path]:
    """Find matrix, features, and barcode files for a sample token."""

    raw_dir = Path(raw_dir)
    matrices = sorted(raw_dir.glob(f"*{sample_token}*matrix.mtx.gz"))
    features = sorted(raw_dir.glob(f"*{sample_token}*features.tsv.gz"))
    barcodes = sorted(raw_dir.glob(f"*{sample_token}*barcodes.tsv.gz"))
    if not (matrices and features and barcodes):
        raise FileNotFoundError(f"Could not find complete 10x triplet for token {sample_token!r} in {raw_dir}")
    return matrices[0], features[0], barcodes[0]

