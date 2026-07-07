"""MOT-MIC: multi-organ optimal transport for metastasis-initiating cells."""

from .core import MOTMIC, MOTMICResult
from .io import downsample_cells, find_10x_triplet, read_10x_mtx
from .simulation import simulate_metastasis_dataset
from .validation import evaluate_binary_labels, summarize_site_scores

__all__ = [
    "MOTMIC",
    "MOTMICResult",
    "downsample_cells",
    "find_10x_triplet",
    "read_10x_mtx",
    "simulate_metastasis_dataset",
    "evaluate_binary_labels",
    "summarize_site_scores",
]
