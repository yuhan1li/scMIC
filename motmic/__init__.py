"""MOT-MIC: multi-organ optimal transport for metastasis-initiating cells."""

from .core import MOTMIC, MOTMICResult
from .simulation import simulate_metastasis_dataset
from .validation import evaluate_binary_labels, summarize_site_scores

__all__ = [
    "MOTMIC",
    "MOTMICResult",
    "simulate_metastasis_dataset",
    "evaluate_binary_labels",
    "summarize_site_scores",
]

