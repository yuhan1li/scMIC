"""MOT-MIC: multi-organ optimal transport for metastasis-initiating cells."""

from .core import MOTMIC, MOTMICResult
from .branch_embedding import BranchPreservingEmbedder, assign_organ_branch, compute_program_scores
from .io import downsample_cells, find_10x_triplet, read_10x_mtx
from .organ_prior import MIC_STATE_PRIORS, ORGAN_TROPISM_PRIORS
from .validation import evaluate_binary_labels, summarize_site_scores

__all__ = [
    "BranchPreservingEmbedder",
    "MIC_STATE_PRIORS",
    "MOTMIC",
    "MOTMICResult",
    "ORGAN_TROPISM_PRIORS",
    "assign_organ_branch",
    "compute_program_scores",
    "downsample_cells",
    "find_10x_triplet",
    "read_10x_mtx",
    "evaluate_binary_labels",
    "summarize_site_scores",
]
