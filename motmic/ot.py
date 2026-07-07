"""Small optimal-transport utilities used by MOT-MIC.

The implementation intentionally avoids a hard dependency on POT so the example
can run in lean server environments. For publication-scale analysis, POT or
moscot should be preferred for speed and numerical diagnostics.
"""

from __future__ import annotations

import numpy as np
from scipy.spatial.distance import cdist


def pairwise_cost(x: np.ndarray, y: np.ndarray, metric: str = "sqeuclidean") -> np.ndarray:
    """Return a non-negative cell-cell cost matrix."""

    cost = cdist(x, y, metric=metric)
    if metric == "euclidean":
        cost = cost**2
    scale = np.median(cost[cost > 0]) if np.any(cost > 0) else 1.0
    return cost / max(scale, 1e-12)


def unbalanced_sinkhorn(
    cost: np.ndarray,
    source_mass: np.ndarray | None = None,
    target_mass: np.ndarray | None = None,
    epsilon: float = 0.05,
    rho: float = 1.0,
    max_iter: int = 1000,
    tol: float = 1e-8,
) -> np.ndarray:
    """Compute an entropic unbalanced OT plan.

    This is the generalized Sinkhorn update for KL-relaxed marginals:

    u = (a / K v) ** tau
    v = (b / K.T u) ** tau
    tau = rho / (rho + epsilon)
    """

    n_source, n_target = cost.shape
    a = np.ones(n_source) / n_source if source_mass is None else source_mass.astype(float)
    b = np.ones(n_target) / n_target if target_mass is None else target_mass.astype(float)
    a = a / max(a.sum(), 1e-12)
    b = b / max(b.sum(), 1e-12)

    kernel = np.exp(-cost / max(epsilon, 1e-12))
    kernel = np.maximum(kernel, 1e-300)
    tau = rho / (rho + epsilon)

    u = np.ones(n_source)
    v = np.ones(n_target)
    for _ in range(max_iter):
        u_prev = u.copy()
        kv = kernel @ v
        u = (a / np.maximum(kv, 1e-300)) ** tau
        ktu = kernel.T @ u
        v = (b / np.maximum(ktu, 1e-300)) ** tau
        if np.linalg.norm(u - u_prev, ord=1) < tol:
            break

    plan = (u[:, None] * kernel) * v[None, :]
    return plan


def topk_filter(plan: np.ndarray, k: int = 1) -> np.ndarray:
    """Keep only the top-k primary-origin weights for each metastatic cell."""

    if k <= 0:
        raise ValueError("k must be positive")
    filtered = np.zeros_like(plan)
    k_eff = min(k, plan.shape[0])
    top_idx = np.argpartition(plan, -k_eff, axis=0)[-k_eff:, :]
    cols = np.arange(plan.shape[1])
    for row_positions in top_idx:
        filtered[row_positions, cols] = plan[row_positions, cols]
    return filtered

