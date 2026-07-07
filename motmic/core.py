"""Core MOT-MIC estimator."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from .ot import pairwise_cost, topk_filter, unbalanced_sinkhorn


@dataclass
class MOTMICResult:
    """Container for MOT-MIC outputs."""

    primary_ids: list[str]
    site_scores: pd.DataFrame
    transport_plans: dict[str, np.ndarray]
    filtered_plans: dict[str, np.ndarray]
    embedding: np.ndarray

    @property
    def pan_score(self) -> pd.Series:
        return self.site_scores.sum(axis=1)

    @property
    def organ_specificity(self) -> pd.Series:
        denom = self.pan_score.replace(0, np.nan)
        return self.site_scores.max(axis=1) / denom

    def to_frame(self) -> pd.DataFrame:
        out = self.site_scores.copy()
        out["pan_MIC_score"] = self.pan_score
        out["organ_specificity"] = self.organ_specificity.fillna(0)
        out["predicted_site"] = self.site_scores.idxmax(axis=1)
        return out


class MOTMIC:
    """Multi-organ optimal-transport MIC estimator."""

    def __init__(
        self,
        n_components: int = 20,
        epsilon: float = 0.05,
        rho: float = 1.0,
        top_k: int = 1,
        random_state: int = 13,
    ) -> None:
        self.n_components = n_components
        self.epsilon = epsilon
        self.rho = rho
        self.top_k = top_k
        self.random_state = random_state

    def fit_predict(
        self,
        primary_expr: pd.DataFrame,
        metastasis_expr_by_site: dict[str, pd.DataFrame],
    ) -> MOTMICResult:
        """Estimate organ-specific MIC scores for primary tumor cells."""

        if not metastasis_expr_by_site:
            raise ValueError("At least one metastatic site is required")

        common_genes = set(primary_expr.columns)
        for expr in metastasis_expr_by_site.values():
            common_genes &= set(expr.columns)
        common_genes = sorted(common_genes)
        if len(common_genes) < 2:
            raise ValueError("Primary and metastatic matrices must share genes")

        matrices = [primary_expr[common_genes]]
        site_order = list(metastasis_expr_by_site)
        matrices.extend(metastasis_expr_by_site[s][common_genes] for s in site_order)
        stacked = pd.concat(matrices, axis=0)

        scaled = StandardScaler().fit_transform(stacked.values)
        n_components = min(self.n_components, scaled.shape[0] - 1, scaled.shape[1])
        embedding = PCA(n_components=n_components, random_state=self.random_state).fit_transform(scaled)

        n_primary = len(primary_expr)
        primary_z = embedding[:n_primary]
        cursor = n_primary

        scores = {}
        plans = {}
        filtered_plans = {}
        for site in site_order:
            n_site = len(metastasis_expr_by_site[site])
            site_z = embedding[cursor : cursor + n_site]
            cursor += n_site
            cost = pairwise_cost(primary_z, site_z)
            plan = unbalanced_sinkhorn(cost, epsilon=self.epsilon, rho=self.rho)
            filtered = topk_filter(plan, k=self.top_k)
            raw_score = filtered.sum(axis=1)
            norm_score = raw_score / max(raw_score.max(), 1e-12)
            scores[site] = norm_score
            plans[site] = plan
            filtered_plans[site] = filtered

        score_frame = pd.DataFrame(scores, index=primary_expr.index)
        return MOTMICResult(
            primary_ids=list(primary_expr.index),
            site_scores=score_frame,
            transport_plans=plans,
            filtered_plans=filtered_plans,
            embedding=embedding,
        )

