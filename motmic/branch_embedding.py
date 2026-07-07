"""SAKURA-inspired branch-preserving embedding utilities.

The goal is to preserve subtle organotropic programs while removing technical
covariates. This module deliberately separates two operations that are often
mixed together:

1. nuisance residualization removes non-biological effects such as library size,
   mitochondrial percentage, or cell-cycle scores;
2. knowledge-guided weighting boosts MIC and organ-prior genes so rare branch
   signals remain visible to downstream scTour / OT / graph learning.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

from .organ_prior import MIC_STATE_PRIORS, ORGAN_TROPISM_PRIORS, match_prior_genes


def _as_frame(expr: pd.DataFrame) -> pd.DataFrame:
    if not isinstance(expr, pd.DataFrame):
        raise TypeError("expr must be a pandas DataFrame with cells as rows and genes as columns")
    if expr.empty:
        raise ValueError("expr is empty")
    return expr


def compute_program_scores(expr: pd.DataFrame, priors: Mapping[str, list[str]]) -> pd.DataFrame:
    """Score prior programs by average log1p expression of matched genes."""

    expr = _as_frame(expr)
    matched = match_prior_genes(list(expr.columns), priors)
    logged = np.log1p(expr)
    scores = {}
    for name, genes in matched.items():
        if genes:
            scores[name] = logged[genes].mean(axis=1)
        else:
            scores[name] = pd.Series(0.0, index=expr.index)
    return pd.DataFrame(scores, index=expr.index)


def residualize_nuisance(expr: pd.DataFrame, covariates: pd.DataFrame | None) -> pd.DataFrame:
    """Regress out nuisance covariates from log1p expression.

    Organ labels, metastatic-site labels, and organ-prior scores should not be
    passed as covariates; those are biological branch signals to preserve.
    """

    expr = _as_frame(expr)
    logged = pd.DataFrame(np.log1p(expr), index=expr.index, columns=expr.columns)
    if covariates is None or covariates.empty:
        return logged

    cov = covariates.loc[expr.index].copy()
    cov = cov.select_dtypes(include=[np.number]).replace([np.inf, -np.inf], np.nan)
    cov = cov.fillna(cov.median(numeric_only=True))
    if cov.shape[1] == 0:
        return logged

    x = StandardScaler().fit_transform(cov.values)
    y = logged.values
    model = LinearRegression().fit(x, y)
    residual = y - model.predict(x)
    residual += y.mean(axis=0, keepdims=True)
    residual = np.maximum(residual, 0.0)
    return pd.DataFrame(residual, index=expr.index, columns=expr.columns)


@dataclass
class BranchEmbeddingResult:
    """Outputs from branch-preserving embedding."""

    embedding: pd.DataFrame
    organ_scores: pd.DataFrame
    mic_program_scores: pd.DataFrame
    gene_weights: pd.Series
    matched_organ_priors: dict[str, list[str]]
    matched_mic_priors: dict[str, list[str]]


class BranchPreservingEmbedder:
    """Knowledge-guided PCA embedding that preserves organotropic branches."""

    def __init__(
        self,
        organ_priors: Mapping[str, list[str]] | None = None,
        mic_priors: Mapping[str, list[str]] | None = None,
        prior_weight: float = 3.0,
        n_components: int = 20,
        random_state: int = 13,
    ) -> None:
        self.organ_priors = dict(organ_priors or ORGAN_TROPISM_PRIORS)
        self.mic_priors = dict(mic_priors or MIC_STATE_PRIORS)
        self.prior_weight = prior_weight
        self.n_components = n_components
        self.random_state = random_state

    def fit_transform(self, expr: pd.DataFrame, nuisance_covariates: pd.DataFrame | None = None) -> BranchEmbeddingResult:
        """Create a branch-preserving embedding from expression data."""

        expr = _as_frame(expr)
        corrected = residualize_nuisance(expr, nuisance_covariates)
        matched_organ = match_prior_genes(list(corrected.columns), self.organ_priors)
        matched_mic = match_prior_genes(list(corrected.columns), self.mic_priors)

        prior_genes = set()
        for genes in [*matched_organ.values(), *matched_mic.values()]:
            prior_genes.update(genes)
        weights = pd.Series(1.0, index=corrected.columns)
        if prior_genes:
            weights.loc[list(prior_genes)] = self.prior_weight

        weighted = corrected.mul(weights, axis=1)
        scaled = StandardScaler().fit_transform(weighted.values)
        n_components = min(self.n_components, scaled.shape[0] - 1, scaled.shape[1])
        z = PCA(n_components=n_components, random_state=self.random_state).fit_transform(scaled)
        embedding = pd.DataFrame(
            z,
            index=expr.index,
            columns=[f"branch_pc_{i + 1}" for i in range(n_components)],
        )
        return BranchEmbeddingResult(
            embedding=embedding,
            organ_scores=compute_program_scores(corrected, self.organ_priors),
            mic_program_scores=compute_program_scores(corrected, self.mic_priors),
            gene_weights=weights,
            matched_organ_priors=matched_organ,
            matched_mic_priors=matched_mic,
        )


def assign_organ_branch(site_scores: pd.DataFrame, min_specificity: float = 0.0) -> pd.Series:
    """Assign each primary cell to the organ branch with the strongest score."""

    if site_scores.empty:
        raise ValueError("site_scores is empty")
    denom = site_scores.sum(axis=1).replace(0, np.nan)
    specificity = site_scores.max(axis=1) / denom
    branch = site_scores.idxmax(axis=1).astype(str)
    branch = branch.where(specificity.fillna(0) >= min_specificity, "shared_MIC_trunk")
    return branch.rename("predicted_organ_branch")

