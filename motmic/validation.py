"""Validation helpers."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


def evaluate_binary_labels(scores: pd.Series, labels: pd.Series) -> dict[str, float]:
    """Evaluate MIC scores against binary labels."""

    y = labels.astype(int).values
    s = scores.loc[labels.index].values
    out = {
        "mean_positive_score": float(np.mean(s[y == 1])),
        "mean_negative_score": float(np.mean(s[y == 0])),
    }
    if len(np.unique(y)) == 2:
        out["AUROC"] = float(roc_auc_score(y, s))
        out["AUPRC"] = float(average_precision_score(y, s))
    return out


def summarize_site_scores(result_frame: pd.DataFrame, truth: pd.DataFrame) -> pd.DataFrame:
    """Summarize predicted site scores by true site label."""

    merged = result_frame.join(truth)
    score_cols = [c for c in result_frame.columns if c.endswith("_MIC_score") is False]
    score_cols = [c for c in score_cols if c in ["lung", "liver", "bone", "brain"]]
    return merged.groupby("true_site")[score_cols + ["pan_MIC_score", "organ_specificity"]].mean()

