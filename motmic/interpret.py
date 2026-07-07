"""Gene prioritization with SHAP or permutation fallback."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.inspection import permutation_importance


def rank_genes_with_shap(
    expr: pd.DataFrame,
    target: pd.Series,
    task: str = "classification",
    random_state: int = 11,
) -> pd.DataFrame:
    """Rank genes by SHAP values when available, otherwise permutation importance."""

    aligned = expr.loc[target.index]
    y = target.values
    if task == "classification":
        model = RandomForestClassifier(n_estimators=300, random_state=random_state, class_weight="balanced")
    else:
        model = RandomForestRegressor(n_estimators=300, random_state=random_state)
    model.fit(aligned, y)

    method = "permutation"
    try:
        import shap

        explainer = shap.TreeExplainer(model)
        shap_values = explainer.shap_values(aligned)
        if isinstance(shap_values, list):
            values = shap_values[-1]
        else:
            values = shap_values
        if values.ndim == 3:
            if values.shape[2] == 2:
                values = values[:, :, 1]
            else:
                values = values.mean(axis=2)
        importance = np.abs(values).mean(axis=0)
        method = "shap"
    except Exception:
        perm = permutation_importance(model, aligned, y, n_repeats=10, random_state=random_state)
        importance = perm.importances_mean

    ranked = pd.DataFrame({"gene": aligned.columns, "importance": importance, "method": method})
    ranked = ranked.sort_values("importance", ascending=False).reset_index(drop=True)
    ranked["rank"] = np.arange(1, len(ranked) + 1)
    return ranked
