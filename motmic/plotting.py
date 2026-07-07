"""Plotting utilities."""

from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd


def plot_site_scores(result_frame: pd.DataFrame, output: str) -> None:
    score_cols = [c for c in ["lung", "liver", "bone", "brain"] if c in result_frame.columns]
    fig, ax = plt.subplots(figsize=(7, 4))
    result_frame[score_cols].mean().sort_values(ascending=False).plot(kind="bar", ax=ax, color="#4C78A8")
    ax.set_ylabel("Mean organ-specific MIC score")
    ax.set_xlabel("Metastatic site")
    ax.set_title("MOT-MIC organotropic scores")
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)

