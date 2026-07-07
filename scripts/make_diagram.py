from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle


def box(ax, xy, text, color):
    patch = FancyBboxPatch(
        xy,
        2.4,
        0.72,
        boxstyle="round,pad=0.03,rounding_size=0.04",
        linewidth=1.2,
        facecolor=color,
        edgecolor="#333333",
    )
    ax.add_patch(patch)
    ax.text(xy[0] + 1.2, xy[1] + 0.36, text, ha="center", va="center", fontsize=9)


def arrow(ax, start, end):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=12, lw=1.2, color="#333333"))


def branch_arrow(ax, start, end, color, label):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            connectionstyle="arc3,rad=0.12",
            arrowstyle="-|>",
            mutation_scale=14,
            lw=2.0,
            color=color,
        )
    )
    ax.text(end[0] + 0.1, end[1], label, ha="left", va="center", fontsize=9, color=color, weight="bold")


def draw_branch_panel(ax):
    root = (1.0, 1.25)
    trunk = [(1.0, 1.25), (2.0, 1.6), (3.1, 1.9), (4.1, 2.0)]
    ax.plot([p[0] for p in trunk], [p[1] for p in trunk], color="#2f3437", lw=3)
    for i, (x, y) in enumerate(trunk):
        ax.add_patch(Circle((x, y), 0.09, facecolor=plt.cm.viridis(i / max(len(trunk) - 1, 1)), edgecolor="white", lw=0.8))
    ax.text(root[0] - 0.15, root[1] - 0.3, "primary\nlow-MIC", ha="right", va="top", fontsize=8)
    ax.text(3.0, 2.0, "shared MIC trunk", ha="center", va="bottom", fontsize=9, weight="bold")

    branch_arrow(ax, trunk[-1], (5.7, 3.45), "#7a9b42", "liver branch")
    branch_arrow(ax, trunk[-1], (5.8, 2.6), "#7b5ea7", "lung branch")
    branch_arrow(ax, trunk[-1], (5.6, 1.8), "#b5893b", "bone branch")
    branch_arrow(ax, trunk[-1], (5.55, 1.1), "#4f9a9a", "lymph-node branch")
    branch_arrow(ax, trunk[-1], (5.4, 0.5), "#bf5f82", "brain branch")

    rng = np.random.default_rng(2)
    colors = ["#92adc9", "#d98674", "#aabf80", "#a896c8"]
    centers = [(1.1, 1.3), (2.4, 1.65), (4.2, 1.95), (5.4, 2.35)]
    for c, color in zip(centers, colors):
        pts = rng.normal(c, (0.18, 0.12), size=(55, 2))
        ax.scatter(pts[:, 0], pts[:, 1], s=9, color=color, alpha=0.35, linewidths=0)


def main():
    root = Path(__file__).resolve().parents[1]
    out = root / "figures" / "motmic_algorithm_schematic.png"
    branch_out = root / "figures" / "branch_preserving_scmic_schematic.png"
    out.parent.mkdir(exist_ok=True)

    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.set_xlim(0, 11)
    ax.set_ylim(0, 5.5)
    ax.axis("off")

    box(ax, (0.4, 4.1), "Paired primary\nscRNA-seq", "#DDEEFF")
    box(ax, (0.4, 2.9), "Metastatic sites\nlung/liver/bone", "#E6F3D8")
    box(ax, (3.0, 3.5), "Tumor-cell filter\nmarkers + CNV", "#FFF0CC")
    box(ax, (5.6, 3.5), "Latent space\nPCA/scVI/AE", "#F5DDFF")
    box(ax, (8.2, 4.1), "UOT per organ\n+ top-k origins", "#FFE1D6")
    box(ax, (8.2, 2.9), "MIC score vector\npan + organ bias", "#DFF6F0")
    box(ax, (5.6, 1.1), "Validation\nlineage/clinical/spatial", "#EDEDED")
    box(ax, (8.2, 1.1), "SHAP ranking\nkey genes", "#FBE4A7")

    arrow(ax, (2.8, 4.45), (3.0, 4.0))
    arrow(ax, (2.8, 3.25), (3.0, 3.85))
    arrow(ax, (5.4, 3.86), (5.6, 3.86))
    arrow(ax, (8.0, 3.86), (8.2, 4.3))
    arrow(ax, (9.4, 4.1), (9.4, 3.62))
    arrow(ax, (8.2, 3.2), (7.9, 1.84))
    arrow(ax, (8.2, 3.2), (9.0, 1.84))

    ax.text(0.4, 0.35, "Gold standard priority: GSE173958 lineage tracing > paired human cohorts > spatial/bulk validation", fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(12.5, 6.4))
    ax.set_xlim(0, 12.5)
    ax.set_ylim(0, 6.4)
    ax.axis("off")

    box(ax, (0.4, 5.25), "Paired primary +\nmetastatic scRNA-seq", "#DDEEFF")
    box(ax, (3.1, 5.25), "Nuisance removal\nQC/cell cycle/batch", "#EFEFEF")
    box(ax, (5.8, 5.25), "SAKURA-inspired\nknowledge priors", "#FFF0CC")
    box(ax, (8.5, 5.25), "Branch-preserving\nscTour latent", "#F5DDFF")

    arrow(ax, (2.8, 5.61), (3.1, 5.61))
    arrow(ax, (5.5, 5.61), (5.8, 5.61))
    arrow(ax, (8.2, 5.61), (8.5, 5.61))

    ax.text(6.95, 4.85, "Do not regress out organ programs", ha="center", fontsize=9, color="#A35A00")
    ax.text(6.95, 4.55, "Boost MIC + liver/lung/bone/LN/brain prior genes", ha="center", fontsize=9, color="#A35A00")

    draw_branch_panel(ax)
    box(ax, (7.2, 3.4), "UOT from primary\ninto each organ branch", "#FFE1D6")
    box(ax, (9.9, 3.4), "Branch assignment\nP(liver/lung/bone/LN/brain)", "#DFF6F0")
    box(ax, (7.2, 2.05), "Validation\nlineage + organ labels", "#E6F3D8")
    box(ax, (9.9, 2.05), "SHAP / DE genes\nper branch", "#FBE4A7")

    arrow(ax, (6.0, 2.05), (7.2, 3.75))
    arrow(ax, (9.6, 3.75), (9.9, 3.75))
    arrow(ax, (8.4, 3.4), (8.4, 2.77))
    arrow(ax, (11.1, 3.4), (11.1, 2.77))

    ax.text(
        0.45,
        0.05,
        "Output: shared MIC trunk + organotropic branches. Remove technical variation, preserve true organ-specific metastatic biology.",
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(branch_out, dpi=240)
    plt.close(fig)
    print(out)
    print(branch_out)


if __name__ == "__main__":
    main()
