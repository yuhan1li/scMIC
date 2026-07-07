from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


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


def main():
    root = Path(__file__).resolve().parents[1]
    out = root / "figures" / "motmic_algorithm_schematic.png"
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
    print(out)


if __name__ == "__main__":
    main()

