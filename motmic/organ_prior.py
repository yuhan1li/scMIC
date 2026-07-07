"""Knowledge priors for organotropic metastasis branches.

The gene lists are intentionally compact and editable. They are not meant to be
final biomarkers; they are anchors that help preserve rare organ-branch signals
when learning a branch-aware latent representation.
"""

from __future__ import annotations

from collections.abc import Mapping


ORGAN_TROPISM_PRIORS: dict[str, list[str]] = {
    "liver": [
        "CXCR4",
        "CXCL12",
        "MET",
        "SPP1",
        "ITGA5",
        "ITGB1",
        "ALDH1A1",
        "CEACAM5",
        "CEACAM6",
        "MUC1",
        "TGFBI",
        "LGALS3",
    ],
    "lung": [
        "SPARC",
        "VCAM1",
        "ITGB1",
        "ITGA3",
        "MMP2",
        "MMP9",
        "ANGPTL4",
        "EREG",
        "COX2",
        "PTGS2",
        "TNC",
        "L1CAM",
    ],
    "bone": [
        "CXCR4",
        "IL11",
        "PTHLH",
        "TGFB1",
        "MMP1",
        "MMP13",
        "SPP1",
        "IBSP",
        "RUNX2",
        "VCAM1",
        "ITGAV",
        "ITGB3",
    ],
    "lymph_node": [
        "CCR7",
        "CCL19",
        "CCL21",
        "CXCR4",
        "SELL",
        "ICAM1",
        "ITGAL",
        "ITGB2",
        "VEGFC",
        "VEGFD",
        "PDPN",
        "LYVE1",
    ],
    "brain": [
        "ST6GALNAC5",
        "L1CAM",
        "SERPINE1",
        "COX2",
        "PTGS2",
        "HBEGF",
        "ANGPTL4",
        "CTSL",
        "MMP2",
        "MMP9",
        "S100A4",
        "VIM",
    ],
}

MIC_STATE_PRIORS: dict[str, list[str]] = {
    "hybrid_emt": [
        "EPCAM",
        "KRT8",
        "KRT18",
        "KRT19",
        "VIM",
        "ZEB1",
        "ZEB2",
        "SNAI1",
        "SNAI2",
        "TWIST1",
    ],
    "survival_stress": [
        "AXL",
        "L1CAM",
        "CD44",
        "ITGA5",
        "ITGB1",
        "HIF1A",
        "VEGFA",
        "BCL2",
        "MCL1",
        "SOD2",
    ],
    "invasion_ecm": [
        "FN1",
        "SPARC",
        "TNC",
        "COL1A1",
        "COL3A1",
        "MMP2",
        "MMP9",
        "MMP14",
        "LOX",
        "TGFBI",
    ],
}

DEFAULT_NUISANCE_COVARIATES = [
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mt",
    "S_score",
    "G2M_score",
]


def normalize_gene_symbol(gene: str) -> str:
    """Normalize gene symbols for case-insensitive matching."""

    return gene.strip().upper()


def match_prior_genes(
    var_names: list[str] | tuple[str, ...],
    priors: Mapping[str, list[str]] = ORGAN_TROPISM_PRIORS,
) -> dict[str, list[str]]:
    """Return prior genes that are present in an expression matrix."""

    lookup = {normalize_gene_symbol(g): g for g in var_names}
    matched: dict[str, list[str]] = {}
    for name, genes in priors.items():
        present = [lookup[normalize_gene_symbol(g)] for g in genes if normalize_gene_symbol(g) in lookup]
        matched[name] = present
    return matched

