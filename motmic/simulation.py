"""Synthetic paired primary/metastasis data for examples and tests."""

from __future__ import annotations

import numpy as np
import pandas as pd


CORE_GENES = [
    "CD44",
    "TACSTD2",
    "S100A14",
    "RHOD",
    "S100A4",
    "S100A8",
    "S100A9",
    "OCIAD2",
    "AXL",
    "L1CAM",
    "VIM",
    "FN1",
    "ZEB1",
    "ZEB2",
    "SNAI1",
    "SNAI2",
    "TWIST1",
    "ITGA5",
    "ITGB1",
    "CXCR4",
]


def simulate_metastasis_dataset(
    n_primary: int = 240,
    n_metastasis_per_site: int = 90,
    n_noise_genes: int = 80,
    random_state: int = 7,
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame], pd.DataFrame]:
    """Create a toy dataset with known MIC and organotropic labels."""

    rng = np.random.default_rng(random_state)
    genes = CORE_GENES + [f"GENE{i:03d}" for i in range(n_noise_genes)]
    sites = ["lung", "liver", "bone"]

    primary = rng.normal(0, 0.7, size=(n_primary, len(genes)))
    labels = np.array(["nonMIC"] * n_primary, dtype=object)
    site_label = np.array(["none"] * n_primary, dtype=object)

    mic_idx = rng.choice(n_primary, size=int(n_primary * 0.22), replace=False)
    split = np.array_split(mic_idx, len(sites))
    for site, idx in zip(sites, split):
        labels[idx] = "MIC"
        site_label[idx] = site

    gene_to_idx = {g: i for i, g in enumerate(genes)}
    pan_genes = ["CD44", "TACSTD2", "S100A14", "RHOD", "VIM", "FN1", "AXL"]
    liver_genes = ["S100A8", "S100A9", "OCIAD2"]
    lung_genes = ["L1CAM", "ITGA5", "ITGB1"]
    bone_genes = ["CXCR4", "ZEB1", "ZEB2", "TWIST1"]
    site_genes = {"lung": lung_genes, "liver": liver_genes, "bone": bone_genes}

    for idx in mic_idx:
        primary[idx, [gene_to_idx[g] for g in pan_genes]] += rng.normal(2.0, 0.25, len(pan_genes))
        genes_for_site = site_genes[site_label[idx]]
        primary[idx, [gene_to_idx[g] for g in genes_for_site]] += rng.normal(1.8, 0.25, len(genes_for_site))

    primary = np.log1p(np.exp(primary))
    primary_expr = pd.DataFrame(primary, index=[f"P{i:03d}" for i in range(n_primary)], columns=genes)

    metastases = {}
    for site in sites:
        origin_idx = np.where(site_label == site)[0]
        sampled = rng.choice(origin_idx, size=n_metastasis_per_site, replace=True)
        expr = primary[sampled] + rng.normal(0.05, 0.25, size=(n_metastasis_per_site, len(genes)))
        genes_for_site = site_genes[site]
        expr[:, [gene_to_idx[g] for g in genes_for_site]] += rng.normal(
            0.7, 0.15, (n_metastasis_per_site, len(genes_for_site))
        )
        expr = np.clip(expr, 0, None)
        metastases[site] = pd.DataFrame(expr, index=[f"{site}_M{i:03d}" for i in range(n_metastasis_per_site)], columns=genes)

    truth = pd.DataFrame(
        {
            "cell_id": primary_expr.index,
            "true_MIC": labels == "MIC",
            "true_site": site_label,
            "aggressive_clone": np.where(labels == "MIC", "aggressive", "non_aggressive"),
        }
    ).set_index("cell_id")
    return primary_expr, metastases, truth
