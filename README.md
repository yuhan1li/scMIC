# MOT-MIC

MOT-MIC, short for Multi-Organ Optimal Transport for Metastasis-Initiating Cells,
is a prototype framework for identifying metastasis-initiating cells (MICs) in
paired primary and metastatic single-cell RNA-seq data.

## Algorithm Overview

![MOT-MIC algorithm schematic](figures/motmic_algorithm_schematic.png)

The method is designed around the literature in `相关文献`:

- scMIC: embedding + unbalanced optimal transport + top-k transport filtering.
- macsGESTALT / GSE173958: lineage tracing as the closest available cell-level
  gold standard for metastatic aggression.
- MetaNet: organotropic risk modeling and SHAP-based feature interpretation.
- SIDISH: linking single-cell states to clinical outcomes and spatial validation.

## What The Algorithm Estimates

For every primary tumor cell, MOT-MIC estimates:

```text
pan_MIC_score
lung_MIC_score
liver_MIC_score
bone_MIC_score
brain_MIC_score
organ_specificity
```

The output is therefore richer than a binary MIC label: it captures both
metastatic competence and metastatic destination bias.

## Recommended Datasets

| Dataset | Role | Gold-standard level | Notes |
|---|---|---:|---|
| GSE173958 | Main validation | Strongest | scRNA-seq plus CRISPR lineage tracing in metastatic PDAC; aggressive clones provide the closest cell-level ground truth. |
| GSE249057 / GSE249058 | Discovery/training | Weak | Multi-timepoint ESCC lung metastasis model; useful for learning early-to-late metastatic state transitions. |
| GSE178318 | Human paired validation | Weak clinical | Matched primary CRC and liver metastasis scRNA-seq; validate MIC burden against liver metastatic burden. |
| GSE277783 | Spatial validation | Weak spatial | PDAC spatial data; validate MIC/spots in spatial niches and treatment-associated regions. |
| OMIX002487 | External human validation | Weak clinical | Human PDAC scRNA-seq used in the scMIC paper. |
| TCGA | Bulk survival validation | Cohort-level only | Test whether derived MIC signatures predict OS/DFS. |

## Repository Layout

```text
motmic/                  Python package
scripts/run_example.py   Reproducible synthetic demonstration
scripts/download_geo.py  GEO supplementary-file downloader
scripts/make_diagram.py  Algorithm schematic
notebooks/               Example notebook
figures/                 Generated schematic
results/                 Example outputs
```

## Quick Start

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/run_example.py
python scripts/make_diagram.py
```

The example uses a synthetic paired primary/multi-metastasis dataset. It is small
enough for GitHub while preserving the expected inputs and outputs.

## Core Workflow

1. Select malignant tumor cells using marker expression and CNV/malignancy scores.
2. Learn a robust latent representation from primary and metastatic cells.
3. For each metastatic organ, run unbalanced optimal transport from primary cells
   to metastatic cells.
4. Apply top-k filtering on transport columns to enforce sparse interpretable
   metastatic origins.
5. Sum filtered transport weights per primary cell to obtain organ-specific MIC
   scores.
6. Validate against lineage, clinical burden, spatial enrichment, and survival.
7. Train an interpretable model and rank genes with SHAP or permutation
   importance.

## SHAP Gene Ranking

MOT-MIC can train a model using primary-cell expression as input and MIC scores or
high/low MIC labels as output. For tree models, SHAP values are summarized as:

```text
importance_gene = mean(abs(SHAP_gene))
```

The recommended final ranking combines:

```text
0.35 * SHAP score
+ 0.25 * MIC vs non-MIC differential expression
+ 0.20 * lineage/aggressive-clone enrichment
+ 0.10 * cross-dataset reproducibility
+ 0.10 * bulk survival association
```

SHAP should be used as a prioritization tool, not as a standalone proof of causal
metastatic function.
