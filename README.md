# MOT-MIC

[![Docs](https://github.com/yuhan1li/scMIC/actions/workflows/docs.yml/badge.svg)](https://github.com/yuhan1li/scMIC/actions/workflows/docs.yml)

<img src="figures/motmic_algorithm_schematic.png" width="520px" align="left">

MOT-MIC is a computational framework for discovering metastasis-initiating
cells from paired primary and metastatic single-cell transcriptomes.

It uses a multi-organ optimal-transport architecture to map primary tumor cells
to metastatic tumor cells across different secondary sites, then assigns each
primary cell both a pan-metastatic score and organ-specific metastatic
propensity scores.

MOT-MIC is designed for method development and validation on public metastasis
datasets, including lineage-traced, time-course, paired human, spatial, and
bulk survival cohorts.

<br clear="left"/>

## Key features

- Multi-organ metastasis-initiating cell scoring from paired primary and
  metastatic scRNA-seq data.
- Unbalanced optimal transport with sparse top-k origin filtering for
  interpretable primary-to-metastasis mapping.
- Organotropic MIC scores for liver, lung, bone, brain, or user-defined
  metastatic sites.
- Lineage-aware validation workflow using `GSE173958`, the strongest available
  public benchmark for metastatic clone aggression.
- scTour-style tutorial notebooks for each example dataset.
- SHAP-based gene prioritization for pan-MIC and organ-specific metastatic
  programs.
- Spatial and bulk-cohort validation templates for MIC signature transfer.

## Installation

Clone the repository and install the required Python packages:

```console
git clone https://github.com/yuhan1li/scMIC.git
cd scMIC
pip install -r requirements.txt
```

Run the lightweight example:

```console
python scripts/run_example.py
python scripts/make_diagram.py
```

## Basic usage

```python
from motmic import MOTMIC, simulate_metastasis_dataset
from motmic.interpret import rank_genes_with_shap

primary, metastases, truth = simulate_metastasis_dataset()

model = MOTMIC(n_components=15, epsilon=0.08, rho=1.2, top_k=1)
result = model.fit_predict(primary, metastases)

scores = result.to_frame()
high_mic = result.pan_score >= result.pan_score.quantile(0.8)
gene_ranking = rank_genes_with_shap(primary, high_mic.astype(int))
```

## Tutorials

The tutorials follow the same teaching rhythm as the scTour basic inference
notebooks:

```text
Dataset -> Data loading -> Model training -> Inference -> Visualization -> Validation -> Robustness
```

| Notebook | Dataset | Analysis shown |
|---|---|---|
| [MOT-MIC quick start](https://scmic.readthedocs.io/en/latest/notebook/01_motmic_example.html) | Synthetic paired data | Fully runnable MOT-MIC quick start with SHAP ranking |
| [GSE173958 lineage validation](https://scmic.readthedocs.io/en/latest/notebook/02_GSE173958_lineage_validation.html) | GSE173958 | Primary-to-liver/lung/macrometastasis inference with lineage validation plan |
| [GSE249057 time-course discovery](https://scmic.readthedocs.io/en/latest/notebook/03_GSE249057_timecourse_discovery.html) | GSE249057 | 0h parental tumor to 6h/2mo/4mo metastatic-state discovery |
| [GSE178318 human CRC liver validation](https://scmic.readthedocs.io/en/latest/notebook/04_GSE178318_human_crc_liver_metastasis.html) | GSE178318 | Human CRC primary-to-liver validation and patient-level aggregation |
| [GSE277783 spatial validation](https://scmic.readthedocs.io/en/latest/notebook/05_GSE277783_spatial_validation.html) | GSE277783 | Spatial validation of MOT-MIC-derived MIC gene programs |

## Recommended datasets

| Dataset | Role | Gold-standard level | Notes |
|---|---|---:|---|
| `GSE173958` | Main validation | Strongest | scRNA-seq plus CRISPR lineage tracing in metastatic PDAC; aggressive clones provide the closest cell-level ground truth. |
| `GSE249057` / `GSE249058` | Discovery/training | Weak | Multi-timepoint ESCC lung metastasis model; useful for learning early-to-late metastatic state transitions. |
| `GSE178318` | Human paired validation | Weak clinical | Matched primary CRC and liver metastasis scRNA-seq; validate MIC burden against liver metastatic burden. |
| `GSE277783` | Spatial validation | Weak spatial | PDAC spatial data; validate MIC/spots in spatial niches and treatment-associated regions. |
| `OMIX002487` | External human validation | Weak clinical | Human PDAC scRNA-seq used in the scMIC manuscript. |
| TCGA | Bulk survival validation | Cohort-level only | Test whether derived MIC signatures predict OS/DFS. |

Inspect available GEO files without downloading large matrices:

```console
python scripts/download_geo.py --dry-run --from-filelist --gse GSE173958 GSE249057
```

## Documentation

- Algorithm details: `ALGORITHM.md`
- Dataset plan: `DATASETS.md`
- Notebook index: `notebooks/README.md`
- Sphinx-style documentation skeleton: `docs/source/index.rst`
- Online documentation: https://scmic.readthedocs.io/en/latest/

## Method summary

For each metastatic site `k`, MOT-MIC estimates a transport plan from primary
tumor cells to metastatic tumor cells:

```text
T_k = UOT(primary_cells, metastatic_cells_k)
```

Sparse top-k filtering is applied to each metastatic cell column, and each
primary cell receives:

```text
site_MIC_score_i,k = sum_j T_filtered_k(i,j)
pan_MIC_score_i = sum_k site_MIC_score_i,k
organ_specificity_i = max_k(site_MIC_score_i,k) / pan_MIC_score_i
```

Candidate metastatic genes are prioritized with SHAP and then cross-checked by
MIC/non-MIC differential expression, lineage enrichment, cross-dataset
reproducibility, and survival association.

## Reference datasets and literature

MOT-MIC is motivated by recent work on scMIC, macsGESTALT lineage tracing,
MetaNet organotropic risk modeling, SIDISH clinical single-cell integration, and
single-cell trajectory/transport methods including scTour and optimal transport
models.
