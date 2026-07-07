# scMIC Algorithm

## Inputs

- `P`: primary tumor-cell expression matrix, cells by genes.
- `M_k`: metastatic tumor-cell expression matrix for organ `k`.
- Optional metadata:
  - clone ID and aggressive clone label,
  - patient ID,
  - metastatic site,
  - spatial coordinates,
  - survival and treatment annotations.

## Step 1: Tumor-Cell Selection

For human tissue samples, restrict analysis to malignant tumor cells using:

- epithelial markers: `EPCAM`, `KRT8`, `KRT18`, `KRT19`;
- CNV inference: CopyKAT, inferCNV, SPATA2 for spatial data;
- malignancy classifiers such as scCancer2/XGBoost when available.

## Step 2: scTour Latent Representation

Concatenate primary and metastatic tumor cells over shared genes and learn a
dynamic latent representation with scTour. This replaces the PCA/AutoEncoder
representation used by the reference scMIC manuscript while keeping the same
biological assumption: metastasis-initiating primary cells should lie closer to
metastatic-state dynamics than ordinary primary tumor cells.

Recommended first-pass gene universe:

- highly variable genes;
- primary tumor subcluster DEGs;
- EMT/hybrid EMT, TGF-beta, hypoxia, apoptosis, glycolysis, MAPK/ERK,
  JAK/STAT3, PI3K/AKT/mTOR, WNT, angiogenesis, ECM-receptor, and focal
  adhesion genes.

The main primary-cell MIC score is the normalized scTour pseudotime / metastatic
state axis among primary tumor cells:

```text
sctour_MIC_score_i = minmax(scTour_time_i)
```

In GSE173958 M1, this score recovers the dominant aggressive lineage group with
AUROC = 0.744 and top-20% Fisher enrichment OR = 4.96.

## Step 3: Organ-Specific UOT

For each metastatic site `k`, compute a cost matrix between primary and metastatic
scTour latent profiles:

```text
C_k(i,j) = distance(z_primary_i, z_metastasis_k_j)
```

Run entropic unbalanced optimal transport:

```text
T_k = UOT(P, M_k; C_k, epsilon, rho)
```

Then retain only the strongest `top_k` origins for each metastatic cell:

```text
T_filtered_k[:,j] = top_k(T_k[:,j])
```

## Step 4: MIC Scoring

OT is used primarily for organotropic propensity and primary-to-metastasis
mapping. In the current GSE173958 M1 run, OT transport mass alone was not a
reliable pan-MIC score, so it is reported separately and used as an auxiliary
mapping signal.

```text
site_MIC_score_i,k = sum_j T_filtered_k(i,j)
transport_pan_score_i = sum_k site_MIC_score_i,k
organ_specificity_i = max_k(site_MIC_score_i,k) / transport_pan_score_i
state_transport_MIC_score_i =
    0.7 * sctour_MIC_score_i + 0.3 * (1 - minmax(transport_pan_score_i))
```

## Step 5: Validation

Validation should be hierarchical:

1. Lineage gold standard: `GSE173958`.
2. Time-course weak labels: `GSE249057/GSE249058`.
3. Human paired metastasis: `GSE178318`.
4. Spatial validation: `GSE277783`.
5. Bulk survival: TCGA.

## Step 6: SHAP-Based Gene Prioritization

Train a classifier or regressor on primary tumor cells:

```text
X = primary-cell gene expression
y = high_sctour_MIC_label or continuous sctour_MIC_score
```

Then compute SHAP:

```text
importance_g = mean(abs(SHAP_g))
```

Recommended final gene score:

```text
0.35 * SHAP
+ 0.25 * MIC/non-MIC DE
+ 0.20 * lineage enrichment
+ 0.10 * cross-dataset reproducibility
+ 0.10 * survival association
```

Candidate genes to track in first-pass analyses:

```text
CD44, TACSTD2, S100A14, RHOD, S100A4, S100A8, S100A9, OCIAD2,
AXL, L1CAM, VIM, FN1, ZEB1, ZEB2, SNAI1, SNAI2, TWIST1, ITGA5,
ITGB1, CXCR4
```
