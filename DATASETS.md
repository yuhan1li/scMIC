# Dataset Plan

## Primary Gold Standard

### GSE173958

Use this as the main validation dataset because it contains both single-cell
transcriptomes and CRISPR lineage information from a metastatic PDAC model.
Aggressive clones and subclonal dissemination can be used as the nearest
available gold standard for MIC predictions.

The GEO file list includes `GSE173958_RAW.tar` and separate 10x files such as
`M1-PT`, `M1-Liver`, `M1-Lung`, `M1-CTC`, `M2-PT`, `M2-Met`, `M2-Liver`, and
`M2-Lung`. This makes it suitable for primary-to-organ-specific metastatic
mapping.

Completed M1 validation run:

- Downloaded `GSE173958_RAW.tar` on the server and analyzed M1 primary tumor,
  M1-Met, M1-Liver, and M1-Lung.
- Parsed macsGESTALT `*.stats.txt.gz` files to recover lineage clone labels from
  10x cell barcodes.
- Dominant metastatic lineage group:
  `1I+48+C|105D+67|105D+67|105D+67|105D+67|NONE|NONE|NONE|NONE|NONE`.
- `sctour_MIC_score` versus dominant aggressive lineage in primary cells:
  AUROC `0.744`, AUPRC `0.117`, top-20% enrichment OR `4.96`,
  Fisher P `2.74e-18`.
- OT transport mass alone was not sufficient as a pan-MIC score in this run
  (AUROC `0.448`), so the current framework uses scTour for the main MIC score
  and OT for organotropic mapping.
- Trajectory visualization should prioritize the 1D scTour pseudotime axis,
  sample-state distributions, and lineage enrichment along pseudotime. The 2D
  scTour latent scatter is only an auxiliary diagnostic and should not be read
  as a complete cancer progression path.

Validation tasks:

- AUROC/AUPRC for `sctour_MIC_score` against aggressive clone labels.
- Enrichment of top 5%, 10%, and 20% predicted MICs in aggressive clones.
- Correlation between organ-specific scores and observed metastatic harvest
  sites.
- Agreement with late-hybrid EMT / pseudoEMT states.

## Discovery And Training

### GSE249057 / GSE249058

Use as a discovery dataset because it contains multi-timepoint ESCC metastasis
states. The MIC-like cluster defined in the associated scMIC manuscript should be
treated as a weak label, not a strict gold standard.

The GEO file list includes `GSE249057_RAW.tar` and 10x matrices for `0h`,
`6h`, `48h`, `2mo`, and `4mo`. The `0h` cells serve as primary/parental tumor
cells; late metastatic cells are useful as UOT targets.

Validation tasks:

- Recover the MIC-like primary cluster.
- Test sensitivity to embedding method, UOT parameters, and top-k filtering.
- Check enrichment of CD44, TACSTD2, S100A14, RHOD, and EMT/hybrid EMT
  programs.

## Human Paired Validation

### GSE178318

Use matched colorectal cancer primary and liver metastasis samples to validate
human transferability.

Validation tasks:

- Compute primary-cell liver MIC scores.
- Aggregate per patient and compare with liver metastasis number/size when
  available.
- Compare treatment-naive and treated samples if metadata supports it.

## Spatial Validation

### GSE277783

Use PDAC spatial transcriptomic data to validate whether MIC signatures are
spatially enriched in metastatic or treatment-associated tumor regions.

Validation tasks:

- Score malignant spots with the MOT-MIC-derived MIC signature.
- Test enrichment near invasive, hypoxic, stromal, or metastatic niches.
- Compare treatment-associated changes in MIC-positive spots.

## Cohort-Level Validation

### TCGA

Use TCGA only for bulk-level clinical relevance, not cell-level MIC validation.

Validation tasks:

- Derive a compact MIC signature from SHAP-ranked and lineage-validated genes.
- Stratify patients by signature score.
- Test OS/DFS association with Kaplan-Meier and Cox regression.

## Download Commands

Inspect available files without downloading large matrices:

```bash
python scripts/download_geo.py --dry-run --from-filelist --gse GSE173958 GSE249057
```

Download listed supplementary files:

```bash
python scripts/download_geo.py --from-filelist --gse GSE173958 GSE249057
```

Large archives should stay under `data/raw/`, which is ignored by git.
