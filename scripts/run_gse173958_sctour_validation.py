"""Run the GSE173958 scTour + OT lineage-validation workflow.

This analysis follows the MIC paper's logic while replacing PCA/AutoEncoder
representation with scTour latent space:

1. read paired primary and metastatic M1 10x matrices,
2. train scTour on tumor-cell expression snapshots,
3. compute primary-to-metastasis UOT plans in scTour latent space,
4. score primary cells with top-1 transport mass,
5. validate high-scoring cells against lineage clones inferred from macsGESTALT
   stats files.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from collections import Counter
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import shap
from scipy import io, sparse, stats
from scipy.spatial.distance import cdist
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.preprocessing import StandardScaler

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from motmic.ot import topk_filter, unbalanced_sinkhorn


SAMPLES = {
    "primary": "GSM5283482_M1-PT",
    "met": "GSM5283484_M1-Met",
    "liver": "GSM5283485_M1-Liver",
    "lung": "GSM5283486_M1-Lung",
}

STATS_FILES = {
    "primary": "GSM5283494_M1-PT.stats.txt.gz",
    "met": "GSM5283496_M1-Met.stats.txt.gz",
    "liver": "GSM5283497_M1-Liver.stats.txt.gz",
    "lung": "GSM5283498_M1-Lung.stats.txt.gz",
}

EPITHELIAL_GENES = ["Epcam", "Krt8", "Krt18", "Krt19", "Cdh1"]
MESENCHYMAL_GENES = ["Vim", "Zeb1", "Zeb2", "Fn1", "Sparc", "Col1a1", "S100a4"]


def read_10x_prefixed(raw_dir: Path, sample_prefix: str) -> ad.AnnData:
    matrix = io.mmread(raw_dir / f"{sample_prefix}-matrix.mtx.gz").tocsr().T
    features = pd.read_csv(raw_dir / f"{sample_prefix}-features.tsv.gz", sep="\t", header=None)
    with gzip.open(raw_dir / f"{sample_prefix}-barcodes.tsv.gz", "rt") as handle:
        barcodes = [line.strip() for line in handle if line.strip()]
    gene_ids = features.iloc[:, 0].astype(str).values
    gene_symbols = features.iloc[:, 1].astype(str).values if features.shape[1] > 1 else gene_ids
    obs = pd.DataFrame(index=barcodes)
    var = pd.DataFrame({"gene_ids": gene_ids}, index=gene_symbols)
    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    adata.var_names_make_unique()
    return adata


def load_m1_data(raw_dir: Path, max_cells_per_sample: int, random_state: int) -> ad.AnnData:
    adatas = []
    rng = np.random.default_rng(random_state)
    for sample, prefix in SAMPLES.items():
        a = read_10x_prefixed(raw_dir, prefix)
        a.obs["sample"] = sample
        a.obs_names = [f"{sample}:{bc}" for bc in a.obs_names]
        if max_cells_per_sample and a.n_obs > max_cells_per_sample:
            keep = rng.choice(a.n_obs, size=max_cells_per_sample, replace=False)
            a = a[sorted(keep)].copy()
        adatas.append(a)
    combined = ad.concat(adatas, join="inner", label=None, index_unique=None)
    combined.var_names_make_unique()
    sc.pp.calculate_qc_metrics(combined, inplace=True, percent_top=None, log1p=False)
    combined = combined[combined.obs["n_genes_by_counts"] >= 300].copy()
    sc.pp.filter_genes(combined, min_cells=10)
    sc.pp.highly_variable_genes(combined, n_top_genes=2000, flavor="cell_ranger", subset=True)
    return combined


def parse_cell_barcode(read_name: str) -> str | None:
    parts = read_name.split("_")
    if len(parts) < 2:
        return None
    match = re.match(r"([ACGTN]{16})", parts[1])
    if not match:
        return None
    return f"{match.group(1)}-1"


def infer_cell_lineages(stats_path: Path) -> pd.Series:
    target_cols = [f"target{i}" for i in range(1, 11)]
    usecols = ["readName", "keep", *target_cols]
    votes: dict[str, Counter[str]] = {}
    for chunk in pd.read_csv(stats_path, sep="\t", usecols=usecols, chunksize=200_000):
        chunk = chunk[chunk["keep"].eq("PASS")]
        for row in chunk.itertuples(index=False):
            barcode = parse_cell_barcode(row.readName)
            if barcode is None:
                continue
            targets = [str(getattr(row, c)) for c in target_cols]
            key = "|".join(targets)
            votes.setdefault(barcode, Counter())[key] += 1
    consensus = {bc: counts.most_common(1)[0][0] for bc, counts in votes.items() if counts}
    return pd.Series(consensus, name="lineage_clone")


def attach_lineage_labels(adata: ad.AnnData, raw_dir: Path) -> tuple[ad.AnnData, str]:
    clone_tables = {}
    for sample, filename in STATS_FILES.items():
        clone_tables[sample] = infer_cell_lineages(raw_dir / filename)

    lineage = []
    for cell_id, sample in zip(adata.obs_names, adata.obs["sample"]):
        barcode = cell_id.split(":", 1)[1]
        lineage.append(clone_tables[str(sample)].get(barcode, np.nan))
    adata.obs["lineage_clone"] = lineage

    metastatic_clones = []
    for sample in ["met", "liver", "lung"]:
        metastatic_clones.extend(clone_tables[sample].dropna().tolist())
    dominant_clone = Counter(metastatic_clones).most_common(1)[0][0]
    adata.obs["aggressive_clone"] = adata.obs["lineage_clone"].eq(dominant_clone)
    return adata, dominant_clone


def train_sctour_embedding(adata: ad.AnnData, n_latent: int, epochs: int, random_state: int) -> np.ndarray:
    import sctour as sct

    if sparse.issparse(adata.X):
        adata.X = adata.X.toarray()
    adata.X = np.asarray(adata.X, dtype=np.float32)

    trainer = sct.train.Trainer(
        adata,
        n_latent=n_latent,
        loss_mode="nb",
        nepoch=epochs,
        batch_size=1024,
        random_state=random_state,
        use_gpu=False,
    )
    trainer.train()
    adata.obs["sctour_pseudotime"] = trainer.get_time()
    mix_zs, _, _ = trainer.get_latentsp(batch_size=1024)
    adata.obsm["X_sctour"] = np.asarray(mix_zs)
    return adata.obsm["X_sctour"]


def scaled_pairwise_cost(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    cost = cdist(x, y, metric="sqeuclidean")
    scale = np.median(cost[cost > 0]) if np.any(cost > 0) else 1.0
    return cost / max(scale, 1e-12)


def run_transport_scores(adata: ad.AnnData, epsilon: float, rho: float, top_k: int) -> pd.DataFrame:
    z = StandardScaler().fit_transform(adata.obsm["X_sctour"])
    sample = adata.obs["sample"].astype(str)
    primary_mask = sample.eq("primary").values
    primary_ids = adata.obs_names[primary_mask]
    primary_z = z[primary_mask]

    site_scores = {}
    for site in ["met", "liver", "lung"]:
        site_mask = sample.eq(site).values
        cost = scaled_pairwise_cost(primary_z, z[site_mask])
        plan = unbalanced_sinkhorn(cost, epsilon=epsilon, rho=rho, max_iter=500, tol=1e-7)
        filtered = topk_filter(plan, k=top_k)
        raw_score = filtered.sum(axis=1)
        site_scores[site] = raw_score / max(raw_score.max(), 1e-12)

    scores = pd.DataFrame(site_scores, index=primary_ids)
    scores["transport_pan_score"] = scores[["met", "liver", "lung"]].sum(axis=1)
    scores["transport_pan_score"] = scores["transport_pan_score"] / max(scores["transport_pan_score"].max(), 1e-12)
    scores["predicted_site"] = scores[["met", "liver", "lung"]].idxmax(axis=1)
    return scores


def score_gene_set(adata: ad.AnnData, genes: list[str]) -> np.ndarray:
    genes = [g for g in genes if g in adata.var_names]
    if not genes:
        return np.zeros(adata.n_obs)
    x = adata[:, genes].X
    x = x.toarray() if sparse.issparse(x) else x
    return np.asarray(np.log1p(x).mean(axis=1)).ravel()


def minmax(values: pd.Series | np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    return (x - np.nanmin(x)) / (np.nanmax(x) - np.nanmin(x) + 1e-12)


def evaluate_primary(adata: ad.AnnData, scores: pd.DataFrame, outdir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    primary = adata[adata.obs["sample"].eq("primary")].copy()
    scores = scores.loc[primary.obs_names].copy()
    scores["lineage_clone"] = primary.obs["lineage_clone"].values
    scores["aggressive_clone"] = primary.obs["aggressive_clone"].fillna(False).astype(bool).values
    scores["sctour_pseudotime"] = primary.obs["sctour_pseudotime"].values
    scores["epithelial_score"] = score_gene_set(primary, EPITHELIAL_GENES)
    scores["mesenchymal_score"] = score_gene_set(primary, MESENCHYMAL_GENES)
    scores["hybrid_emt_score"] = np.minimum(scores["epithelial_score"], scores["mesenchymal_score"])
    scores["sctour_MIC_score"] = minmax(scores["sctour_pseudotime"])
    scores["transport_inverse_score"] = 1.0 - minmax(scores["transport_pan_score"])
    scores["state_transport_MIC_score"] = (
        0.7 * scores["sctour_MIC_score"] + 0.3 * scores["transport_inverse_score"]
    )

    valid = scores["lineage_clone"].notna()
    y = scores.loc[valid, "aggressive_clone"].astype(int)
    pred = scores.loc[valid, "sctour_MIC_score"]
    threshold = scores["sctour_MIC_score"].quantile(0.8)
    high = scores.loc[valid, "sctour_MIC_score"].ge(threshold)
    table = pd.crosstab(high, y).reindex(index=[False, True], columns=[0, 1], fill_value=0)
    odds, pvalue = stats.fisher_exact(table.values)
    transport_pred = scores.loc[valid, "transport_pan_score"]
    state_transport_pred = scores.loc[valid, "state_transport_MIC_score"]
    metrics = pd.DataFrame(
        [
            {
                "dataset": "GSE173958_M1",
                "primary_cells": len(scores),
                "lineage_labeled_primary_cells": int(valid.sum()),
                "aggressive_clone_primary_cells": int(y.sum()),
                "top20_percent_threshold": float(threshold),
                "sctour_MIC_AUROC": float(roc_auc_score(y, pred)) if y.nunique() > 1 else np.nan,
                "sctour_MIC_AUPRC": float(average_precision_score(y, pred)) if y.nunique() > 1 else np.nan,
                "transport_mass_AUROC": float(roc_auc_score(y, transport_pred)) if y.nunique() > 1 else np.nan,
                "state_transport_MIC_AUROC": (
                    float(roc_auc_score(y, state_transport_pred)) if y.nunique() > 1 else np.nan
                ),
                "fisher_odds_ratio_top20_vs_aggressive_clone": float(odds),
                "fisher_pvalue_top20_vs_aggressive_clone": float(pvalue),
            }
        ]
    )

    outdir.mkdir(parents=True, exist_ok=True)
    scores.to_csv(outdir / "GSE173958_m1_sctour_ot_mic_scores.csv")
    metrics.to_csv(outdir / "GSE173958_lineage_validation_metrics.csv", index=False)
    return scores, metrics


def rank_genes_with_shap(adata: ad.AnnData, scores: pd.DataFrame, outdir: Path, random_state: int) -> pd.DataFrame:
    primary = adata[adata.obs["sample"].eq("primary")].copy()
    y = scores.loc[primary.obs_names, "sctour_MIC_score"].ge(scores["sctour_MIC_score"].quantile(0.8)).astype(int)
    x = primary.X
    x = x.toarray() if sparse.issparse(x) else x
    rng = np.random.default_rng(random_state)
    if x.shape[0] > 2500:
        idx = np.sort(rng.choice(x.shape[0], size=2500, replace=False))
        x = x[idx]
        y = y.iloc[idx]
    model = RandomForestClassifier(n_estimators=300, max_depth=6, random_state=random_state, n_jobs=-1)
    model.fit(x, y)
    explainer = shap.TreeExplainer(model)
    values = explainer.shap_values(x)
    if isinstance(values, list):
        values = values[1]
    if values.ndim == 3:
        values = values[:, :, 1]
    ranking = pd.DataFrame(
        {
            "gene": primary.var_names,
            "mean_abs_shap": np.abs(values).mean(axis=0),
            "rf_importance": model.feature_importances_,
        }
    ).sort_values("mean_abs_shap", ascending=False)
    ranking.to_csv(outdir / "GSE173958_shap_gene_ranking.csv", index=False)
    return ranking


def make_figures(adata: ad.AnnData, scores: pd.DataFrame, ranking: pd.DataFrame, metrics: pd.DataFrame, figdir: Path) -> None:
    figdir.mkdir(parents=True, exist_ok=True)
    primary = adata[adata.obs["sample"].eq("primary")].copy()
    coords = primary.obsm["X_sctour"][:, :2]
    score = scores.loc[primary.obs_names, "sctour_MIC_score"]

    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    sca = ax.scatter(coords[:, 0], coords[:, 1], c=score, s=8, cmap="viridis", linewidths=0)
    ax.set_xlabel("scTour latent 1")
    ax.set_ylabel("scTour latent 2")
    ax.set_title("GSE173958 M1 primary scTour-MIC score")
    fig.colorbar(sca, ax=ax, label="scTour-MIC score")
    fig.tight_layout()
    fig.savefig(figdir / "GSE173958_sctour_latent_mic_score.png", dpi=220)
    plt.close(fig)

    threshold = scores["sctour_MIC_score"].quantile(0.8)
    plot_df = scores.assign(high_MIC=scores["sctour_MIC_score"].ge(threshold))
    enrich = plot_df.groupby("high_MIC")["aggressive_clone"].mean().reindex([False, True])
    fig, ax = plt.subplots(figsize=(4.2, 3.6))
    ax.bar(["other primary", "top 20% MIC"], enrich.values, color=["#8aa1b4", "#ca6f57"])
    ax.set_ylabel("Aggressive-clone fraction")
    ax.set_ylim(0, max(0.05, float(enrich.max()) * 1.35))
    ax.set_title(f"Lineage validation AUROC={metrics.loc[0, 'sctour_MIC_AUROC']:.3f}")
    fig.tight_layout()
    fig.savefig(figdir / "GSE173958_lineage_enrichment.png", dpi=220)
    plt.close(fig)

    top = ranking.head(20).iloc[::-1]
    fig, ax = plt.subplots(figsize=(5.5, 5.2))
    ax.barh(top["gene"], top["mean_abs_shap"], color="#5f8f7b")
    ax.set_xlabel("Mean |SHAP|")
    ax.set_title("Top MIC-associated genes")
    fig.tight_layout()
    fig.savefig(figdir / "GSE173958_shap_top_genes.png", dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", default="data/raw/GSE173958")
    parser.add_argument("--outdir", default="data/processed")
    parser.add_argument("--figdir", default="figures")
    parser.add_argument("--max-cells-per-sample", type=int, default=4000)
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--n-latent", type=int, default=10)
    parser.add_argument("--epsilon", type=float, default=0.08)
    parser.add_argument("--rho", type=float, default=1.2)
    parser.add_argument("--top-k", type=int, default=1)
    parser.add_argument("--random-state", type=int, default=13)
    args = parser.parse_args()

    raw_dir = Path(args.raw_dir)
    outdir = Path(args.outdir)
    figdir = Path(args.figdir)

    adata = load_m1_data(raw_dir, args.max_cells_per_sample, args.random_state)
    adata, dominant_clone = attach_lineage_labels(adata, raw_dir)
    train_sctour_embedding(adata, args.n_latent, args.epochs, args.random_state)
    scores = run_transport_scores(adata, args.epsilon, args.rho, args.top_k)
    scores, metrics = evaluate_primary(adata, scores, outdir)
    ranking = rank_genes_with_shap(adata, scores, outdir, args.random_state)
    make_figures(adata, scores, ranking, metrics, figdir)

    summary = {
        "dominant_metastatic_clone_key": dominant_clone,
        "metrics": metrics.iloc[0].to_dict(),
        "top_shap_genes": ranking.head(20)["gene"].tolist(),
    }
    (outdir / "GSE173958_analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
