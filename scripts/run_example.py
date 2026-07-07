from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from motmic import MOTMIC, evaluate_binary_labels, simulate_metastasis_dataset
from motmic.interpret import rank_genes_with_shap
from motmic.plotting import plot_site_scores


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    (root / "results").mkdir(exist_ok=True)
    (root / "figures").mkdir(exist_ok=True)
    (root / "data" / "example").mkdir(parents=True, exist_ok=True)

    primary, metastases, truth = simulate_metastasis_dataset()
    primary.to_csv(root / "data" / "example" / "primary_expression.csv")
    truth.to_csv(root / "data" / "example" / "primary_truth.csv")
    for site, expr in metastases.items():
        expr.to_csv(root / "data" / "example" / f"{site}_metastasis_expression.csv")

    model = MOTMIC(n_components=15, epsilon=0.08, rho=1.2, top_k=1)
    result = model.fit_predict(primary, metastases)
    result_frame = result.to_frame()
    result_frame.to_csv(root / "results" / "motmic_scores.csv")

    metrics = evaluate_binary_labels(result.pan_score, truth["true_MIC"])
    pd.Series(metrics).to_csv(root / "results" / "lineage_like_validation_metrics.csv")

    high_mic = result.pan_score >= result.pan_score.quantile(0.8)
    ranking = rank_genes_with_shap(primary, high_mic.astype(int), task="classification")
    ranking.to_csv(root / "results" / "shap_gene_ranking.csv", index=False)

    plot_site_scores(result_frame, str(root / "figures" / "example_site_scores.png"))
    print("Wrote results to", root / "results")
    print(pd.Series(metrics).to_string())
    print(ranking.head(15).to_string(index=False))


if __name__ == "__main__":
    main()
