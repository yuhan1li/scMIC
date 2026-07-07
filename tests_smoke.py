import pandas as pd

from motmic import MOTMIC, evaluate_binary_labels


def test_smoke():
    genes = ["EPCAM", "KRT19", "VIM", "AXL", "ITGB1", "MKI67"]
    primary = pd.DataFrame(
        [
            [5.0, 4.8, 0.3, 0.2, 0.5, 1.0],
            [4.9, 4.6, 0.4, 0.3, 0.6, 1.1],
            [3.7, 3.3, 1.8, 2.1, 2.0, 1.4],
            [3.5, 3.0, 2.0, 2.4, 2.2, 1.3],
            [2.0, 1.9, 4.1, 4.3, 4.0, 0.7],
            [1.8, 1.6, 4.4, 4.5, 4.2, 0.8],
        ],
        index=[f"primary_{i}" for i in range(6)],
        columns=genes,
    )
    metastases = {
        "liver": pd.DataFrame(
            [[3.6, 3.2, 1.9, 2.2, 2.1, 1.3], [3.4, 3.1, 2.1, 2.3, 2.0, 1.5]],
            index=["liver_0", "liver_1"],
            columns=genes,
        ),
        "lung": pd.DataFrame(
            [[1.9, 1.7, 4.2, 4.4, 4.1, 0.7], [1.7, 1.5, 4.5, 4.6, 4.3, 0.8]],
            index=["lung_0", "lung_1"],
            columns=genes,
        ),
    }
    labels = pd.Series([0, 0, 1, 1, 1, 1], index=primary.index)

    result = MOTMIC(n_components=3, epsilon=0.1).fit_predict(primary, metastases)
    frame = result.to_frame()
    assert frame.shape[0] == primary.shape[0]
    assert "pan_MIC_score" in frame
    metrics = evaluate_binary_labels(result.pan_score, labels)
    assert "AUROC" in metrics


if __name__ == "__main__":
    test_smoke()
    print("smoke test passed")
