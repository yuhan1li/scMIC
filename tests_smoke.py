from motmic import MOTMIC, evaluate_binary_labels, simulate_metastasis_dataset


def test_smoke():
    primary, metastases, truth = simulate_metastasis_dataset(n_primary=60, n_metastasis_per_site=20, n_noise_genes=10)
    result = MOTMIC(n_components=8, epsilon=0.1).fit_predict(primary, metastases)
    frame = result.to_frame()
    assert frame.shape[0] == primary.shape[0]
    assert "pan_MIC_score" in frame
    metrics = evaluate_binary_labels(result.pan_score, truth["true_MIC"])
    assert "AUROC" in metrics


if __name__ == "__main__":
    test_smoke()
    print("smoke test passed")

