import pandas as pd
from sklearn.metrics import roc_curve

from helpers.deg_analysis.classifier_metrics import get_metrics_at_threshold


def test_get_metrics_at_threshold():
    df = pd.DataFrame(
        {
            "score": [0.1, 0.2, 0.3, 0.4, 0.5],
            "label": [0, 0, 1, 1, 1],
        }
    ).set_index("score")
    tpr, fpr, thresholds = roc_curve(
        y_true=df["label"],
        y_score=df["score"],
    )
    df_curves = pd.DataFrame(
        {
            "tpr": tpr,
            "fpr": fpr,
            "threshold": thresholds,
        }
    )
    result = get_metrics_at_threshold(
        df_curves=df_curves,
        groupby=[],
        score="score",
        threshold=0.4,
        metrics=["tpr", "fpr"],
    )
    assert result["tpr"] == 0.5
    assert result["fpr"] == 0.0
