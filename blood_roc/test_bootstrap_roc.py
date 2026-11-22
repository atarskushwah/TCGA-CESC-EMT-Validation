"""Synthetic-only unit test; no participant data are stored in GitHub."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location("roc", Path(__file__).with_name("bootstrap_roc.py"))
roc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(roc)


def test_synthetic_data():
    rng = np.random.default_rng(42)
    y = np.r_[np.zeros(28, dtype=int), np.ones(62, dtype=int)]
    x = rng.standard_normal((90, 4))
    x[y == 1, 3] -= 1.0
    df = pd.DataFrame(x, columns=roc.PREDICTORS)
    df.insert(0, "CODE", y)
    result = roc.analyze(df, repeats=12, seed=5)
    assert result["sample_n"] == 90
    assert result["cases_n"] == 62
    assert 0 <= result["apparent_auc"] <= 1
    assert 0 <= result["fixed_panel_optimism_corrected_auc"] <= 1


if __name__ == "__main__":
    test_synthetic_data()
    print("Synthetic smoke test passed")
