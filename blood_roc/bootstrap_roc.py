#!/usr/bin/env python3
"""Reproduce the fixed-panel bootstrap ROC analysis for the exploratory EMT blood cohort.

No clinical/qPCR source data are distributed with this public repository.
The user supplies a local CSV with columns CODE, CD44, CTNNB1, MMP9, TIMP2.
CODE: 1 = cervical cancer, 0 = HPV-negative healthy control.
All four predictors must be participant-level log2 fold changes.

This is a *reconstructed*, near-unpenalized sklearn logistic model, not
a verified replication of the original SPSS/Prism model. Predictors were
selected on the full cohort before resampling, so the corrected AUC does
not account for predictor-selection optimism.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve

PREDICTORS = ["CD44", "CTNNB1", "MMP9", "TIMP2"]


def fit(X, y):
    # C=1e6 gives a near-unpenalized fit; it is not identical to unpenalized MLE.
    model = LogisticRegression(C=1e6, solver="lbfgs", max_iter=1000)
    return model.fit(X, y)


def validate_data(df):
    cols = ["CODE"] + PREDICTORS
    missing = [col for col in cols if col not in df.columns]
    if missing:
        raise ValueError("Missing required columns: " + ", ".join(missing))
    data = df[cols].copy().dropna(how="all")
    if data.isna().any().any():
        raise ValueError("Incomplete participant data: no imputation was used.")
    y = data.CODE.to_numpy(dtype=int)
    if not np.isin(y, [0, 1]).all():
        raise ValueError("CODE must be 0 for control or 1 for cancer.")
    X = data[PREDICTORS].to_numpy(dtype=float)
    if not np.isfinite(X).all():
        raise ValueError("Nonfinite expression values.")
    if (y == 1).sum() == 0 or (y == 0).sum() == 0:
        raise ValueError("Both case and control groups are required.")
    return X, y


def analyze(df, repeats=1000, seed=20261122):
    X, y = validate_data(df)
    model = fit(X, y)
    predicted = model.predict_proba(X)[:, 1]
    apparent = float(roc_auc_score(y, predicted))
    fpr, tpr, thresholds = roc_curve(y, predicted)
    best = int(np.argmax(tpr - fpr))
    rng = np.random.default_rng(seed)
    group_ids = [np.where(y == group)[0] for group in (0, 1)]
    optimism = []
    for _ in range(repeats):
        idx = np.concatenate([rng.choice(g, size=len(g), replace=True) for g in group_ids])
        bootstrap_model = fit(X[idx], y[idx])
        auc_training = roc_auc_score(y[idx], bootstrap_model.predict_proba(X[idx])[:, 1])
        auc_original = roc_auc_score(y, bootstrap_model.predict_proba(X)[:, 1])
        optimism.append(float(auc_training - auc_original))
    mean_optimism = float(np.mean(optimism))
    return {
        "sample_n": int(len(y)),
        "cases_n": int((y == 1).sum()),
        "controls_n": int((y == 0).sum()),
        "outcome_coding": "CODE: 1 cancer; 0 HPV-negative healthy control",
        "predictors": PREDICTORS,
        "algorithm": "sklearn LogisticRegression(lbfgs,C=1000000,max_iter=1000)",
        "apparent_auc": apparent,
        "original_dataset_youden_sensitivity": float(tpr[best]),
        "original_dataset_youden_specificity": float(1 - fpr[best]),
        "original_dataset_youden_probability_cutoff": float(thresholds[best]),
        "coefficients": {g: float(c) for g, c in zip(PREDICTORS, model.coef_[0])},
        "intercept": float(model.intercept_[0]),
        "bootstrap_repetitions": repeats,
        "stratified": True,
        "seed": seed,
        "mean_fitting_optimism": mean_optimism,
        "fixed_panel_optimism_corrected_auc": apparent - mean_optimism,
        "validation_scope": "fitting optimism conditional on full-cohort-selected fixed panel; not selection-adjusted, external or clinical validation"
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="Local patient-level CSV (never commit)")
    p.add_argument("--output", required=True, help="Local output JSON file")
    p.add_argument("--repeats", type=int, default=1000)
    p.add_argument("--seed", type=int, default=20261122)
    args = p.parse_args()
    if args.repeats < 1:
        p.error("--repeats must be positive")
    result = analyze(pd.read_csv(args.input), args.repeats, args.seed)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print("Apparent AUC:", round(result["apparent_auc"], 5))
    print("Mean fitting optimism:", round(result["mean_fitting_optimism"], 5))
    print("Fixed-panel corrected AUC:", round(result["fixed_panel_optimism_corrected_auc"], 5))


if __name__ == "__main__":
    main()
