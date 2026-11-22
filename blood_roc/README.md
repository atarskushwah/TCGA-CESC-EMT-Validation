# Exploratory whole-blood four-transcript ROC: reproducibility protocol

**Status:** reconstructed internal assessment, **not** a verified copy of the original SPSS/Prism model, and **not** clinical validation.

This procedure supplements, but is analytically distinct from, the tumor-tissue TCGA-CESC contextualization scripts elsewhere in the repository. **Do not mix participant-level blood and TCGA tumor data.**

## Inputs (local only)

Export the `Case_Control` sheet of the original `CODED_DATA.xlsx` study workbook to a CSV with columns `CODE,CD44,CTNNB1,MMP9,TIMP2`, retaining the original participant row order. There should be 62 cases (CODE=1) and 28 HPV-negative healthy controls (CODE=0), each with four log2-fold-change measurements. **Never commit participant-level or potentially identifying records.**

## Reproduce the fixed-panel bootstrap

```bash
python -m pip install numpy pandas scipy scikit-learn
python blood_roc/bootstrap_roc.py --input /private/path/emt_case_control.csv --output /private/path/bootstrap_roc.json --repeats 1000 --seed 20261122
```

Reconstructed output with the uploaded workbook at the time of analysis:

| Quantity | Estimate |
|---|---:|
| Cases / controls | 62 / 28 |
| Full-cohort apparent AUC | 0.86924 |
| Bootstrap resamples | 1,000 (stratified with replacement) |
| Average AUC optimism | 0.02785 |
| Fixed-panel corrected AUC | 0.84139 |

The four predictors are **CD44, CTNNB1, MMP9, TIMP2**. A near-unpenalized logistic regression (`sklearn.LogisticRegression`, `C=1e6`, `solver=lbfgs`, `max_iter=1000`) is fitted to the complete sample and then refitted in each bootstrap sample, using the same four genes. In bootstrap sample *b*, fitting optimism equals AUC in the resample minus AUC on the original cohort using the model fitted in that resample. Mean optimism is subtracted from the apparent full-cohort AUC.

**Limits:** The four genes were originally chosen after testing all seven genes in this cohort. The gene selection step is held fixed during bootstrapping, so selection optimism is not corrected. These values are not estimates of real-world screening performance, since the comparator group is HPV-negative healthy controls. The original 80.6% sensitivity and 75.0% specificity are not reproduced by the reconstructed model's Youden-optimal operating point (77.42% and 85.71%); do not interpret the latter as a correction of the original before inspecting its original model probabilities and threshold rules. Do not claim a bootstrap confidence interval from the distribution of optimism.

A separate exploratory analysis with Welch/BH gene screening produced a different value but did **not** reproduce the documented original feature-selection protocol and is therefore not included as a validated result.

## GitHub Actions

The `blood_roc.yml` workflow runs a **synthetic-data smoke test** only. It tests that the pipeline executes and returns plausible values. It cannot confirm the manuscript's numerical results because the individual-level clinical data are not public. Those numbers were obtained by a local run and should be independently verified in an authorized research environment before journal submission.

## TCGA-CESC distinction

The TCGA-CESC cohort provides independent **tumor-tissue transcriptomic contextualization**. It is not external validation of this peripheral whole-blood diagnostic classifier, and its transcript abundance cannot establish functional EMT activity.
