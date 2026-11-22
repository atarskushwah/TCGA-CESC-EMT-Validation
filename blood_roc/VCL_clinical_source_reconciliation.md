# Source-data reconciliation protocol: VCL normalization and clinical outcomes

This document separates source-data observations from previously reported manuscript results. **No individual-level clinical or qPCR records are included in this public repository.**

## VCL sign convention (22 November 2026 source audit)

In the local final qPCR workbook, the VCL (vinculin) sheet records fold change as `2^(+DeltaDeltaCt)` and `log2FC = +DeltaDeltaCt`; the six other genes use `2^(-DeltaDeltaCt)` and `log2FC = -DeltaDeltaCt`. To represent all seven targets on the same stated normalization scale:

```python
corrected_vcl_log2fc = -source_vcl_log2fc
corrected_vcl_foldchange = 2.0 ** corrected_vcl_log2fc
```

Among 62 cases, source VCL grouping gives 41 positive/21 negative observations, whereas the corrected grouping gives 21 positive/41 negative. Two-sided case–control P values do not change under global sign reversal, but correlations reverse sign and continuous Cox HRs become reciprocals (with inverse/reordered CI bounds). The four-transcript ROC panel excludes VCL and therefore is not numerically altered.

## Clinical survival findings

The separately supplied local 62-case SPSS follow-up dataset contains 28 PFS events, nine OS deaths, and an observed follow-up median of 20 months (range 11–29). The adjusted Cox source analysis used a binary FIGO stage code (early I–II, n=33; advanced III–IV, n=29) with age and parity. The MMP9 adjusted PFS HR is approximately 0.880 (95% CI 0.798–0.971; P=0.011), reproducible with a Breslow-tied Cox fit. Exact time origin, censoring calendar date, and RECIST assessment timing still require site records.

The survival categories in the submitted manuscript do not all match the stated `log2FC > 0` / `log2FC < 0` definition in the uploaded follow-up file. In particular, TIMP2 log-rank estimates differ from the historical figure. Do not claim the Kaplan–Meier analyses are reconciled until their original SPSS grouping derivation and input syntax are reviewed.

## ROC threshold

The uploaded original Youden-index calculation identifies maximum index about 0.631 at predicted-probability threshold about 0.673525, with apparent sensitivity 48/62 = 77.4% and specificity 24/28 = 85.7%. This conflicts with historical annotations of 80.6%/75.0%. Figure operating-point labels require revision or a documented alternate threshold. This is separate from the previously documented fixed-panel bootstrap-corrected AUC.

## Confidentiality

This protocol describes only formulas and cohort-level aggregates. Raw `.sav` files, individual-level observations, and workbooks are **not** committed or uploaded to GitHub. The existing GitHub Actions test uses synthetic data and is not evidence that privately held patient-level results have been independently reproduced.
