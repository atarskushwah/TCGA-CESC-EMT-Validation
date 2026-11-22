# TCGA-CESC EMT-Panel Tumor-Tissue Contextualization

This repository provides reproducible code for complementary tumor-tissue transcriptomic contextualization of seven EMT-associated genes in TCGA-CESC. The TCGA tumor profiles do not externally validate the separate peripheral whole-blood RT-qPCR assay.

## Gene panel

- CD44
- CTNNB1
- MMP7
- MMP9
- VCL
- TIMP1
- TIMP2

The script also calculates a composite EMT score using the mean of z-score-normalized expression values across the panel.

## Required input files

1. TCGA-CESC normalized RNA-seq matrix  
   Example: `CESC_normalized.csv`

2. TCGA-CESC clinical/survival metadata  
   Example: `cesc_tcga_pan_can_atlas_2018_clinical_data.tsv`

The expression matrix can use Ensembl IDs or gene symbols. Internal Ensembl mapping is included.

## Install dependencies

```bash
pip install -r requirements.txt
```

## Run analysis

```bash
python tcga_cesc_emt_validation_complete.py \
  --expression CESC_normalized.csv \
  --clinical cesc_tcga_pan_can_atlas_2018_clinical_data.tsv \
  --output tcga_cesc_emt_outputs
```

## Outputs

Figures:
- `Figure_1_expression_landscape.png`
- `Figure_2_spearman_clustered_heatmap.png`
- `Figure_3_all_gene_KM_survival_panels.png`
- `Figure_3_all_gene_KM_survival_panels.pdf`

Tables:
- `Table_1_expression_summary.csv`
- `Table_2_spearman_correlation_matrix.csv`
- `Table_3_spearman_pvalue_matrix.csv`
- `Table_4_all_gene_KM_Cox_summary.csv`
- `TCGA_CESC_EMT_merged_analysis_dataset.csv`

Text:
- `analysis_summary.txt`
- `All_gene_KM_interpretation.txt`

## Manuscript-ready methods text

Normalized TCGA-CESC RNA-seq data were integrated with clinical survival metadata from the TCGA PanCancer Atlas cervical cancer cohort. Expression values for CD44, CTNNB1, MMP7, MMP9, VCL, TIMP1 and TIMP2 were extracted using gene symbols or Ensembl identifiers and log2-transformed where required. Samples lacking overall survival time or event information were excluded. Spearman correlation was used to evaluate inter-gene association. A composite EMT score was calculated as the mean of z-score-normalized expression values across the seven-gene panel. Patients were stratified into high- and low-expression groups using the median value for each gene or score. Kaplan-Meier analysis with two-sided log-rank testing was used to compare overall survival, and univariate Cox proportional hazards regression was used to estimate hazard ratios and 95% confidence intervals.

## Interpretation note

TCGA-CESC tumor transcriptomics provides complementary tumor-tissue expression and co-expression context, not functional validation of EMT activity or validation of peripheral whole-blood RNA measurements. Weak or non-significant survival associations should not be overinterpreted. Static baseline tumor EMT expression may reflect tumor heterogeneity, stromal-matrix remodeling and epithelial-mesenchymal plasticity more strongly than dynamic therapy-induced EMT states.

## Data availability wording

The code used for TCGA-CESC tumor-tissue contextualization is available in this GitHub repository. TCGA-CESC RNA-seq and clinical data are publicly available from the Genomic Data Commons, UCSC Xena and cBioPortal TCGA PanCancer Atlas resources. Processed outputs include expression summary statistics, Spearman correlation matrices, Kaplan-Meier survival results, Cox regression tables and publication-ready figures.

## Separate whole-blood ROC internal assessment

The [exploratory four-transcript blood ROC reproducibility protocol](blood_roc/README.md) documents a fixed-panel bootstrap analysis and a separate synthetic-data CI smoke test. No participant-level qPCR or clinical data are distributed in this public repository. This analysis is distinct from TCGA tumor-tissue contextualization and cannot be described as independent clinical validation.
