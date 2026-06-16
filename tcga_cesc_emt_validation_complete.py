#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
TCGA-CESC external validation of EMT-related gene-expression panel.

Panel:
CD44, CTNNB1, MMP7, MMP9, VCL, TIMP1, TIMP2

Analyses:
- Expression landscape
- Spearman correlation heatmap
- Composite EMT score
- Kaplan-Meier survival for EMT score and all panel genes
- Univariate Cox proportional hazards analysis
- Manuscript-ready tables and figures
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import chi2
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

warnings.filterwarnings("ignore")

GENES = ["CD44", "CTNNB1", "MMP7", "MMP9", "VCL", "TIMP1", "TIMP2"]

ENSEMBL_TO_SYMBOL = {
    "ENSG00000026508": "CD44",
    "ENSG00000168036": "CTNNB1",
    "ENSG00000137673": "MMP7",
    "ENSG00000100985": "MMP9",
    "ENSG00000035403": "VCL",
    "ENSG00000102265": "TIMP1",
    "ENSG00000035862": "TIMP2",
}


def read_table(path):
    path = Path(path)
    if path.suffix.lower() in [".tsv", ".txt"]:
        return pd.read_csv(path, sep="\t", low_memory=False)
    return pd.read_csv(path, low_memory=False)


def clean_tcga_id(x, length=15):
    if pd.isna(x):
        return np.nan
    return str(x).strip().replace(".", "-")[:length]


def extract_expression_panel(expr):
    expr = expr.copy()
    first_col = expr.columns[0]
    first_values = expr[first_col].astype(str)

    if first_values.str.startswith("ENSG").mean() > 0.5:
        expr[first_col] = first_values.str.replace(r"\.\d+$", "", regex=True)
        expr = expr.set_index(first_col)
    else:
        possible_gene_cols = [
            c for c in expr.columns
            if c.lower() in ["gene", "genes", "gene_symbol", "symbol", "hugo_symbol",
                             "id", "ensembl_id", "ensembl_gene_id"]
        ]
        if possible_gene_cols:
            gene_col = possible_gene_cols[0]
            expr[gene_col] = expr[gene_col].astype(str).str.replace(r"\.\d+$", "", regex=True)
            expr = expr.set_index(gene_col)

    expr.index = expr.index.astype(str).str.replace(r"\.\d+$", "", regex=True)

    if set(GENES).issubset(expr.index):
        panel = expr.loc[GENES].T
        panel.columns = GENES
    elif all(k in expr.index for k in ENSEMBL_TO_SYMBOL):
        panel = expr.loc[list(ENSEMBL_TO_SYMBOL.keys())].T
        panel.columns = [ENSEMBL_TO_SYMBOL[k] for k in ENSEMBL_TO_SYMBOL]
        panel = panel[GENES]
    elif set(GENES).issubset(expr.columns):
        panel = expr[GENES].copy()
        sample_cols = [c for c in expr.columns if "sample" in c.lower() or "barcode" in c.lower()]
        if sample_cols:
            panel.index = expr[sample_cols[0]]
    else:
        raise ValueError("Could not detect the EMT gene panel in the expression matrix.")

    panel = panel.apply(pd.to_numeric, errors="coerce")

    if np.nanmax(panel.values) > 50:
        panel = np.log2(panel + 1)

    panel.index = panel.index.astype(str).str.replace(".", "-", regex=False).str[:15]
    return panel


def detect_clinical_columns(clin):
    id_candidates = ["Sample ID", "Patient ID", "PATIENT_ID", "SAMPLE_ID",
                     "case_submitter_id", "submitter_id", "bcr_patient_barcode"]
    time_candidates = ["Overall Survival (Months)", "OS_MONTHS", "OS.time",
                       "days_to_death", "days_to_last_follow_up", "days_to_last_followup"]
    status_candidates = ["Overall Survival Status", "OS_STATUS", "OS", "vital_status"]

    id_col = next((c for c in id_candidates if c in clin.columns), None)
    time_col = next((c for c in time_candidates if c in clin.columns), None)
    status_col = next((c for c in status_candidates if c in clin.columns), None)

    if id_col is None:
        tmp = [c for c in clin.columns if "sample" in c.lower() or "patient" in c.lower() or "barcode" in c.lower()]
        id_col = tmp[0] if tmp else None

    if time_col is None:
        tmp = [c for c in clin.columns if "survival" in c.lower() and "month" in c.lower()]
        time_col = tmp[0] if tmp else None

    if status_col is None:
        tmp = [c for c in clin.columns if "survival status" in c.lower() or "vital" in c.lower()]
        status_col = tmp[0] if tmp else None

    if id_col is None or time_col is None or status_col is None:
        raise ValueError("Could not identify clinical ID, OS time, or OS status columns.")

    return id_col, time_col, status_col


def parse_event_status(series):
    s = series.astype(str).str.upper()
    return s.str.contains("DECEASED|DEAD|1:").astype(int)


def km_curve(times, events):
    df = pd.DataFrame({"time": times, "event": events}).dropna().sort_values("time")
    event_times = np.sort(df.loc[df["event"] == 1, "time"].unique())

    tvals = [0.0]
    surv = [1.0]
    s = 1.0

    for t in event_times:
        at_risk = (df["time"] >= t).sum()
        d = ((df["time"] == t) & (df["event"] == 1)).sum()
        if at_risk > 0:
            s *= (1 - d / at_risk)
        tvals.append(t)
        surv.append(s)

    return np.array(tvals), np.array(surv)


def logrank_p(df, group_col):
    a = df[df[group_col] == "High"]
    b = df[df[group_col] == "Low"]
    times = np.sort(np.unique(df.loc[df["event"] == 1, "time"].values))

    O1 = E1 = V1 = 0.0

    for t in times:
        n1 = (a["time"] >= t).sum()
        n2 = (b["time"] >= t).sum()
        d1 = ((a["time"] == t) & (a["event"] == 1)).sum()
        d2 = ((b["time"] == t) & (b["event"] == 1)).sum()

        n = n1 + n2
        d = d1 + d2

        if n <= 1:
            continue

        e1 = d * n1 / n
        v1 = n1 * n2 * d * (n - d) / (n * n * (n - 1)) if n > 1 else 0

        O1 += d1
        E1 += e1
        V1 += v1

    if V1 <= 0:
        return np.nan

    chi_square = (O1 - E1) ** 2 / V1
    return 1 - chi2.cdf(chi_square, 1)


def cox_univariate_breslow(df, variable):
    d = df[["time", "event", variable]].dropna().sort_values("time")
    x = d[variable].values.astype(float)
    t = d["time"].values.astype(float)
    e = d["event"].values.astype(int)

    xs = (x - np.nanmean(x)) / np.nanstd(x)
    beta = 0.0

    for _ in range(100):
        expbx = np.exp(beta * xs)
        score = 0.0
        info = 0.0

        for ti in np.unique(t[e == 1]):
            idx_event = (t == ti) & (e == 1)
            idx_risk = t >= ti
            d_count = idx_event.sum()

            sum_exp = expbx[idx_risk].sum()
            mean_x = (expbx[idx_risk] * xs[idx_risk]).sum() / sum_exp
            mean_x2 = (expbx[idx_risk] * xs[idx_risk] ** 2).sum() / sum_exp

            score += xs[idx_event].sum() - d_count * mean_x
            info += d_count * (mean_x2 - mean_x ** 2)

        if info <= 0:
            return np.nan, np.nan, np.nan, np.nan

        step = score / info
        beta += step

        if abs(step) < 1e-8:
            break

    se = np.sqrt(1 / info)
    z = beta / se
    p = 2 * (1 - stats.norm.cdf(abs(z)))

    return np.exp(beta), np.exp(beta - 1.96 * se), np.exp(beta + 1.96 * se), p


def plot_expression_landscape(data, output):
    fig, ax = plt.subplots(figsize=(9, 5.2))
    values = [data[g].dropna().values for g in GENES]

    ax.violinplot(values, showmedians=True, showextrema=True)
    ax.boxplot(values, widths=0.16, patch_artist=True, showfliers=False,
               boxprops={"facecolor": "white", "alpha": 0.65})

    ax.set_xticks(np.arange(1, len(GENES) + 1))
    ax.set_xticklabels(GENES, rotation=35, ha="right")
    ax.set_ylabel("log2-transformed normalized expression")
    ax.set_title("Expression landscape of EMT-related genes in TCGA-CESC")
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)


def plot_correlation_heatmap(corr, output):
    dist = 1 - corr
    np.fill_diagonal(dist.values, 0)

    try:
        order = leaves_list(linkage(squareform(dist.values, checks=False), method="average"))
    except Exception:
        order = np.arange(len(corr))

    ordered = corr.iloc[order, order]

    fig, ax = plt.subplots(figsize=(6.7, 5.8))
    im = ax.imshow(ordered.values, vmin=-1, vmax=1)

    ax.set_xticks(np.arange(len(ordered)))
    ax.set_yticks(np.arange(len(ordered)))
    ax.set_xticklabels(ordered.columns, rotation=45, ha="right")
    ax.set_yticklabels(ordered.index)

    for i in range(len(ordered)):
        for j in range(len(ordered)):
            ax.text(j, i, f"{ordered.iloc[i, j]:.2f}", ha="center", va="center", fontsize=8)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Spearman rho")
    ax.set_title("Clustered EMT-panel Spearman correlation")

    fig.tight_layout()
    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)


def plot_all_gene_km(data, results, output_png, output_pdf):
    variables = ["EMT_score"] + GENES

    fig, axes = plt.subplots(2, 4, figsize=(18, 9), sharey=True)
    axes = axes.flatten()

    for ax, var in zip(axes, variables):
        tmp = data[["time", "event", var]].dropna().copy()
        cutoff = tmp[var].median()
        tmp["group"] = np.where(tmp[var] >= cutoff, "High", "Low")

        row = results[results["Variable"] == var].iloc[0]

        for group in ["Low", "High"]:
            d = tmp[tmp["group"] == group]
            t, s = km_curve(d["time"], d["event"])
            ax.step(t, s, where="post", label=f"{group} ({len(d)})", linewidth=1.9)

        ax.set_title(var, fontsize=12, fontweight="bold")
        ax.set_xlabel("Overall survival (months)")
        ax.set_ylim(0, 1.05)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(alpha=0.18, linewidth=0.5)
        ax.legend(frameon=False, fontsize=8, loc="lower left")

        stat_text = (
            f"log-rank p = {row['Log_rank_p']:.3f}\n"
            f"HR/SD = {row['Cox_HR_per_1SD']:.2f} "
            f"({row['Cox_95CI_lower']:.2f}-{row['Cox_95CI_upper']:.2f})"
        )

        ax.text(
            0.04, 0.20, stat_text,
            transform=ax.transAxes,
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="0.65", alpha=0.95)
        )

    axes[0].set_ylabel("Overall survival probability")
    axes[4].set_ylabel("Overall survival probability")
    fig.suptitle("Kaplan-Meier overall survival analysis for EMT-panel genes in TCGA-CESC",
                 fontsize=15, fontweight="bold", y=1.02)

    fig.tight_layout()
    fig.savefig(output_png, dpi=600, bbox_inches="tight")
    fig.savefig(output_pdf, bbox_inches="tight")
    plt.close(fig)


def run_analysis(expression_file, clinical_file, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    expr = read_table(expression_file)
    clin = read_table(clinical_file)

    panel = extract_expression_panel(expr)

    id_col, time_col, status_col = detect_clinical_columns(clin)

    clin["sample_short"] = clin[id_col].apply(clean_tcga_id)
    clin["time"] = pd.to_numeric(clin[time_col], errors="coerce")
    clin["event"] = parse_event_status(clin[status_col])

    merged = panel.merge(
        clin[["sample_short", "time", "event"]],
        left_index=True,
        right_on="sample_short",
        how="inner"
    )

    merged = merged.dropna(subset=["time", "event"]).drop_duplicates(subset=["sample_short"])
    merged = merged.dropna(subset=GENES)

    zpanel = merged[GENES].apply(lambda x: stats.zscore(x, nan_policy="omit"))
    merged["EMT_score"] = zpanel.mean(axis=1)
    merged["EMT_group"] = np.where(
        merged["EMT_score"] >= merged["EMT_score"].median(),
        "High EMT",
        "Low EMT"
    )

    expression_summary = merged[GENES].describe().T[["mean", "50%", "std", "min", "max"]]
    expression_summary.columns = ["Mean", "Median", "SD", "Minimum", "Maximum"]
    expression_summary.to_csv(output_dir / "Table_1_expression_summary.csv")

    corr = merged[GENES].corr(method="spearman")
    corr.to_csv(output_dir / "Table_2_spearman_correlation_matrix.csv")

    pmat = pd.DataFrame(index=GENES, columns=GENES, dtype=float)
    for g1 in GENES:
        for g2 in GENES:
            _, p = stats.spearmanr(merged[g1], merged[g2], nan_policy="omit")
            pmat.loc[g1, g2] = p
    pmat.to_csv(output_dir / "Table_3_spearman_pvalue_matrix.csv")

    variables = ["EMT_score"] + GENES
    survival_rows = []

    for var in variables:
        tmp = merged[["sample_short", "time", "event", var]].dropna().copy()
        cutoff = tmp[var].median()
        tmp["group"] = np.where(tmp[var] >= cutoff, "High", "Low")

        p_logrank = logrank_p(tmp, "group")
        hr, lo, hi, p_cox = cox_univariate_breslow(tmp, var)

        survival_rows.append({
            "Variable": var,
            "Cutoff": "Median",
            "Cutoff_value": cutoff,
            "Low_n": int((tmp["group"] == "Low").sum()),
            "High_n": int((tmp["group"] == "High").sum()),
            "Events_total": int(tmp["event"].sum()),
            "Log_rank_p": p_logrank,
            "Cox_HR_per_1SD": hr,
            "Cox_95CI_lower": lo,
            "Cox_95CI_upper": hi,
            "Cox_p": p_cox,
        })

    survival_results = pd.DataFrame(survival_rows)
    survival_results.to_csv(output_dir / "Table_4_all_gene_KM_Cox_summary.csv", index=False)

    merged.to_csv(output_dir / "TCGA_CESC_EMT_merged_analysis_dataset.csv", index=False)

    plot_expression_landscape(merged, output_dir / "Figure_1_expression_landscape.png")
    plot_correlation_heatmap(corr, output_dir / "Figure_2_spearman_clustered_heatmap.png")
    plot_all_gene_km(
        merged,
        survival_results,
        output_dir / "Figure_3_all_gene_KM_survival_panels.png",
        output_dir / "Figure_3_all_gene_KM_survival_panels.pdf"
    )

    lines = []
    for _, r in survival_results.iterrows():
        sig = "statistically significant" if r["Log_rank_p"] < 0.05 else ("trend-level" if r["Log_rank_p"] < 0.10 else "non-significant")
        lines.append(
            f"{r['Variable']}: log-rank p={r['Log_rank_p']:.3f}; "
            f"Cox HR per 1 SD={r['Cox_HR_per_1SD']:.2f} "
            f"(95% CI {r['Cox_95CI_lower']:.2f}-{r['Cox_95CI_upper']:.2f}), "
            f"Cox p={r['Cox_p']:.3f}; {sig}."
        )

    interpretation = (
        "TCGA-CESC all-gene Kaplan-Meier survival analysis\n"
        "=================================================\n\n"
        f"Matched cohort size: {merged.shape[0]}\n"
        f"OS events: {int(merged['event'].sum())}\n"
        "Survival endpoint: Overall survival in months\n"
        "Stratification: Median cutoff for each gene/score\n\n"
        + "\n".join(lines)
        + "\n\nManuscript interpretation:\n"
        "Median-stratified Kaplan-Meier analysis was performed for the composite EMT score and each EMT-panel gene. "
        "The overall pattern indicates limited prognostic separation for most individual EMT transcripts in TCGA-CESC bulk tumor data. "
        "This suggests that baseline tumor EMT-related transcript abundance may capture tumor-state heterogeneity, adhesion remodeling, "
        "stromal-matrix interaction and epithelial-mesenchymal plasticity more strongly than stable overall-survival determinants. "
        "This is important because the primary clinical study measured peripheral blood expression in patients receiving cisplatin-based chemoradiotherapy, "
        "where EMT-related signals may reflect dynamic systemic and therapy-induced biology that is not fully represented in static pretreatment TCGA tumor transcriptomes.\n"
    )

    (output_dir / "All_gene_KM_interpretation.txt").write_text(interpretation, encoding="utf-8")

    summary = (
        "TCGA-CESC EMT validation completed.\n"
        f"Expression file: {expression_file}\n"
        f"Clinical file: {clinical_file}\n"
        f"Clinical ID column: {id_col}\n"
        f"Survival time column: {time_col}\n"
        f"Survival status column: {status_col}\n"
        f"Matched cohort size: {merged.shape[0]}\n"
        f"OS events: {int(merged['event'].sum())}\n"
    )
    (output_dir / "analysis_summary.txt").write_text(summary, encoding="utf-8")

    return merged, expression_summary, corr, survival_results


def main():
    parser = argparse.ArgumentParser(description="TCGA-CESC EMT-panel external validation workflow.")
    parser.add_argument("--expression", required=True, help="Path to TCGA-CESC normalized RNA-seq expression matrix.")
    parser.add_argument("--clinical", required=True, help="Path to TCGA-CESC clinical/survival metadata.")
    parser.add_argument("--output", default="tcga_cesc_emt_outputs", help="Output directory.")

    args = parser.parse_args()
    run_analysis(args.expression, args.clinical, args.output)


if __name__ == "__main__":
    main()
