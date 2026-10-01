#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mapping_analysis_full.py

PHẦN A — So sánh thống kê (paired): #1–8
    mapping_rate, properly_paired_rate, mapq_mean

PHẦN B — Thống kê MÔ TẢ (không kiểm định):
    - Coverage (graph: depth, stddev, cv; linear: depth, breadth)
    - Paired plot cho depth (mô tả xu hướng)
    - Bar chart so sánh MAPQ median (thay cho histogram vô nghĩa)
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

# ==================== CẤU HÌNH ====================
INPUT_CSV = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/Mapping_and_variant_assessment/linear_reference/read_mapping_comparison/mapping_comparison.csv"

OUTPUT_DIR = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/Mapping_and_variant_assessment/linear_reference/read_mapping_comparison"

OUTPUT_STATS_CSV = os.path.join(OUTPUT_DIR, "statistical_results.csv")
OUTPUT_DESCRIPTIVE_CSV = os.path.join(OUTPUT_DIR, "descriptive_stats.csv")
OUTPUT_PLOT_DIR = os.path.join(OUTPUT_DIR, "plots")
OUTPUT_DESCRIPTIVE_PLOT_DIR = os.path.join(OUTPUT_DIR, "descriptive_plots")

# --- PHẦN A: các cặp chỉ số so sánh (có kiểm định) ---
COMPARISON_METRICS = [
    ("graph_mapping_rate",         "linear_mapping_rate",         "Mapping rate (%)"),
    ("graph_properly_paired_rate", "linear_properly_paired_rate", "Properly paired (%)"),
    ("graph_mapq_mean",            "linear_mapq_mean",            "MAPQ mean"),
]

# --- PHẦN B: các cột cần thống kê mô tả (bỏ MAPQ median — hằng số) ---
DESCRIPTIVE_COLS = [
    ("graph_mean_depth",            "Graph: mean depth"),
    ("graph_depth_stddev",          "Graph: depth stddev"),
    ("graph_depth_cv",              "Graph: depth CV"),
    ("linear_mean_depth",           "Linear: mean depth"),
    ("linear_breadth_coverage_pct", "Linear: breadth coverage (%)"),
]

# --- PHẦN B: các cặp muốn vẽ paired plot MÔ TẢ (không kiểm định) ---
DESCRIPTIVE_PAIRED_COLS = [
    ("graph_mean_depth", "linear_mean_depth", "Mean depth"),
]

N_BOOTSTRAP = 10000
ALPHA = 0.05
BONFERRONI_ALPHA = ALPHA / len(COMPARISON_METRICS)
RANDOM_SEED = 42
# ==================================================


# ============================================================
# PHẦN A — SO SÁNH THỐNG KÊ (PAIRED)
# ============================================================

def cohens_d_paired(differences):
    return np.mean(differences) / np.std(differences, ddof=1)


def ci_95_paired(differences):
    n = len(differences)
    mean_diff = np.mean(differences)
    se = np.std(differences, ddof=1) / np.sqrt(n)
    t_crit = stats.t.ppf(0.975, n - 1)
    return mean_diff - t_crit * se, mean_diff + t_crit * se


def bootstrap_ci(differences, n_boot=N_BOOTSTRAP, seed=RANDOM_SEED):
    rng = np.random.default_rng(seed)
    n = len(differences)
    boot_means = np.array([
        np.mean(rng.choice(differences, size=n, replace=True))
        for _ in range(n_boot)
    ])
    return np.percentile(boot_means, [2.5, 97.5])


def leave_one_out(graph, linear):
    n = len(graph)
    p_values = []
    for i in range(n):
        idx = [j for j in range(n) if j != i]
        _, p = stats.ttest_rel(graph[idx], linear[idx])
        p_values.append(p)
    return np.min(p_values), np.max(p_values), p_values


def analyze_metric(graph, linear, metric_name):
    graph = np.asarray(graph, dtype=float)
    linear = np.asarray(linear, dtype=float)
    diff = graph - linear
    n = len(diff)

    results = {
        "metric": metric_name,
        "n_samples": n,
        "graph_mean": round(np.mean(graph), 3),
        "graph_sd": round(np.std(graph, ddof=1), 3),
        "linear_mean": round(np.mean(linear), 3),
        "linear_sd": round(np.std(linear, ddof=1), 3),
        "mean_diff": round(np.mean(diff), 3),
        "sd_diff": round(np.std(diff, ddof=1), 3),
    }

    shapiro_stat, shapiro_p = stats.shapiro(diff)
    results["shapiro_W"] = round(shapiro_stat, 4)
    results["shapiro_p"] = round(shapiro_p, 5)
    results["is_normal"] = "Yes" if shapiro_p > 0.05 else "No"

    t_stat, t_p = stats.ttest_rel(graph, linear)
    results["t_stat"] = round(t_stat, 4)
    results["t_p"] = t_p

    try:
        w_stat, w_p = stats.wilcoxon(graph, linear)
        results["wilcoxon_stat"] = round(w_stat, 2)
        results["wilcoxon_p"] = w_p
    except ValueError:
        results["wilcoxon_stat"] = None
        results["wilcoxon_p"] = None

    results["t_p_bonferroni"] = "significant" if t_p < BONFERRONI_ALPHA else "not significant"
    if results["wilcoxon_p"] is not None:
        results["wilcoxon_p_bonferroni"] = (
            "significant" if results["wilcoxon_p"] < BONFERRONI_ALPHA else "not significant"
        )

    d = cohens_d_paired(diff)
    results["cohens_d"] = round(d, 4)
    if abs(d) < 0.2:
        results["effect_size"] = "negligible"
    elif abs(d) < 0.5:
        results["effect_size"] = "small"
    elif abs(d) < 0.8:
        results["effect_size"] = "medium"
    else:
        results["effect_size"] = "large"

    ci_low, ci_high = ci_95_paired(diff)
    results["ci95_low"] = round(ci_low, 3)
    results["ci95_high"] = round(ci_high, 3)

    boot_low, boot_high = bootstrap_ci(diff)
    results["bootstrap_ci_low"] = round(boot_low, 3)
    results["bootstrap_ci_high"] = round(boot_high, 3)

    loo_min, loo_max, _ = leave_one_out(graph, linear)
    results["loo_p_min"] = round(loo_min, 5)
    results["loo_p_max"] = round(loo_max, 5)
    results["loo_robust"] = "Yes" if loo_max < BONFERRONI_ALPHA else "No"

    return results


def plot_comparison_metric(graph, linear, metric_name, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    for g, l in zip(graph, linear):
        ax.plot([0, 1], [g, l], marker="o", color="gray", alpha=0.5, linewidth=0.8)
    ax.plot([0, 1], [np.mean(graph), np.mean(linear)],
            marker="s", color="red", linewidth=2.5, markersize=10, label="Mean")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Graph", "Linear"])
    ax.set_ylabel(metric_name)
    ax.set_title(f"Paired comparison\n{metric_name}")
    ax.legend()
    ax.grid(alpha=0.3)

    ax = axes[1]
    bp = ax.boxplot([graph, linear], tick_labels=["Graph", "Linear"],
                    patch_artist=True, widths=0.5)
    bp["boxes"][0].set_facecolor("#4C72B0")
    bp["boxes"][1].set_facecolor("#DD8452")
    ax.scatter([1, 2], [np.mean(graph), np.mean(linear)],
               marker="D", color="red", s=60, zorder=3, label="Mean")
    ax.set_ylabel(metric_name)
    ax.set_title(f"Distribution\n{metric_name}")
    ax.legend()
    ax.grid(alpha=0.3, axis="y")

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"📊 Đã lưu biểu đồ: {out_path}")


def run_comparison_analysis(df):
    print("\n" + "=" * 70)
    print("PHẦN A — SO SÁNH THỐNG KÊ (PAIRED)")
    print("=" * 70)

    os.makedirs(OUTPUT_PLOT_DIR, exist_ok=True)
    all_results = []

    for graph_col, linear_col, metric_name in COMPARISON_METRICS:
        if graph_col not in df.columns or linear_col not in df.columns:
            print(f"⚠️ Bỏ qua {metric_name} — thiếu cột")
            continue

        print(f"\n{'='*70}")
        print(f"📈 {metric_name}")
        print(f"{'='*70}")

        sub = df[[graph_col, linear_col]].dropna()
        graph = sub[graph_col].values
        linear = sub[linear_col].values

        if len(graph) < 2:
            print("⚠️ Không đủ mẫu — bỏ qua")
            continue

        results = analyze_metric(graph, linear, metric_name)
        all_results.append(results)

        print(f"  Graph:  {results['graph_mean']} ± {results['graph_sd']}")
        print(f"  Linear: {results['linear_mean']} ± {results['linear_sd']}")
        print(f"  Chênh lệch: {results['mean_diff']} ± {results['sd_diff']}")
        print(f"\n  #1 Shapiro-Wilk:      W={results['shapiro_W']}, p={results['shapiro_p']} "
              f"→ {'Chuẩn' if results['is_normal']=='Yes' else 'Không chuẩn'}")
        print(f"  #2 Paired t-test:     t={results['t_stat']}, p={results['t_p']:.6f}")
        print(f"  #3 Wilcoxon:          W={results['wilcoxon_stat']}, p={results['wilcoxon_p']:.6f}")
        print(f"  #4 Bonferroni (α={BONFERRONI_ALPHA:.4f}):")
        print(f"       t-test:     {results['t_p_bonferroni']}")
        print(f"       Wilcoxon:   {results.get('wilcoxon_p_bonferroni', 'N/A')}")
        print(f"  #5 Cohen's d:         {results['cohens_d']} ({results['effect_size']})")
        print(f"  #6 95% CI:            [{results['ci95_low']}, {results['ci95_high']}]")
        print(f"  #7 Bootstrap 95% CI:  [{results['bootstrap_ci_low']}, {results['bootstrap_ci_high']}]")
        print(f"  #8 Leave-one-out:     p ∈ [{results['loo_p_min']}, {results['loo_p_max']}] "
              f"→ {'Robust' if results['loo_robust']=='Yes' else 'Không robust'}")

        safe_name = metric_name.replace(" ", "_").replace("(%)", "").replace("(", "").replace(")", "")
        plot_path = os.path.join(OUTPUT_PLOT_DIR, f"{safe_name}.png")
        plot_comparison_metric(graph, linear, metric_name, plot_path)

    if all_results:
        df_results = pd.DataFrame(all_results)
        df_results.to_csv(OUTPUT_STATS_CSV, index=False)
        print(f"\n🎉 Đã lưu kết quả thống kê: {OUTPUT_STATS_CSV}")


# ============================================================
# PHẦN B — THỐNG KÊ MÔ TẢ
# ============================================================

def describe_series(values, label):
    values = np.asarray(values, dtype=float)
    values = values[~np.isnan(values)]
    n = len(values)

    if n == 0:
        return {"metric": label, "n": 0}

    return {
        "metric": label,
        "n": n,
        "mean": round(np.mean(values), 3),
        "sd": round(np.std(values, ddof=1), 3) if n > 1 else 0.0,
        "min": round(np.min(values), 3),
        "q1": round(np.percentile(values, 25), 3),
        "median": round(np.median(values), 3),
        "q3": round(np.percentile(values, 75), 3),
        "max": round(np.max(values), 3),
        "iqr": round(np.percentile(values, 75) - np.percentile(values, 25), 3),
    }


def plot_distribution(values, label, out_path):
    values = np.asarray(values, dtype=float)
    values = values[~np.isnan(values)]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    ax = axes[0]
    ax.hist(values, bins=min(10, len(values)), color="#4C72B0",
            edgecolor="black", alpha=0.8)
    ax.axvline(np.mean(values), color="red", linestyle="--",
               linewidth=2, label=f"Mean = {np.mean(values):.2f}")
    ax.axvline(np.median(values), color="green", linestyle=":",
               linewidth=2, label=f"Median = {np.median(values):.2f}")
    ax.set_xlabel(label)
    ax.set_ylabel("Số mẫu")
    ax.set_title(f"Phân bố\n{label}")
    ax.legend()
    ax.grid(alpha=0.3)

    ax = axes[1]
    bp = ax.boxplot([values], tick_labels=[""], patch_artist=True, widths=0.4)
    bp["boxes"][0].set_facecolor("#4C72B0")
    ax.scatter([1], [np.mean(values)], marker="D", color="red",
               s=60, zorder=3, label=f"Mean = {np.mean(values):.2f}")
    ax.set_ylabel(label)
    ax.set_title(f"Boxplot\n{label}")
    ax.legend()
    ax.grid(alpha=0.3, axis="y")

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"📊 Đã lưu: {out_path}")


def plot_descriptive_paired(graph, linear, metric_name, out_path):
    """Vẽ paired line plot + boxplot MÔ TẢ (không kiểm định)."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    for g, l in zip(graph, linear):
        ax.plot([0, 1], [g, l], marker="o", color="gray", alpha=0.5, linewidth=0.8)
    ax.plot([0, 1], [np.mean(graph), np.mean(linear)],
            marker="s", color="red", linewidth=2.5, markersize=10, label="Mean")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Graph", "Linear"])
    ax.set_ylabel(metric_name)
    ax.set_title(f"Paired comparison (descriptive)\n{metric_name}\n"
                 f"⚠️ Không kiểm định — không gian đo khác nhau")
    ax.legend()
    ax.grid(alpha=0.3)

    ax = axes[1]
    bp = ax.boxplot([graph, linear], tick_labels=["Graph", "Linear"],
                    patch_artist=True, widths=0.5)
    bp["boxes"][0].set_facecolor("#4C72B0")
    bp["boxes"][1].set_facecolor("#DD8452")
    ax.scatter([1, 2], [np.mean(graph), np.mean(linear)],
               marker="D", color="red", s=60, zorder=3, label="Mean")
    ax.set_ylabel(metric_name)
    ax.set_title(f"Distribution\n{metric_name}")
    ax.legend()
    ax.grid(alpha=0.3, axis="y")

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"📊 Đã lưu: {out_path}")


def plot_mapq_median_bar(df, out_path):
    """
    Bar chart so sánh MAPQ median giữa graph và linear.
    Vì cả hai đều = 60 (hằng số) → histogram vô nghĩa, dùng bar chart.
    """
    graph_vals = df["graph_mapq_median"].dropna().values
    linear_vals = df["linear_mapq_median"].dropna().values

    if len(graph_vals) == 0 or len(linear_vals) == 0:
        print("⚠️ Không đủ dữ liệu MAPQ median để vẽ")
        return

    fig, ax = plt.subplots(figsize=(7, 5))

    x = [0, 1]
    means = [np.mean(graph_vals), np.mean(linear_vals)]
    sds = [np.std(graph_vals, ddof=1), np.std(linear_vals, ddof=1)]

    bars = ax.bar(x, means, yerr=sds, capsize=10,
                  color=["#4C72B0", "#DD8452"], edgecolor="black", width=0.5)

    # Ghi giá trị lên mỗi bar
    for bar, val in zip(bars, means):
        ax.text(bar.get_x() + bar.get_width() / 2, val + 1,
                f"{val:.1f}", ha="center", fontsize=12, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([f"Graph\n(n={len(graph_vals)})",
                        f"Linear\n(n={len(linear_vals)})"])
    ax.set_ylabel("MAPQ median")
    ax.set_ylim(0, 70)
    ax.set_title("MAPQ median — cả hai đều đạt giá trị tối đa\n"
                 "→ Chất lượng mapping tương đương")
    ax.axhline(60, color="green", linestyle="--", alpha=0.5,
               label="MAPQ max = 60")
    ax.legend()
    ax.grid(alpha=0.3, axis="y")

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"📊 Đã lưu: {out_path}")


def run_descriptive_analysis(df):
    print("\n" + "=" * 70)
    print("PHẦN B — THỐNG KÊ MÔ TẢ (COVERAGE + MAPQ median)")
    print("=" * 70)

    os.makedirs(OUTPUT_DESCRIPTIVE_PLOT_DIR, exist_ok=True)

    rows = []
    for col, label in DESCRIPTIVE_COLS:
        if col not in df.columns:
            print(f"⚠️ Bỏ qua {col} — thiếu cột")
            continue
        rows.append(describe_series(df[col].values, label))

    # Thêm MAPQ median vào bảng mô tả (dù không vẽ histogram)
    for col, label in [("graph_mapq_median",  "Graph: MAPQ median"),
                       ("linear_mapq_median", "Linear: MAPQ median")]:
        if col in df.columns:
            rows.append(describe_series(df[col].values, label))

    if not rows:
        print("⚠️ Không có cột nào để mô tả.")
        return

    df_stats = pd.DataFrame(rows)

    print()
    print(df_stats.to_string(index=False))
    print()

    print("📌 Kiểm tra MAPQ median:")
    for col, label in [("graph_mapq_median",  "Graph: MAPQ median"),
                       ("linear_mapq_median", "Linear: MAPQ median")]:
        if col in df.columns:
            vals = df[col].dropna().values
            unique_vals = np.unique(vals)
            print(f"   {label}: giá trị duy nhất = {unique_vals.tolist()} "
                  f"→ {'✅ Tất cả = 60' if np.all(vals == 60) else '⚠️ Có giá trị khác 60'}")
    print()

    df_stats.to_csv(OUTPUT_DESCRIPTIVE_CSV, index=False)
    print(f"🎉 Đã lưu: {OUTPUT_DESCRIPTIVE_CSV}\n")

    # --- Vẽ histogram cho các chỉ số coverage ---
    for col, label in DESCRIPTIVE_COLS:
        if col not in df.columns:
            continue
        safe_name = col.replace(" ", "_")
        plot_path = os.path.join(OUTPUT_DESCRIPTIVE_PLOT_DIR, f"{safe_name}.png")
        plot_distribution(df[col].values, label, plot_path)

    # --- Vẽ paired plot MÔ TẢ cho depth ---
    print("\n📊 Vẽ paired plot MÔ TẢ cho depth (không kiểm định)...")
    for graph_col, linear_col, label in DESCRIPTIVE_PAIRED_COLS:
        if graph_col not in df.columns or linear_col not in df.columns:
            print(f"⚠️ Bỏ qua {label} — thiếu cột")
            continue

        sub = df[[graph_col, linear_col]].dropna()
        graph = sub[graph_col].values
        linear = sub[linear_col].values

        safe_name = f"paired_{graph_col}_vs_{linear_col}"
        plot_path = os.path.join(OUTPUT_DESCRIPTIVE_PLOT_DIR, f"{safe_name}.png")
        plot_descriptive_paired(graph, linear, label, plot_path)

    # --- Vẽ bar chart MAPQ median ---
    if "graph_mapq_median" in df.columns and "linear_mapq_median" in df.columns:
        print("\n📊 Vẽ bar chart so sánh MAPQ median...")
        mapq_plot_path = os.path.join(OUTPUT_DESCRIPTIVE_PLOT_DIR,
                                      "mapq_median_comparison.png")
        plot_mapq_median_bar(df, mapq_plot_path)


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("PHÂN TÍCH TỔNG HỢP: GRAPH vs LINEAR MAPPING")
    print("=" * 70)

    df = pd.read_csv(INPUT_CSV)
    print(f"\n📊 Số mẫu: {len(df)}")
    print(f"📋 Cột: {list(df.columns)}\n")

    run_comparison_analysis(df)
    run_descriptive_analysis(df)

    print("\n" + "=" * 70)
    print("🎉 HOÀN THÀNH TẤT CẢ")
    print("=" * 70)
    print(f"\n📂 Output:")
    print(f"   - {OUTPUT_STATS_CSV}")
    print(f"   - {OUTPUT_DESCRIPTIVE_CSV}")
    print(f"   - {OUTPUT_PLOT_DIR}/")
    print(f"   - {OUTPUT_DESCRIPTIVE_PLOT_DIR}/")


if __name__ == "__main__":
    main()
