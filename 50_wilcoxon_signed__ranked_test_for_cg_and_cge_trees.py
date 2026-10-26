#!/usr/bin/env python3
"""
Wilcoxon signed-rank test + biểu đồ phân bố Δ.
Bản v5: cột Δ=0 trong histogram có cùng width với các cột khác
        (không bị đè, không nhỏ hơn), bỏ rug plot, boxplot giữ nguyên.

Input : bootstrap_both.tsv
Output:
  - wilcoxon_test.txt              (báo cáo văn bản)
  - wilcoxon_test_summary.csv      (bảng tóm tắt)
  - wilcoxon_test.png              (histogram Δ có dấu + CI)
  - wilcoxon_distribution_check.png (boxplot 1x2 — giữ y nguyên)
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import skew, probplot
import matplotlib.pyplot as plt

# ============================================================
# CẤU HÌNH
# ============================================================
INPUT_TSV = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/comparison/bootstrap_both.tsv"
OUTPUT_DIR = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/comparison"
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_TXT = os.path.join(OUTPUT_DIR, "wilcoxon_test.txt")
OUTPUT_CSV = os.path.join(OUTPUT_DIR, "wilcoxon_test_summary.csv")
OUTPUT_PNG_MAIN = os.path.join(OUTPUT_DIR, "wilcoxon_test.png")
OUTPUT_PNG_DIST = os.path.join(OUTPUT_DIR, "wilcoxon_distribution_check.png")

REFERENCE_NAMES = ['BJAB07104', 'NC_021726.1']
N_BOOTSTRAP = 10000
RANDOM_SEED = 42

# ============================================================
# ĐỌC DỮ LIỆU
# ============================================================
df = pd.read_csv(INPUT_TSV, sep='\t')
df = df[~df['sample'].isin(REFERENCE_NAMES)]

print("=" * 78)
print("WILCOXON SIGNED-RANK TEST + KIỂM TRA PHÂN PHỐI")
print("=" * 78)
print(f"Số mẫu sau khi loại genome tham chiếu: {len(df)}")
print()
print("Tiêu chí chọn mẫu:")
print("  Loại mẫu mà CẢ HAI cây đều đạt bootstrap = 100")
print("  (bs_cg = 100 VÀ bs_cge = 100). Mọi mẫu còn lại — bao gồm")
print("  cả mẫu core genome GIẢM từ 100 xuống <100 — đều được giữ.")
print()

# ============================================================
# ĐỊNH NGHĨA NHÁNH
# ============================================================
BRANCHES = [
    {
        'key': 'closest',
        'label': 'Smallest internal branch',
        'bs_cg':  'bs_closest_cg',
        'bs_cge': 'bs_closest_cge',
        'diff':   'diff_closest',
    },
    {
        'key': 'parent',
        'label': 'Parent internal branch',
        'bs_cg':  'bs_parent_cg',
        'bs_cge': 'bs_parent_cge',
        'diff':   'diff_parent',
    },
]

# ============================================================
# PHÂN TÍCH
# ============================================================
results = {}

for br in BRANCHES:
    mask_keep = ~((df[br['bs_cg']] == 100) & (df[br['bs_cge']] == 100))
    df_keep = df[mask_keep]
    deltas = df_keep[br['diff']].to_numpy(dtype=int)

    deltas_nz = deltas[deltas != 0]
    n_total_keep = len(deltas)
    n_used = len(deltas_nz)
    n_dropped_zero = n_total_keep - n_used

    stat_two, p_two = stats.wilcoxon(deltas_nz, alternative='two-sided')
    stat_greater, p_greater = stats.wilcoxon(deltas_nz,
                                              alternative='greater')

    n_pos = int((deltas_nz > 0).sum())
    n_neg = int((deltas_nz < 0).sum())
    median_all = float(np.median(deltas))
    mean_all = float(deltas.mean())
    min_val = int(deltas.min())
    max_val = int(deltas.max())
    skew_val = float(skew(deltas_nz))

    rng = np.random.default_rng(RANDOM_SEED)
    boot_medians = np.array([
        np.median(rng.choice(deltas, size=len(deltas), replace=True))
        for _ in range(N_BOOTSTRAP)
    ])
    ci_lo, ci_hi = np.percentile(boot_medians, [2.5, 97.5])

    results[br['key']] = {
        'label': br['label'],
        'n_total_keep': n_total_keep,
        'n_used': n_used,
        'n_dropped_zero': n_dropped_zero,
        'n_pos': n_pos,
        'n_neg': n_neg,
        'median': median_all,
        'mean': mean_all,
        'min': min_val,
        'max': max_val,
        'skewness': skew_val,
        'median_ci_lo': float(ci_lo),
        'median_ci_hi': float(ci_hi),
        'W_two': float(stat_two),
        'p_two': float(p_two),
        'W_greater': float(stat_greater),
        'p_greater': float(p_greater),
        'deltas': deltas,
        'deltas_nz': deltas_nz,
    }

    print(f"\n{br['label']}:")
    print(f"  Nhóm test chính: {n_total_keep} mẫu "
          f"(loại {n_dropped_zero} mẫu Δ=0 khỏi Wilcoxon)")
    print(f"  Số mẫu dùng trong test: {n_used}")
    print(f"    Δ > 0: {n_pos}")
    print(f"    Δ < 0: {n_neg}")
    print(f"  Phân bố Δ:")
    print(f"    min = {min_val:+d}, max = {max_val:+d}")
    print(f"    mean = {mean_all:+.2f}, median = {median_all:+.2f}")
    print(f"    skewness = {skew_val:+.2f}")
    if abs(skew_val) < 0.5:
        print(f"    → Gần đối xứng: Wilcoxon rất đáng tin.")
    elif abs(skew_val) < 1.0:
        print(f"    → Lệch vừa: Wilcoxon ổn, nên báo cáo kèm CI.")
    else:
        print(f"    → Lệch mạnh: cần diễn giải thận trọng, ưu tiên CI.")
    print(f"  Bootstrap 95% CI cho median: "
          f"[{ci_lo:+.2f}, {ci_hi:+.2f}]")
    print(f"  Wilcoxon W = {stat_two:.1f}")
    print(f"    p-value (two-sided): {p_two:.4g}")
    print(f"    p-value (greater)  : {p_greater:.4g}")

    if p_two < 0.05:
        direction = "CẢI THIỆN" if median_all > 0 else "GIẢM"
        print(f"  → Median Δ khác 0 có ý nghĩa (p < 0.05), "
              f"hướng {direction}.")
    else:
        print(f"  → Median Δ KHÔNG khác 0 có ý nghĩa (p >= 0.05).")
    if ci_lo > 0 or ci_hi < 0:
        print(f"  → 95% CI của median KHÔNG chứa 0 — xu hướng có ý nghĩa "
              f"theo bootstrap.")
    else:
        print(f"  → 95% CI của median CHỨA 0 — chưa đủ bằng chứng về "
              f"xu hướng.")

# ============================================================
# LƯU TEXT
# ============================================================
with open(OUTPUT_TXT, "w", encoding="utf-8") as f:
    f.write("WILCOXON SIGNED-RANK TEST + KIỂM TRA PHÂN PHỐI\n")
    f.write("=" * 78 + "\n\n")
    f.write(f"Số mẫu sau khi loại genome tham chiếu: {len(df)}\n\n")
    f.write("Tiêu chí chọn mẫu:\n")
    f.write("  Loại mẫu mà CẢ HAI cây đều đạt bootstrap = 100.\n")
    f.write("  Mọi mẫu còn lại — bao gồm cả mẫu core genome GIẢM từ\n")
    f.write("  100 xuống <100 — đều được giữ.\n\n")
    f.write("Wilcoxon xét cả DẤU và ĐỘ LỚN của Δ (qua xếp hạng |Δ|).\n")
    f.write("Kèm bootstrap CI cho median — kiểm định không phụ thuộc\n")
    f.write("giả định đối xứng.\n\n")

    for br in BRANCHES:
        r = results[br['key']]
        f.write("=" * 78 + "\n")
        f.write(f"{r['label']}\n")
        f.write("=" * 78 + "\n\n")
        f.write(f"  Nhóm test chính: {r['n_total_keep']} mẫu "
                f"(loại {r['n_dropped_zero']} mẫu Δ=0 khỏi Wilcoxon)\n")
        f.write(f"  Số mẫu dùng trong test: {r['n_used']}\n")
        f.write(f"    Δ > 0: {r['n_pos']}\n")
        f.write(f"    Δ < 0: {r['n_neg']}\n\n")
        f.write(f"  Phân bố Δ:\n")
        f.write(f"    min = {r['min']:+d}, max = {r['max']:+d}\n")
        f.write(f"    mean = {r['mean']:+.2f}, "
                f"median = {r['median']:+.2f}\n")
        f.write(f"    skewness = {r['skewness']:+.2f}\n")
        if abs(r['skewness']) < 0.5:
            f.write(f"    → Gần đối xứng: Wilcoxon rất đáng tin.\n")
        elif abs(r['skewness']) < 1.0:
            f.write(f"    → Lệch vừa: Wilcoxon ổn, nên báo cáo kèm CI.\n")
        else:
            f.write(f"    → Lệch mạnh: cần diễn giải thận trọng.\n")
        f.write(f"\n")
        f.write(f"  Bootstrap 95% CI cho median: "
                f"[{r['median_ci_lo']:+.2f}, "
                f"{r['median_ci_hi']:+.2f}]\n")
        f.write(f"  Wilcoxon W = {r['W_two']:.1f}\n")
        f.write(f"    p-value (two-sided): {r['p_two']:.4g}\n")
        f.write(f"    p-value (greater)  : {r['p_greater']:.4g}\n\n")
        if r['p_two'] < 0.05:
            direction = "cải thiện" if r['median'] > 0 else "giảm"
            f.write(f"  → Median Δ khác 0 có ý nghĩa (p < 0.05), "
                    f"hướng {direction}.\n\n")
        else:
            f.write(f"  → Median Δ KHÔNG khác 0 có ý nghĩa (p >= 0.05).\n\n")

print(f"\n✅ Đã lưu text: {OUTPUT_TXT}")

# ============================================================
# LƯU CSV
# ============================================================
rows = []
for br in BRANCHES:
    r = results[br['key']]
    rows.append({
        'branch': br['key'],
        'n_total_keep': r['n_total_keep'],
        'n_used': r['n_used'],
        'n_dropped_zero': r['n_dropped_zero'],
        'n_pos': r['n_pos'],
        'n_neg': r['n_neg'],
        'min': r['min'],
        'max': r['max'],
        'median': f"{r['median']:.2f}",
        'mean': f"{r['mean']:.2f}",
        'skewness': f"{r['skewness']:.2f}",
        'median_ci_lo': f"{r['median_ci_lo']:.2f}",
        'median_ci_hi': f"{r['median_ci_hi']:.2f}",
        'W': f"{r['W_two']:.1f}",
        'p_two': f"{r['p_two']:.4g}",
        'p_greater': f"{r['p_greater']:.4g}",
    })
pd.DataFrame(rows).to_csv(OUTPUT_CSV, sep='\t', index=False)
print(f"✅ Đã lưu CSV: {OUTPUT_CSV}")

# ============================================================
# BIỂU ĐỒ 1 — HISTOGRAM Δ CÓ DẤU + CI
# (cột Δ=0 cùng width với các cột khác, không bị đè, không rug plot)
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))

for ax, br in zip(axes, BRANCHES):
    r = results[br['key']]
    deltas = r['deltas']

    neg_vals = deltas[deltas < 0]
    pos_vals = deltas[deltas > 0]
    zero_count = int((deltas == 0).sum())

    span = max(abs(r['min']), abs(r['max']))
    # Bin rộng 10: vùng âm kết thúc ở -5, vùng dương bắt đầu từ +5
    neg_bins = np.arange(-span - 5, -4.999, 10)   # ... -25, -15, -5
    pos_bins = np.arange(5, span + 15, 10)         # 5, 15, 25 ...

    # Histogram Δ < 0
    if len(neg_vals) > 0 and len(neg_bins) > 1:
        ax.hist(neg_vals, bins=neg_bins,
                color='#C44E52', alpha=0.8,
                edgecolor='black', linewidth=0.5,
                label='Δ < 0 (giảm)')
    # Histogram Δ > 0
    if len(pos_vals) > 0 and len(pos_bins) > 1:
        ax.hist(pos_vals, bins=pos_bins,
                color='#55A868', alpha=0.8,
                edgecolor='black', linewidth=0.5,
                label='Δ > 0 (cải thiện)')
    # Cột Δ = 0: width = 10 (cùng width với histogram), từ -5 đến +5
    if zero_count > 0:
        ax.bar([0], [zero_count], width=10,
               color='#BDBDBD', alpha=0.9,
               edgecolor='black', linewidth=0.8,
               label='Δ = 0 (không đổi)')

    # Đường tham chiếu
    ax.axvline(0, color='black', linestyle='--', linewidth=1)
    ax.axvline(r['median'], color='blue', linestyle=':',
               linewidth=2,
               label=f"median = {r['median']:+.1f} "
                     f"(95% CI [{r['median_ci_lo']:+.0f}, "
                     f"{r['median_ci_hi']:+.0f}])")

    p_label = f"Wilcoxon p = {r['p_two']:.4g}"
    color_p = '#2E7D32' if r['p_two'] < 0.05 else '#C62828'
    ax.text(0.02, 0.97, p_label, transform=ax.transAxes,
            fontsize=11, fontweight='bold', color=color_p, va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    ax.set_xlabel('Δ bootstrap (core genome − core gene)', fontsize=11)
    ax.set_ylabel('Số mẫu', fontsize=11)
    ax.set_title(r['label'], fontsize=12, fontweight='bold')
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(axis='y', alpha=0.3)

plt.suptitle('Wilcoxon Signed-Rank Test — Xu hướng cải thiện Bootstrap\n'
             '(nhóm test chính, đã loại mẫu cả hai đều =100)',
             fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(OUTPUT_PNG_MAIN, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"✅ Đã lưu biểu đồ chính: {OUTPUT_PNG_MAIN}")

# ============================================================
# BIỂU ĐỒ 2 — BOXPLOT (GIỮ Y NGUYÊN NHƯ SCRIPT GỐC)
# ============================================================
fig_dist, axes_dist = plt.subplots(2, 2, figsize=(14, 10))

for i, br in enumerate(BRANCHES):
    r = results[br['key']]
    deltas = r['deltas_nz']
    sk = r['skewness']

    # QQ plot
    ax_qq = axes_dist[i, 0]
    probplot(deltas, dist="norm", plot=ax_qq)
    ax_qq.set_title(f"QQ plot — {r['label']}\n"
                    f"(skewness = {sk:+.2f})",
                    fontsize=11, fontweight='bold')
    ax_qq.set_xlabel('Theoretical quantiles (normal)', fontsize=10)
    ax_qq.set_ylabel('Sample quantiles (Δ)', fontsize=10)
    ax_qq.grid(alpha=0.3)

    # Boxplot — GIỮ Y NGUYÊN
    ax_bx = axes_dist[i, 1]
    box = ax_bx.boxplot([deltas], vert=True, patch_artist=True,
                        labels=[r['label']])
    for patch in box['boxes']:
        patch.set_facecolor('#7FB3D5')
        patch.set_alpha(0.7)
    ax_bx.axhline(0, color='red', linestyle='--', linewidth=1.2,
                  label='Δ = 0')
    ax_bx.axhline(r['median'], color='blue', linestyle=':',
                  linewidth=1.5,
                  label=f"median = {r['median']:+.1f}")
    ax_bx.axhline(r['mean'], color='orange', linestyle=':',
                  linewidth=1.5,
                  label=f"mean = {r['mean']:+.1f}")
    ax_bx.set_ylabel('Δ bootstrap', fontsize=10)
    ax_bx.set_title(f"Boxplot — {r['label']}", fontsize=11,
                    fontweight='bold')
    ax_bx.legend(fontsize=8)
    ax_bx.grid(alpha=0.3)

plt.suptitle('Kiểm tra phân phối Δ (đánh giá giả định đối xứng '
             'của Wilcoxon)',
             fontsize=13, fontweight='bold', y=1.00)
plt.tight_layout()
plt.savefig(OUTPUT_PNG_DIST, dpi=300, bbox_inches='tight',
            facecolor='white')
plt.close()
print(f"✅ Đã lưu biểu đồ kiểm tra phân phối: {OUTPUT_PNG_DIST}")

# ============================================================
# TÓM TẮT
# ============================================================
print("\n" + "=" * 78)
print("TÓM TẮT")
print("=" * 78)
print(f"\n{'Loại nhánh':<28} {'n':>4} {'median':>8} {'mean':>8} "
      f"{'skew':>6} {'CI 95% median':>18} {'p (2-side)':>12}")
print("-" * 78)
for br in BRANCHES:
    r = results[br['key']]
    ci_str = f"[{r['median_ci_lo']:+.1f}, {r['median_ci_hi']:+.1f}]"
    print(f"{r['label']:<28} {r['n_used']:>4} "
          f"{r['median']:>+7.2f} {r['mean']:>+7.2f} "
          f"{r['skewness']:>+5.2f} {ci_str:>18} "
          f"{r['p_two']:>12.4g}")

print("\nKết luận:")
for br in BRANCHES:
    r = results[br['key']]
    if r['p_two'] < 0.05:
        direction = "cải thiện" if r['median'] > 0 else "giảm"
        ci_check = ("CI không chứa 0" if (r['median_ci_lo'] > 0 or
                                          r['median_ci_hi'] < 0)
                    else "CI chứa 0")
        print(f"  ✅ {r['label']}: median Δ khác 0 có ý nghĩa "
              f"(p = {r['p_two']:.4g}, {ci_check}), hướng {direction}")
    else:
        print(f"  ❌ {r['label']}: median Δ không khác 0 "
              f"(p = {r['p_two']:.4g})")

print()
