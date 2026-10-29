#!/usr/bin/env python3
"""
Wilcoxon signed-rank test + biểu đồ phân bố Δ + thống kê tách Δ dương/âm.
Bản v6.1: giữ nguyên toàn bộ phân tích của v5, thêm thống kê tách
          Δ > 0 / Δ < 0 cho cả hai nhánh, sửa trục Y histogram chia
          1 đơn vị mẫu cho cả hai panel.

Input : bootstrap_both.tsv
Output:
  - wilcoxon_test.txt                 (báo cáo văn bản)
  - wilcoxon_test_summary.csv         (bảng tóm tắt)
  - wilcoxon_test.png                 (histogram Δ có dấu + CI)
  - wilcoxon_distribution_check.png   (boxplot 1x2 — giữ y nguyên)
  - wilcoxon_split_stats.txt          (thống kê tách Δ dương/âm)
  - wilcoxon_split_stats.csv          (bảng tách dạng CSV)
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import skew, probplot
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

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
OUTPUT_SPLIT_TXT = os.path.join(OUTPUT_DIR, "wilcoxon_split_stats.txt")
OUTPUT_SPLIT_CSV = os.path.join(OUTPUT_DIR, "wilcoxon_split_stats.csv")

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
# PHÂN TÍCH CHÍNH
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
# LƯU TEXT CHÍNH
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
# LƯU CSV TÓM TẮT
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
print(f"✅ Đã lưu CSV tóm tắt: {OUTPUT_CSV}")

# ============================================================
# THỐNG KÊ MÔ TẢ TÁCH Δ DƯƠNG / Δ ÂM
# ============================================================
def describe_group(values, group_label):
    """Tính thống kê mô tả đầy đủ cho một nhóm Δ."""
    a = np.asarray(values, dtype=int)
    if len(a) == 0:
        return None
    return {
        'group': group_label,
        'n': len(a),
        'min': int(a.min()),
        'max': int(a.max()),
        'mean': float(a.mean()),
        'median': float(np.median(a)),
        'q1': float(np.percentile(a, 25)),
        'q3': float(np.percentile(a, 75)),
        'sum': int(a.sum()),
    }


split_results = {}
print("\n" + "=" * 78)
print("THỐNG KÊ TÁCH Δ DƯƠNG / Δ ÂM")
print("=" * 78)

for br in BRANCHES:
    r = results[br['key']]
    deltas = r['deltas']

    pos_vals = deltas[deltas > 0]
    neg_vals = deltas[deltas < 0]

    d_pos = describe_group(pos_vals, 'improved')
    d_neg = describe_group(neg_vals, 'worsened')
    split_results[br['key']] = {
        'label': r['label'],
        'pos': d_pos,
        'neg': d_neg,
    }

    print(f"\n{r['label']}:")

    if d_pos:
        print(f"  Δ > 0 (cải thiện): n = {d_pos['n']}")
        print(f"    min = {d_pos['min']:+d}, max = {d_pos['max']:+d}")
        print(f"    mean = {d_pos['mean']:+.2f}, "
              f"median = {d_pos['median']:+.2f}")
        print(f"    Q1 = {d_pos['q1']:+.1f}, Q3 = {d_pos['q3']:+.1f}")
        print(f"    Tổng Δ = {d_pos['sum']:+d}")

    if d_neg:
        print(f"  Δ < 0 (giảm): n = {d_neg['n']}")
        print(f"    min = {d_neg['min']:+d}, max = {d_neg['max']:+d}")
        print(f"    mean = {d_neg['mean']:+.2f}, "
              f"median = {d_neg['median']:+.2f}")
        print(f"    Q1 = {d_neg['q1']:+.1f}, Q3 = {d_neg['q3']:+.1f}")
        print(f"    Tổng Δ = {d_neg['sum']:+d}")

    if d_pos and d_neg:
        print(f"  So sánh:")
        print(f"    Tỷ lệ n (cải thiện/giảm) = "
              f"{d_pos['n']}/{d_neg['n']} = "
              f"{d_pos['n']/d_neg['n']:.2f}")
        print(f"    Tỷ lệ |mean| = {abs(d_pos['mean'])/abs(d_neg['mean']):.2f}")
        print(f"    Tỷ lệ |median| = "
              f"{abs(d_pos['median'])/abs(d_neg['median']):.2f}")
        print(f"    Tỷ lệ tổng Δ (dương/âm) = "
              f"{d_pos['sum']}/{abs(d_neg['sum'])} = "
              f"{d_pos['sum']/abs(d_neg['sum']):.2f}")

# ============================================================
# LƯU SPLIT TXT
# ============================================================
with open(OUTPUT_SPLIT_TXT, "w", encoding="utf-8") as f:
    f.write("THỐNG KÊ MÔ TẢ Δ TÁCH THEO HƯỚNG THAY ĐỔI\n")
    f.write("=" * 78 + "\n\n")
    f.write("Mục đích: làm nổi bật hiệu quả thực sự của core genome bằng\n")
    f.write("cách tách riêng nhóm cải thiện (Δ > 0) và nhóm giảm (Δ < 0),\n")
    f.write("thay vì chỉ xem median/mean của toàn bộ nhóm (bị kéo xuống\n")
    f.write("bởi các mẫu giảm).\n\n")

    for br in BRANCHES:
        sr = split_results[br['key']]
        f.write("=" * 78 + "\n")
        f.write(f"{sr['label']}\n")
        f.write("=" * 78 + "\n\n")

        for d, title in [(sr['pos'], 'Δ > 0 (cải thiện)'),
                          (sr['neg'], 'Δ < 0 (giảm)')]:
            if d is None:
                continue
            f.write(f"  {title}:\n")
            f.write(f"    n = {d['n']}\n")
            f.write(f"    min = {d['min']:+d}, max = {d['max']:+d}\n")
            f.write(f"    mean = {d['mean']:+.2f}\n")
            f.write(f"    median = {d['median']:+.2f}\n")
            f.write(f"    Q1 = {d['q1']:+.1f}, Q3 = {d['q3']:+.1f}\n")
            f.write(f"    Tổng Δ = {d['sum']:+d}\n\n")

        if sr['pos'] and sr['neg']:
            dp, dn = sr['pos'], sr['neg']
            f.write(f"  So sánh:\n")
            f.write(f"    Tỷ lệ n = {dp['n']}/{dn['n']} = "
                    f"{dp['n']/dn['n']:.2f}\n")
            f.write(f"    Tỷ lệ |mean| = "
                    f"{abs(dp['mean'])/abs(dn['mean']):.2f}\n")
            f.write(f"    Tỷ lệ |median| = "
                    f"{abs(dp['median'])/abs(dn['median']):.2f}\n")
            f.write(f"    Tỷ lệ tổng Δ = "
                    f"{dp['sum']}/{abs(dn['sum'])} = "
                    f"{dp['sum']/abs(dn['sum']):.2f}\n\n")

print(f"\n✅ Đã lưu split stats text: {OUTPUT_SPLIT_TXT}")

# ============================================================
# LƯU SPLIT CSV
# ============================================================
split_rows = []
for br in BRANCHES:
    sr = split_results[br['key']]
    for d, group_label in [(sr['pos'], 'Δ > 0 (cải thiện)'),
                            (sr['neg'], 'Δ < 0 (giảm)')]:
        if d is None:
            continue
        split_rows.append({
            'branch': br['key'],
            'group': group_label,
            'n': d['n'],
            'min': d['min'],
            'max': d['max'],
            'mean': f"{d['mean']:.2f}",
            'median': f"{d['median']:.2f}",
            'q1': f"{d['q1']:.1f}",
            'q3': f"{d['q3']:.1f}",
            'sum': d['sum'],
        })
pd.DataFrame(split_rows).to_csv(OUTPUT_SPLIT_CSV, sep='\t', index=False)
print(f"✅ Đã lưu split stats CSV: {OUTPUT_SPLIT_CSV}")

# ============================================================
# BIỂU ĐỒ 1 — HISTOGRAM Δ CÓ DẤU + CI
# (cột Δ=0 cùng width, không rug plot, trục Y chia 1 đơn vị mẫu)
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))

for ax, br in zip(axes, BRANCHES):
    r = results[br['key']]
    deltas = r['deltas']

    neg_vals = deltas[deltas < 0]
    pos_vals = deltas[deltas > 0]
    zero_count = int((deltas == 0).sum())

    span = max(abs(r['min']), abs(r['max']))
    neg_bins = np.arange(-span - 5, -4.999, 10)
    pos_bins = np.arange(5, span + 15, 10)

    if len(neg_vals) > 0 and len(neg_bins) > 1:
        ax.hist(neg_vals, bins=neg_bins,
                color='#C44E52', alpha=0.8,
                edgecolor='black', linewidth=0.5,
                label='Δ < 0 (giảm)')
    if len(pos_vals) > 0 and len(pos_bins) > 1:
        ax.hist(pos_vals, bins=pos_bins,
                color='#55A868', alpha=0.8,
                edgecolor='black', linewidth=0.5,
                label='Δ > 0 (cải thiện)')
    if zero_count > 0:
        ax.bar([0], [zero_count], width=10,
               color='#BDBDBD', alpha=0.9,
               edgecolor='black', linewidth=0.8,
               label='Δ = 0 (không đổi)')

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

    # Chia vạch trục Y theo 1 đơn vị mẫu cho cả hai panel
    ax.yaxis.set_major_locator(
        MaxNLocator(nbins=20, integer=True, steps=[1]))

plt.suptitle('Wilcoxon Signed-Rank Test — Xu hướng cải thiện Bootstrap\n'
             '(nhóm test chính, đã loại mẫu cả hai đều =100)',
             fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(OUTPUT_PNG_MAIN, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"✅ Đã lưu biểu đồ chính: {OUTPUT_PNG_MAIN}")

# ============================================================
# BIỂU ĐỒ 2 — BOXPLOT (GIỮ Y NGUYÊN)
# ============================================================
fig_dist, axes_dist = plt.subplots(2, 2, figsize=(14, 10))

for i, br in enumerate(BRANCHES):
    r = results[br['key']]
    deltas = r['deltas_nz']
    sk = r['skewness']

    ax_qq = axes_dist[i, 0]
    probplot(deltas, dist="norm", plot=ax_qq)
    ax_qq.set_title(f"QQ plot — {r['label']}\n"
                    f"(skewness = {sk:+.2f})",
                    fontsize=11, fontweight='bold')
    ax_qq.set_xlabel('Theoretical quantiles (normal)', fontsize=10)
    ax_qq.set_ylabel('Sample quantiles (Δ)', fontsize=10)
    ax_qq.grid(alpha=0.3)

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
# TÓM TẮT CHÍNH
# ============================================================
print("\n" + "=" * 78)
print("TÓM TẮT CHÍNH")
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

# ============================================================
# TÓM TẮT SPLIT
# ============================================================
print("\n" + "=" * 78)
print("TÓM TẮT THỐNG KÊ TÁCH Δ")
print("=" * 78)
print(f"\n{'Nhánh':<28} {'Nhóm':<20} {'n':>4} {'median':>8} "
      f"{'mean':>8} {'Q1':>6} {'Q3':>6}")
print("-" * 78)
for br in BRANCHES:
    sr = split_results[br['key']]
    for d, gname in [(sr['pos'], 'Δ > 0 (cải thiện)'),
                      (sr['neg'], 'Δ < 0 (giảm)')]:
        if d is None:
            continue
        print(f"{br['label']:<28} {gname:<20} {d['n']:>4} "
              f"{d['median']:>+7.2f} {d['mean']:>+7.2f} "
              f"{d['q1']:>+6.1f} {d['q3']:>+6.1f}")

print()
print("=" * 60)
print("✅ HOÀN THÀNH")
print("=" * 60)
print(f"Files đã xuất trong {OUTPUT_DIR}:")
print("  - wilcoxon_test.txt")
print("  - wilcoxon_test_summary.csv")
print("  - wilcoxon_test.png")
print("  - wilcoxon_distribution_check.png")
print("  - wilcoxon_split_stats.txt")
print("  - wilcoxon_split_stats.csv")
