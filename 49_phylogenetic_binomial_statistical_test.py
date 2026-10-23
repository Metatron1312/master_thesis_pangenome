import os
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np

# ============================================================
# CẤU HÌNH
# ============================================================

INPUT_TSV = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/comparison/bootstrap_both.tsv"
OUTPUT_DIR = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/comparison"
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_TXT = os.path.join(OUTPUT_DIR, "binomial_test.txt")
OUTPUT_PNG = os.path.join(OUTPUT_DIR, "binomial_test.png")

# ============================================================
# ĐỌC DỮ LIỆU
# ============================================================

df = pd.read_csv(INPUT_TSV, sep='\t')
REFERENCE_NAMES = ['BJAB07104', 'NC_021726.1']
df = df[~df['sample'].isin(REFERENCE_NAMES)]

print("=" * 70)
print("BINOMIAL TEST — TỶ LỆ CẢI THIỆN BOOTSTRAP")
print("=" * 70)
print(f"Số mẫu: {len(df)}")

# ============================================================
# TÍNH TOÁN
# ============================================================

results = {}

for col, label in [('diff_closest', 'Smallest internal branch'),
                    ('diff_parent', 'Parent internal branch')]:
    n_total = len(df)
    n_improved = (df[col] > 0).sum()
    n_not_improved = n_total - n_improved
    percent = n_improved / n_total * 100
    p_value = stats.binomtest(n_improved, n_total, 0.5).pvalue
    
    results[col] = {
        'label': label,
        'n_total': n_total,
        'n_improved': n_improved,
        'n_not_improved': n_not_improved,
        'percent': percent,
        'p_value': p_value,
    }
    
    print(f"\n{label}:")
    print(f"  Tổng số mẫu: {n_total}")
    print(f"  Số mẫu cải thiện (Δ > 0): {n_improved} ({percent:.1f}%)")
    print(f"  Số mẫu không cải thiện (Δ ≤ 0): {n_not_improved} ({100-percent:.1f}%)")
    print(f"  Binomial test p-value: {p_value:.4f}")
    
    if p_value < 0.05:
        print(f"  → Tỷ lệ cải thiện KHÁC 50% có ý nghĩa thống kê (p < 0.05)")
    else:
        print(f"  → Tỷ lệ cải thiện KHÔNG khác 50% (p >= 0.05)")

# ============================================================
# LƯU KẾT QUẢ
# ============================================================

with open(OUTPUT_TXT, "w") as f:
    f.write("BINOMIAL TEST — TỶ LỆ CẢI THIỆN BOOTSTRAP\n")
    f.write("=" * 70 + "\n\n")
    f.write(f"Số mẫu: {len(df)}\n\n")
    
    for col, r in results.items():
        f.write(f"{r['label']}:\n")
        f.write(f"  Tổng số mẫu: {r['n_total']}\n")
        f.write(f"  Số mẫu cải thiện: {r['n_improved']} ({r['percent']:.1f}%)\n")
        f.write(f"  Số mẫu không cải thiện: {r['n_not_improved']} ({100-r['percent']:.1f}%)\n")
        f.write(f"  Binomial test p-value: {r['p_value']:.4f}\n\n")

print(f"\nĐã lưu kết quả: {OUTPUT_TXT}")

# ============================================================
# VẼ BIỂU ĐỒ (LEGEND HÀNG TRÊN)
# ============================================================

fig, axes = plt.subplots(1, 2, figsize=(14, 8))

# Biểu đồ 1: Tỷ lệ cải thiện (stacked bar)
ax1 = axes[0]

categories = ['Smallest\ninternal branch', 'Parent\ninternal branch']
improved_pct = [results['diff_closest']['percent'],
                results['diff_parent']['percent']]
not_improved_pct = [100 - p for p in improved_pct]

x = np.arange(len(categories))
width = 0.5

bars1 = ax1.bar(x, improved_pct, width, label='Improved (Δ > 0)', color='#55A868')
bars2 = ax1.bar(x, not_improved_pct, width, bottom=improved_pct,
                label='Not improved (Δ ≤ 0)', color='#C44E52')

# Thêm giá trị %
for i, (imp, not_imp) in enumerate(zip(improved_pct, not_improved_pct)):
    ax1.text(i, imp/2, f'{imp:.1f}%', ha='center', va='center',
             fontsize=12, color='white', fontweight='bold')
    ax1.text(i, imp + not_imp/2, f'{not_imp:.1f}%', ha='center', va='center',
             fontsize=12, color='white', fontweight='bold')

# Đường 50% reference
ax1.axhline(y=50, color='black', linestyle='--', linewidth=2, label='50% reference')

ax1.set_ylabel('Tỷ lệ (%)', fontsize=12)
ax1.set_title('Tỷ lệ mẫu cải thiện Bootstrap', fontsize=13, fontweight='bold')
ax1.set_xticks(x)
ax1.set_xticklabels(categories, fontsize=11)
# LEGEND HÀNG TRÊN, GIỮA
ax1.legend(fontsize=10, loc='upper center', bbox_to_anchor=(0.5, 1.15),
           ncol=3, frameon=True)
ax1.set_ylim(0, 110)
ax1.grid(axis='y', alpha=0.3)

# Thêm p-value
for i, col in enumerate(['diff_closest', 'diff_parent']):
    p = results[col]['p_value']
    color = '#55A868' if p < 0.05 else '#C44E52'
    ax1.text(i, 105, f'p = {p:.4f}', ha='center', va='bottom',
             fontsize=11, fontweight='bold', color=color)

# Biểu đồ 2: Số lượng mẫu
ax2 = axes[1]

n_improved = [results['diff_closest']['n_improved'],
              results['diff_parent']['n_improved']]
n_not_improved = [results['diff_closest']['n_not_improved'],
                  results['diff_parent']['n_not_improved']]

x = np.arange(len(categories))
width = 0.35

bars1 = ax2.bar(x - width/2, n_improved, width, label='Improved (Δ > 0)', color='#55A868')
bars2 = ax2.bar(x + width/2, n_not_improved, width, label='Not improved (Δ ≤ 0)', color='#C44E52')

for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        ax2.annotate(f'{height:.0f}',
                    xy=(bar.get_x() + bar.get_width()/2, height),
                    xytext=(0, 3), textcoords="offset points",
                    ha='center', va='bottom', fontsize=11)

ax2.set_ylabel('Số mẫu', fontsize=12)
ax2.set_title('Số lượng mẫu cải thiện / không cải thiện', fontsize=13, fontweight='bold')
ax2.set_xticks(x)
ax2.set_xticklabels(categories, fontsize=11)
# LEGEND HÀNG TRÊN, GIỮA
ax2.legend(fontsize=10, loc='upper center', bbox_to_anchor=(0.5, 1.15),
           ncol=2, frameon=True)
ax2.set_ylim(0, max(n_not_improved) * 1.2)
ax2.grid(axis='y', alpha=0.3)

plt.suptitle('Binomial Test — Tỷ lệ cải thiện Bootstrap\n(Core Genome vs Core Gene)',
             fontsize=14, fontweight='bold', y=1.05)
plt.tight_layout()
plt.savefig(OUTPUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

print(f"Đã lưu biểu đồ: {OUTPUT_PNG}")

# ============================================================
# TÓM TẮT
# ============================================================

print("\n" + "=" * 70)
print("TÓM TẮT")
print("=" * 70)

print(f"\n{'Loại nhánh':<30} {'Số cải thiện':>15} {'Tỷ lệ':>10} {'p-value':>10}")
print("-" * 70)

for col, r in results.items():
    print(f"{r['label']:<30} {r['n_improved']:>15}/{r['n_total']} {r['percent']:>9.1f}% {r['p_value']:>10.4f}")

print("\nKết luận:")
for col, r in results.items():
    if r['p_value'] < 0.05:
        print(f"  ✅ {r['label']}: Tỷ lệ cải thiện có ý nghĩa thống kê (p < 0.05)")
    else:
        print(f"  ❌ {r['label']}: Tỷ lệ cải thiện không có ý nghĩa (p >= 0.05)")
