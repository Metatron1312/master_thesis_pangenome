#!/usr/bin/env python3
"""
Bước 4: Tính McNemar's Test
Input: summary_statistics.csv (từ bước 1)
Output: 
  - tables/mcnemar_results.csv
  - figures/mcnemar_formula.png
  - figures/mcnemar_results_table.png
  - summary/mcnemar_summary.txt
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from datetime import datetime
from scipy.stats import chi2

# ============================================================================
# 0. CẤU HÌNH ĐƯỜNG DẪN
# ============================================================================

BASE_PATH = Path("/Users/shinra/gene_detection_v2/validate_database")
INPUT_FILE = BASE_PATH / "validation_charts_and_statistic/tables/summary_statistics.csv"
OUTPUT_DIR = BASE_PATH / "validation_charts_and_statistic"
TABLES_DIR = OUTPUT_DIR / "tables"
FIGURES_DIR = OUTPUT_DIR / "figures"
SUMMARY_DIR = OUTPUT_DIR / "summary"

# Tạo thư mục nếu chưa có
TABLES_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
SUMMARY_DIR.mkdir(parents=True, exist_ok=True)

print(f"✓ Input file: {INPUT_FILE}")
print(f"✓ Output tables: {TABLES_DIR}")
print(f"✓ Output figures: {FIGURES_DIR}")
print(f"✓ Output summary: {SUMMARY_DIR}")
print("="*60)

# ============================================================================
# 1. ĐỌC DỮ LIỆU
# ============================================================================

print("\n📊 Đang đọc dữ liệu từ file summary_statistics.csv...")
df = pd.read_csv(INPUT_FILE)
print(f"✓ Đọc thành công {len(df)} mẫu")

# ============================================================================
# 2. TÍNH MCNEMAR'S TEST CHO TỪNG MẪU
# ============================================================================

print("\n📊 Đang tính McNemar's Test...")

mcnemar_data = []

for _, row in df.iterrows():
    sample = row['sample']
    location = row['location']
    a = row['common']           # Common (cả 2 phát hiện)
    b = row['pipeline_only']    # Pipeline only
    c = row['rgi_only']         # RGI only
    d = 0                       # d = 0 (union set)
    
    # Tổng số discordant pairs
    discordant = b + c
    
    if discordant > 0:
        # McNemar's Chi-square
        chi2_val = ((b - c) ** 2) / discordant
        
        # P-value từ chi-square distribution (df=1)
        p_value = chi2.sf(chi2_val, df=1)
        
        # Kiểm tra ý nghĩa thống kê
        if p_value < 0.05:
            significant = 'Yes'
            if b > c:
                direction = 'Pipeline > RGI'
            elif c > b:
                direction = 'RGI > Pipeline'
            else:
                direction = 'Equal'
        else:
            significant = 'No'
            direction = 'Not significant'
    else:
        chi2_val = float('nan')
        p_value = float('nan')
        significant = 'N/A'
        direction = 'No discordant pairs'
    
    mcnemar_data.append({
        'sample': sample,
        'location': location,
        'a (common)': a,
        'b (pipeline_only)': b,
        'c (rgi_only)': c,
        'discordant (b+c)': discordant,
        'chi_square': round(chi2_val, 4) if not pd.isna(chi2_val) else 'N/A',
        'p_value': round(p_value, 6) if not pd.isna(p_value) else 'N/A',
        'significant_p<0.05': significant,
        'direction': direction
    })

df_mcnemar = pd.DataFrame(mcnemar_data)

# Lưu file CSV
df_mcnemar.to_csv(TABLES_DIR / "mcnemar_results.csv", index=False)
print(f"✓ Đã lưu: tables/mcnemar_results.csv")

# ============================================================================
# 3. TÍNH MCNEMAR'S TEST GỘP (POOLED ANALYSIS)
# ============================================================================

print("\n📊 Đang tính McNemar's Test gộp...")

# Tổng hợp từ tất cả mẫu
total_b = df['pipeline_only'].sum()
total_c = df['rgi_only'].sum()
total_a = df['common'].sum()
total_discordant = total_b + total_c

if total_discordant > 0:
    pool_chi2 = ((total_b - total_c) ** 2) / total_discordant
    pool_p_value = chi2.sf(pool_chi2, df=1)
    
    if pool_p_value < 0.05:
        pool_significant = 'Yes'
        if total_b > total_c:
            pool_direction = 'Pipeline detects significantly more genes than RGI'
        elif total_c > total_b:
            pool_direction = 'RGI detects significantly more genes than Pipeline'
        else:
            pool_direction = 'Equal'
    else:
        pool_significant = 'No'
        pool_direction = 'No significant difference'
else:
    pool_chi2 = float('nan')
    pool_p_value = float('nan')
    pool_significant = 'N/A'
    pool_direction = 'No discordant pairs'

print(f"\n  Pooled McNemar's Test:")
print(f"    Total a (common): {total_a}")
print(f"    Total b (pipeline only): {total_b}")
print(f"    Total c (RGI only): {total_c}")
print(f"    Chi-square: {pool_chi2:.4f}")
print(f"    P-value: {pool_p_value:.6f}")
print(f"    Significant: {pool_significant}")
print(f"    Direction: {pool_direction}")

# ============================================================================
# 4. TẠO SUMMARY REPORT
# ============================================================================

print("\n📊 Tạo Summary Report...")

summary_lines = []
summary_lines.append("="*70)
summary_lines.append("MCNEMAR'S TEST SUMMARY REPORT")
summary_lines.append("="*70)
summary_lines.append(f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
summary_lines.append(f"Number of samples: {len(df_mcnemar)}")
summary_lines.append("")

# ====== 1. Pooled Analysis ======
summary_lines.append("-"*70)
summary_lines.append("1. POOLED MCNEMAR'S TEST (All samples combined)")
summary_lines.append("-"*70)
summary_lines.append("")
summary_lines.append(f"  Total a (common): {total_a}")
summary_lines.append(f"  Total b (pipeline only): {total_b}")
summary_lines.append(f"  Total c (RGI only): {total_c}")
summary_lines.append(f"  Discordant pairs (b+c): {total_discordant}")
summary_lines.append("")
summary_lines.append(f"  Chi-square: {pool_chi2:.4f}")
summary_lines.append(f"  P-value: {pool_p_value:.6f}")
summary_lines.append("")

if pool_p_value < 0.05:
    summary_lines.append(f"  ✓ Significant difference (p < 0.05)")
    summary_lines.append(f"  ✓ {pool_direction}")
else:
    summary_lines.append(f"  ✗ No significant difference (p >= 0.05)")
    summary_lines.append(f"  ✗ {pool_direction}")

# ====== 2. Individual Sample Analysis ======
summary_lines.append("")
summary_lines.append("-"*70)
summary_lines.append("2. INDIVIDUAL SAMPLE ANALYSIS")
summary_lines.append("-"*70)
summary_lines.append("")
summary_lines.append(f"{'Sample':<15} {'Location':<12} {'b':<8} {'c':<8} {'χ²':<10} {'p-value':<12} {'Significant'}")
summary_lines.append("-"*70)

for _, row in df_mcnemar.iterrows():
    summary_lines.append(f"{row['sample']:<15} {row['location']:<12} {row['b (pipeline_only)']:<8} {row['c (rgi_only)']:<8} {row['chi_square']:<10} {row['p_value']:<12} {row['significant_p<0.05']}")

# ====== 3. Interpretation ======
summary_lines.append("")
summary_lines.append("-"*70)
summary_lines.append("3. INTERPRETATION")
summary_lines.append("-"*70)
summary_lines.append("")

if pool_p_value < 0.05:
    summary_lines.append("  ✓ McNemar's Test shows a SIGNIFICANT difference between")
    summary_lines.append(f"    the two methods (p = {pool_p_value:.6f})")
    summary_lines.append("")
    if total_b > total_c:
        summary_lines.append("  ✓ Pipeline Pangenome detects significantly MORE AMR genes")
        summary_lines.append(f"    than RGI (CARD): {total_b} vs {total_c} additional genes")
    else:
        summary_lines.append("  ✓ RGI (CARD) detects significantly MORE AMR genes")
        summary_lines.append(f"    than Pipeline Pangenome: {total_c} vs {total_b} additional genes")
else:
    summary_lines.append("  ✗ McNemar's Test shows NO significant difference between")
    summary_lines.append(f"    the two methods (p = {pool_p_value:.6f})")
    summary_lines.append("  ✗ The two methods have comparable performance")

summary_lines.append("")
summary_lines.append("="*70)
summary_lines.append("END OF MCNEMAR'S TEST SUMMARY REPORT")
summary_lines.append("="*70)

# Lưu summary
summary_content = "\n".join(summary_lines)
with open(SUMMARY_DIR / "mcnemar_summary.txt", 'w') as f:
    f.write(summary_content)

print(f"✓ Đã lưu: summary/mcnemar_summary.txt")

# ============================================================================
# 5. VẼ ẢNH CÔNG THỨC MCNEMAR
# ============================================================================

print("\n📊 Vẽ ảnh công thức McNemar's Test...")

fig1, ax1 = plt.subplots(figsize=(10, 8))
ax1.axis('off')

text_content = """


                    MCNEMAR'S TEST


    ┌─────────────────────────────────────────────────────────────┐
    │                                                             │
    │         (b - c)²                                           │
    │    χ² = ──────────                                         │
    │          (b + c)                                           │
    │                                                             │
    └─────────────────────────────────────────────────────────────┘

    df = 1


    ─────────────────────────────────────────────────────────────


    BẢNG 2x2:

                    RGI (CARD)
                Present   Absent
    Pipeline    Present   a         b
    Pangenome   Absent    c         d

    Trong đó:

    • a = Common (cả 2 phát hiện)
    • b = Pipeline only (CHỈ pipeline phát hiện)
    • c = RGI only (CHỈ RGI phát hiện)
    • d = Cả 2 đều không phát hiện


    ─────────────────────────────────────────────────────────────


    GIẢI THÍCH:

    • McNemar's Test chỉ quan tâm đến sự khác biệt giữa
      hai phương pháp (ô b và ô c)

    • b > c: Pipeline phát hiện nhiều hơn RGI
    • b < c: RGI phát hiện nhiều hơn Pipeline
    • b = c: Hai phương pháp tương đương


    ─────────────────────────────────────────────────────────────


    QUYẾT ĐỊNH:

    • p < 0.05: Sự khác biệt có ý nghĩa thống kê
    • p ≥ 0.05: Không có sự khác biệt có ý nghĩa

"""

ax1.text(0.5, 0.5, text_content, 
        transform=ax1.transAxes,
        fontsize=13, 
        verticalalignment='center',
        horizontalalignment='center',
        fontfamily='monospace',
        linespacing=1.8)

plt.tight_layout()
fig1.savefig(FIGURES_DIR / "mcnemar_formula.png", dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Đã lưu: mcnemar_formula.png")

# ============================================================================
# 6. TẠO BẢNG KẾT QUẢ DẠNG ẢNH
# ============================================================================

print("\n📊 Tạo bảng kết quả McNemar dạng ảnh...")

fig2, ax2 = plt.subplots(figsize=(14, max(6, len(df_mcnemar) * 0.4 + 2)))
ax2.axis('tight')
ax2.axis('off')

# Chuẩn bị dữ liệu cho bảng
table_data = []
for _, row in df_mcnemar.iterrows():
    table_data.append([
        row['sample'],
        row['location'],
        row['b (pipeline_only)'],
        row['c (rgi_only)'],
        row['discordant (b+c)'],
        str(row['chi_square']),
        str(row['p_value']),
        row['significant_p<0.05']
    ])

# Thêm hàng Pooled
table_data.append([
    'POOLED',
    '',
    total_b,
    total_c,
    total_discordant,
    f"{pool_chi2:.4f}",
    f"{pool_p_value:.6f}",
    'Yes' if pool_p_value < 0.05 else 'No'
])

columns = ['Sample', 'Location', 'b (Pipeline only)', 'c (RGI only)', 'b+c', 'χ²', 'p-value', 'Significant (p<0.05)']

# Tạo bảng
table = ax2.table(cellText=table_data, colLabels=columns,
                  cellLoc='center', loc='center',
                  colWidths=[0.12, 0.10, 0.10, 0.08, 0.08, 0.10, 0.12, 0.12])
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1, 1.5)

# Tô màu header
for i, col in enumerate(columns):
    table[(0, i)].set_facecolor('#2E86AB')
    table[(0, i)].set_text_props(weight='bold', color='white')

# Tô màu hàng Pooled (hàng cuối)
for j in range(len(columns)):
    table[(len(table_data), j)].set_facecolor('#F18F01')
    table[(len(table_data), j)].set_text_props(weight='bold')

# Tô màu xen kẽ các hàng
for i in range(1, len(table_data)):
    for j in range(len(columns)):
        if i % 2 == 1:
            table[(i, j)].set_facecolor('#f5f5f5')

plt.tight_layout()
fig2.savefig(FIGURES_DIR / "mcnemar_results_table.png", dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Đã lưu: mcnemar_results_table.png")

# ============================================================================
# 7. IN KẾT QUẢ RA MÀN HÌNH
# ============================================================================

print("\n" + "="*60)
print("✅ HOÀN THÀNH BƯỚC 4: MCNEMAR'S TEST")
print("="*60)

print("\n📁 Các file đã tạo:")
print(f"  tables/mcnemar_results.csv          - Bảng dữ liệu McNemar")
print(f"  figures/mcnemar_formula.png         - Ảnh công thức McNemar")
print(f"  figures/mcnemar_results_table.png   - Bảng kết quả McNemar")
print(f"  summary/mcnemar_summary.txt         - Summary report")

print("\n📊 KẾT QUẢ CHÍNH (POOLED):")
print(f"  - Total b (pipeline only): {total_b}")
print(f"  - Total c (RGI only): {total_c}")
print(f"  - Chi-square: {pool_chi2:.4f}")
print(f"  - P-value: {pool_p_value:.6f}")
print(f"  - Significant: {pool_significant}")
print(f"  - Direction: {pool_direction}")

print("\n📝 GHI CHÚ:")
print(f"  - McNemar's Test: p = {pool_p_value:.6f} {'(SIGNIFICANT)' if pool_p_value < 0.05 else '(NOT significant)'}")
print(f"  - Summary report: summary/mcnemar_summary.txt")
