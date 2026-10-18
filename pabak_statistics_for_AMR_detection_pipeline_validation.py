#!/usr/bin/env python3
"""
Bước 3: Tính PABAK (Prevalence-Adjusted Bias-Adjusted Kappa)
Input: summary_statistics.csv (từ bước 1)
Output: 
  - tables/pabak_results.csv
  - figures/pabak_formula.png
  - summary/pabak_summary.txt
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from datetime import datetime

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
# 2. TÍNH PABAK CHO TỪNG MẪU
# ============================================================================

print("\n📊 Đang tính PABAK...")

pabak_data = []

for _, row in df.iterrows():
    sample = row['sample']
    location = row['location']
    a = row['common']           # Common (cả 2 phát hiện)
    b = row['pipeline_only']    # Pipeline only
    c = row['rgi_only']         # RGI only
    
    # Tổng số gene union
    union = a + b + c
    
    # Po (Observed agreement)
    if union > 0:
        po = a / union
    else:
        po = 0
    
    # PABAK = 2 * Po - 1
    pabak = 2 * po - 1
    
    # Đánh giá mức độ đồng thuận (Landis & Koch)
    if pabak >= 0.81:
        level = 'Almost perfect'
    elif pabak >= 0.61:
        level = 'Substantial'
    elif pabak >= 0.41:
        level = 'Moderate'
    elif pabak >= 0.21:
        level = 'Fair'
    else:
        level = 'Slight'
    
    pabak_data.append({
        'sample': sample,
        'location': location,
        'common (a)': a,
        'pipeline_only (b)': b,
        'rgi_only (c)': c,
        'union (a+b+c)': union,
        'Po (a/union)': round(po, 4),
        'PABAK (2*Po-1)': round(pabak, 4),
        'agreement_level': level
    })

df_pabak = pd.DataFrame(pabak_data)

# Lưu file CSV
df_pabak.to_csv(TABLES_DIR / "pabak_results.csv", index=False)
print(f"✓ Đã lưu: tables/pabak_results.csv")

# ============================================================================
# 3. TÍNH THỐNG KÊ TỔNG HỢP
# ============================================================================

print("\n📊 Đang tính thống kê tổng hợp...")

mean_pabak = df_pabak['PABAK (2*Po-1)'].mean()
std_pabak = df_pabak['PABAK (2*Po-1)'].std()
min_pabak = df_pabak['PABAK (2*Po-1)'].min()
max_pabak = df_pabak['PABAK (2*Po-1)'].max()
median_pabak = df_pabak['PABAK (2*Po-1)'].median()

mean_po = df_pabak['Po (a/union)'].mean()
std_po = df_pabak['Po (a/union)'].std()

# Đánh giá mức độ đồng thuận tổng thể
if mean_pabak >= 0.81:
    overall_level = 'Almost perfect'
elif mean_pabak >= 0.61:
    overall_level = 'Substantial'
elif mean_pabak >= 0.41:
    overall_level = 'Moderate'
elif mean_pabak >= 0.21:
    overall_level = 'Fair'
else:
    overall_level = 'Slight'

# ============================================================================
# 4. TẠO SUMMARY REPORT
# ============================================================================

print("\n📊 Tạo Summary Report...")

summary_lines = []
summary_lines.append("="*70)
summary_lines.append("PABAK SUMMARY REPORT")
summary_lines.append("="*70)
summary_lines.append(f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
summary_lines.append(f"Number of samples: {len(df_pabak)}")
summary_lines.append("")

# ====== 1. Overall Statistics ======
summary_lines.append("-"*70)
summary_lines.append("1. OVERALL PABAK STATISTICS")
summary_lines.append("-"*70)
summary_lines.append("")
summary_lines.append(f"  Mean PABAK:   {mean_pabak:.4f}")
summary_lines.append(f"  Median PABAK: {median_pabak:.4f}")
summary_lines.append(f"  Std PABAK:    {std_pabak:.4f}")
summary_lines.append(f"  Min PABAK:    {min_pabak:.4f}")
summary_lines.append(f"  Max PABAK:    {max_pabak:.4f}")
summary_lines.append("")
summary_lines.append(f"  Mean Po:      {mean_po:.4f}")
summary_lines.append(f"  Std Po:       {std_po:.4f}")
summary_lines.append("")
summary_lines.append(f"  Overall agreement level: {overall_level}")

# ====== 2. Detailed Table ======
summary_lines.append("")
summary_lines.append("-"*70)
summary_lines.append("2. PABAK BY SAMPLE")
summary_lines.append("-"*70)
summary_lines.append("")
summary_lines.append(f"{'Sample':<15} {'Location':<12} {'a':<6} {'b':<6} {'c':<6} {'Po':<10} {'PABAK':<10} {'Agreement'}")
summary_lines.append("-"*70)

for _, row in df_pabak.iterrows():
    summary_lines.append(f"{row['sample']:<15} {row['location']:<12} {row['common (a)']:<6} {row['pipeline_only (b)']:<6} {row['rgi_only (c)']:<6} {row['Po (a/union)']:<10.4f} {row['PABAK (2*Po-1)']:<10.4f} {row['agreement_level']}")

# ====== 3. Interpretation ======
summary_lines.append("")
summary_lines.append("-"*70)
summary_lines.append("3. INTERPRETATION")
summary_lines.append("-"*70)
summary_lines.append("")
summary_lines.append("  PABAK (Prevalence-Adjusted Bias-Adjusted Kappa) was used")
summary_lines.append("  instead of Cohen's Kappa due to the high prevalence of AMR genes")
summary_lines.append("  in our dataset (prevalence paradox).")
summary_lines.append("")
summary_lines.append(f"  Mean PABAK = {mean_pabak:.4f} -> {overall_level} agreement")
summary_lines.append("")
summary_lines.append("  Interpretation of PABAK values (Landis & Koch):")
summary_lines.append("    - 0.81-1.00: Almost perfect")
summary_lines.append("    - 0.61-0.80: Substantial")
summary_lines.append("    - 0.41-0.60: Moderate")
summary_lines.append("    - 0.21-0.40: Fair")
summary_lines.append("    - 0.00-0.20: Slight")
summary_lines.append("    - < 0.00: Poor")

summary_lines.append("")
summary_lines.append("="*70)
summary_lines.append("END OF PABAK SUMMARY REPORT")
summary_lines.append("="*70)

# Lưu summary
summary_content = "\n".join(summary_lines)
with open(SUMMARY_DIR / "pabak_summary.txt", 'w') as f:
    f.write(summary_content)

print(f"✓ Đã lưu: summary/pabak_summary.txt")

# ============================================================================
# 5. VẼ ẢNH CÔNG THỨC PABAK
# ============================================================================

print("\n📊 Vẽ ảnh công thức PABAK...")

fig1, ax1 = plt.subplots(figsize=(10, 8))
ax1.axis('off')

text_content = """


                         PABAK FORMULA


    ┌─────────────────────────────────────────────────────────────┐
    │                                                             │
    │    PABAK = 2 × Po - 1                                      │
    │                                                             │
    └─────────────────────────────────────────────────────────────┘


    Trong đó:

    • Po = a / (a + b + c)
    
    • a = Common: gene được phát hiện bởi CẢ HAI phương pháp
    • b = Pipeline only: gene CHỈ pipeline phát hiện
    • c = RGI only: gene CHỈ RGI phát hiện


    ─────────────────────────────────────────────────────────────


    SO SÁNH VỚI COHEN'S KAPPA:

    • Cohen's Kappa bị ảnh hưởng bởi prevalence paradox
      (khi tỷ lệ gene xuất hiện quá cao hoặc quá thấp)

    • PABAK được điều chỉnh để khắc phục vấn đề này,
      phù hợp với dữ liệu có prevalence cao


    ─────────────────────────────────────────────────────────────


    ĐÁNH GIÁ MỨC ĐỘ ĐỒNG THUẬN (Landis & Koch):

    • 0.81 - 1.00: Almost perfect
    • 0.61 - 0.80: Substantial
    • 0.41 - 0.60: Moderate
    • 0.21 - 0.40: Fair
    • 0.00 - 0.20: Slight
    • < 0.00:      Poor

"""

ax1.text(0.5, 0.5, text_content, 
        transform=ax1.transAxes,
        fontsize=13, 
        verticalalignment='center',
        horizontalalignment='center',
        fontfamily='monospace',
        linespacing=1.8)

plt.tight_layout()
fig1.savefig(FIGURES_DIR / "pabak_formula.png", dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Đã lưu: pabak_formula.png")

# ============================================================================
# 6. TẠO BẢNG KẾT QUẢ DẠNG ẢNH (KHÔNG CẦN BIỂU ĐỒ)
# ============================================================================

print("\n📊 Tạo bảng kết quả PABAK dạng ảnh...")

fig2, ax2 = plt.subplots(figsize=(12, max(6, len(df_pabak) * 0.4 + 2)))
ax2.axis('tight')
ax2.axis('off')

# Chuẩn bị dữ liệu cho bảng
table_data = []
for _, row in df_pabak.iterrows():
    table_data.append([
        row['sample'],
        row['location'],
        row['common (a)'],
        row['pipeline_only (b)'],
        row['rgi_only (c)'],
        f"{row['Po (a/union)']:.4f}",
        f"{row['PABAK (2*Po-1)']:.4f}",
        row['agreement_level']
    ])

# Thêm hàng Mean
table_data.append([
    'MEAN',
    '',
    '',
    '',
    '',
    f"{mean_po:.4f}",
    f"{mean_pabak:.4f}",
    overall_level
])

columns = ['Sample', 'Location', 'a (Common)', 'b (Pipeline only)', 'c (RGI only)', 'Po', 'PABAK', 'Agreement Level']

# Tạo bảng
table = ax2.table(cellText=table_data, colLabels=columns,
                  cellLoc='center', loc='center',
                  colWidths=[0.15, 0.12, 0.08, 0.10, 0.08, 0.08, 0.10, 0.15])
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1, 1.5)

# Tô màu header
for i, col in enumerate(columns):
    table[(0, i)].set_facecolor('#2E86AB')
    table[(0, i)].set_text_props(weight='bold', color='white')

# Tô màu hàng Mean (hàng cuối)
for j in range(len(columns)):
    table[(len(table_data), j)].set_facecolor('#F18F01')
    table[(len(table_data), j)].set_text_props(weight='bold')

# Tô màu xen kẽ các hàng
for i in range(1, len(table_data)):
    for j in range(len(columns)):
        if i % 2 == 1:
            table[(i, j)].set_facecolor('#f5f5f5')

plt.tight_layout()
fig2.savefig(FIGURES_DIR / "pabak_results_table.png", dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Đã lưu: pabak_results_table.png")

# ============================================================================
# 7. IN KẾT QUẢ RA MÀN HÌNH
# ============================================================================

print("\n" + "="*60)
print("✅ HOÀN THÀNH BƯỚC 3: PABAK")
print("="*60)

print("\n📁 Các file đã tạo:")
print(f"  tables/pabak_results.csv          - Bảng dữ liệu PABAK")
print(f"  figures/pabak_formula.png         - Ảnh công thức PABAK")
print(f"  figures/pabak_results_table.png   - Bảng kết quả PABAK")
print(f"  summary/pabak_summary.txt         - Summary report")

print("\n📊 KẾT QUẢ CHÍNH:")
print(f"  - Mean PABAK: {mean_pabak:.4f} ({overall_level})")
print(f"  - Mean Po: {mean_po:.4f}")
print(f"  - PABAK range: {min_pabak:.4f} - {max_pabak:.4f}")
print(f"  - Số mẫu: {len(df_pabak)}")

print("\n📝 GHI CHÚ:")
print(f"  - PABAK = {mean_pabak:.4f} -> {overall_level} agreement between methods")
print(f"  - Summary report: summary/pabak_summary.txt")