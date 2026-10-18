#!/usr/bin/env python3
"""
Bước 2: Tính Recovery Rate, vẽ Dot Plot, Công thức và Summary Report
Input: summary_statistics.csv (từ bước 1)
Output: 
  - tables/recovery_rate.csv
  - figures/recovery_rate_dotplot.png
  - figures/recovery_rate_formula.png
  - summary/recovery_rate_summary.txt
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
# 2. TÍNH RECOVERY RATE
# ============================================================================

print("\n📊 Đang tính Recovery Rate...")

recovery_data = []

for _, row in df.iterrows():
    sample = row['sample']
    location = row['location']
    common = row['common']
    rgi_only = row['rgi_only']
    pipeline_only = row['pipeline_only']
    
    total_rgi = common + rgi_only
    total_pipeline = common + pipeline_only
    union_count = common + rgi_only + pipeline_only
    
    # Pipeline Recovery Rate (phần trăm gene RGI mà pipeline phát hiện lại được)
    if total_rgi > 0:
        pipeline_recovery_rate = (common / total_rgi) * 100
    else:
        pipeline_recovery_rate = float('nan')
    
    # RGI Recovery Rate (phần trăm gene pipeline mà RGI phát hiện lại được)
    if total_pipeline > 0:
        rgi_recovery_rate = (common / total_pipeline) * 100
    else:
        rgi_recovery_rate = float('nan')
    
    # Common proportion
    if union_count > 0:
        common_proportion = (common / union_count) * 100
    else:
        common_proportion = 0
    
    # Missed rate (phần trăm gene RGI bị pipeline bỏ sót)
    if total_rgi > 0:
        missed_rgi_rate = (rgi_only / total_rgi) * 100
    else:
        missed_rgi_rate = float('nan')
    
    # Extra rate (phần trăm gene pipeline phát hiện thêm so với RGI)
    if total_rgi > 0:
        extra_pipeline_rate = (pipeline_only / total_rgi) * 100
    else:
        extra_pipeline_rate = float('nan')
    
    recovery_data.append({
        'sample': sample,
        'location': location,
        'total_rgi': total_rgi,
        'total_pipeline': total_pipeline,
        'common': common,
        'rgi_only': rgi_only,
        'pipeline_only': pipeline_only,
        'union_count': union_count,
        'pipeline_recovery_rate': round(pipeline_recovery_rate, 2),
        'rgi_recovery_rate': round(rgi_recovery_rate, 2),
        'common_proportion': round(common_proportion, 2),
        'missed_rgi_rate': round(missed_rgi_rate, 2) if not pd.isna(missed_rgi_rate) else 'N/A',
        'extra_pipeline_rate': round(extra_pipeline_rate, 2) if not pd.isna(extra_pipeline_rate) else 'N/A'
    })

df_recovery = pd.DataFrame(recovery_data)

# Lưu file CSV
df_recovery.to_csv(TABLES_DIR / "recovery_rate.csv", index=False)
print(f"✓ Đã lưu: tables/recovery_rate.csv")

# ============================================================================
# 3. TÍNH THỐNG KÊ TỔNG HỢP
# ============================================================================

print("\n📊 Đang tính thống kê tổng hợp...")

pipeline_mean = df_recovery['pipeline_recovery_rate'].mean()
pipeline_std = df_recovery['pipeline_recovery_rate'].std()
pipeline_min = df_recovery['pipeline_recovery_rate'].min()
pipeline_max = df_recovery['pipeline_recovery_rate'].max()
pipeline_median = df_recovery['pipeline_recovery_rate'].median()

common_mean = df_recovery['common_proportion'].mean()
common_std = df_recovery['common_proportion'].std()
common_min = df_recovery['common_proportion'].min()
common_max = df_recovery['common_proportion'].max()
common_median = df_recovery['common_proportion'].median()

location_stats = df_recovery.groupby('location').agg({
    'pipeline_recovery_rate': ['mean', 'std', 'min', 'max', 'count'],
    'common_proportion': ['mean', 'std', 'min', 'max']
}).round(2)

# ============================================================================
# 4. TẠO SUMMARY REPORT
# ============================================================================

print("\n📊 Tạo Summary Report...")

summary_lines = []
summary_lines.append("="*70)
summary_lines.append("RECOVERY RATE SUMMARY REPORT")
summary_lines.append("="*70)
summary_lines.append(f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
summary_lines.append(f"Number of samples: {len(df_recovery)}")
summary_lines.append("")

# ====== 1. Pipeline Recovery Rate (vs RGI) ======
summary_lines.append("-"*70)
summary_lines.append("1. PIPELINE RECOVERY RATE (vs RGI)")
summary_lines.append("-"*70)
summary_lines.append("")
summary_lines.append(f"  Mean:   {pipeline_mean:.2f}%")
summary_lines.append(f"  Median: {pipeline_median:.2f}%")
summary_lines.append(f"  Std:    {pipeline_std:.2f}%")
summary_lines.append(f"  Min:    {pipeline_min:.2f}%")
summary_lines.append(f"  Max:    {pipeline_max:.2f}%")
summary_lines.append("")
summary_lines.append(f"  Interpretation:")
summary_lines.append(f"  Pipeline recovers {pipeline_mean:.2f}% of RGI-detected genes on average")
summary_lines.append(f"  Pipeline misses {100 - pipeline_mean:.2f}% of RGI-detected genes on average")
summary_lines.append("")

# ====== 2. RGI Recovery Rate (vs Pipeline) ======
summary_lines.append("-"*70)
summary_lines.append("2. RGI RECOVERY RATE (vs Pipeline)")
summary_lines.append("-"*70)
summary_lines.append("")

rgi_mean = df_recovery['rgi_recovery_rate'].mean()
rgi_std = df_recovery['rgi_recovery_rate'].std()
rgi_min = df_recovery['rgi_recovery_rate'].min()
rgi_max = df_recovery['rgi_recovery_rate'].max()
rgi_median = df_recovery['rgi_recovery_rate'].median()

summary_lines.append(f"  Mean:   {rgi_mean:.2f}%")
summary_lines.append(f"  Median: {rgi_median:.2f}%")
summary_lines.append(f"  Std:    {rgi_std:.2f}%")
summary_lines.append(f"  Min:    {rgi_min:.2f}%")
summary_lines.append(f"  Max:    {rgi_max:.2f}%")
summary_lines.append("")
summary_lines.append(f"  Interpretation:")
summary_lines.append(f"  RGI recovers {rgi_mean:.2f}% of Pipeline-detected genes on average")
summary_lines.append(f"  RGI misses {100 - rgi_mean:.2f}% of Pipeline-detected genes on average")
summary_lines.append("")

# ====== 3. Common Proportion ======
summary_lines.append("-"*70)
summary_lines.append("3. COMMON PROPORTION (genes detected by both methods)")
summary_lines.append("-"*70)
summary_lines.append("")
summary_lines.append(f"  Mean:   {common_mean:.2f}%")
summary_lines.append(f"  Median: {common_median:.2f}%")
summary_lines.append(f"  Std:    {common_std:.2f}%")
summary_lines.append(f"  Min:    {common_min:.2f}%")
summary_lines.append(f"  Max:    {common_max:.2f}%")
summary_lines.append("")
summary_lines.append(f"  Interpretation:")
if common_mean > 90:
    summary_lines.append(f"  ✓ High agreement between methods ({common_mean:.1f}% of genes detected by both)")
elif common_mean > 70:
    summary_lines.append(f"  ✓ Moderate agreement between methods ({common_mean:.1f}% of genes detected by both)")
else:
    summary_lines.append(f"  ✗ Low agreement between methods ({common_mean:.1f}% of genes detected by both)")
summary_lines.append("")

# ====== 4. By Location ======
summary_lines.append("-"*70)
summary_lines.append("4. PIPELINE RECOVERY RATE BY LOCATION")
summary_lines.append("-"*70)
summary_lines.append("")
summary_lines.append(f"{'Location':<15} {'n':<5} {'Mean (%)':<12} {'Std':<10} {'Min':<10} {'Max':<10}")
summary_lines.append("-"*70)

for loc in sorted(location_stats.index):
    row = location_stats.loc[loc]
    n = int(row[('pipeline_recovery_rate', 'count')])
    mean = row[('pipeline_recovery_rate', 'mean')]
    std = row[('pipeline_recovery_rate', 'std')]
    min_val = row[('pipeline_recovery_rate', 'min')]
    max_val = row[('pipeline_recovery_rate', 'max')]
    summary_lines.append(f"{loc:<15} {n:<5} {mean:<12.2f} {std:<10.2f} {min_val:<10.2f} {max_val:<10.2f}")

summary_lines.append("")

# ====== 5. Summary Table ======
summary_lines.append("-"*70)
summary_lines.append("5. DETAILED RECOVERY RATE BY SAMPLE")
summary_lines.append("-"*70)
summary_lines.append("")
summary_lines.append(f"{'Sample':<15} {'Location':<12} {'Pipeline Rec %':<15} {'RGI Rec %':<15} {'Common %':<10}")
summary_lines.append("-"*70)

for _, row in df_recovery.iterrows():
    summary_lines.append(f"{row['sample']:<15} {row['location']:<12} {row['pipeline_recovery_rate']:<15.2f} {row['rgi_recovery_rate']:<15.2f} {row['common_proportion']:<10.2f}")

summary_lines.append("")

# ====== 6. Conclusion ======
summary_lines.append("="*70)
summary_lines.append("6. CONCLUSION")
summary_lines.append("="*70)
summary_lines.append("")

if pipeline_mean > 90 and common_mean > 90:
    summary_lines.append("  ✓ The pangenome pipeline shows HIGH sensitivity")
    summary_lines.append(f"  ✓ Pipeline recovers {pipeline_mean:.1f}% of RGI-detected genes")
    summary_lines.append(f"  ✓ High agreement between methods ({common_mean:.1f}% common genes)")
    summary_lines.append("  ✓ The pipeline is suitable for AMR gene detection")
elif pipeline_mean > 80:
    summary_lines.append("  ✓ The pangenome pipeline shows GOOD sensitivity")
    summary_lines.append(f"  ✓ Pipeline recovers {pipeline_mean:.1f}% of RGI-detected genes")
    summary_lines.append(f"  ✓ Agreement between methods: {common_mean:.1f}% common genes")
else:
    summary_lines.append("  ✗ The pangenome pipeline shows LOW sensitivity")
    summary_lines.append(f"  ✗ Pipeline recovers only {pipeline_mean:.1f}% of RGI-detected genes")
    summary_lines.append(f"  ✗ Agreement between methods: {common_mean:.1f}% common genes")

summary_lines.append("")
summary_lines.append("="*70)
summary_lines.append("END OF RECOVERY RATE SUMMARY REPORT")
summary_lines.append("="*70)

# Lưu summary
summary_content = "\n".join(summary_lines)
with open(SUMMARY_DIR / "recovery_rate_summary.txt", 'w') as f:
    f.write(summary_content)

print(f"✓ Đã lưu: summary/recovery_rate_summary.txt")

# ============================================================================
# 5. VẼ BIỂU ĐỒ DOT PLOT - KHÔNG HIỂN THỊ GIÁ TRỊ %
# ============================================================================

print("\n📊 Vẽ biểu đồ Dot Plot...")

# Sắp xếp theo location
location_order = {'Japan': 1, 'Taiwan': 2, 'Singapore': 3, 'Myanmar': 4, 
                  'Vietnam': 5, 'Thailand': 6, 'Indonesia': 7}
df_recovery['location_order'] = df_recovery['location'].map(location_order).fillna(99)
df_plot = df_recovery.sort_values(['location_order', 'sample'])

fig1, ax1 = plt.subplots(figsize=(12, 7))

# Tạo màu cho từng location
locations = sorted(df_plot['location'].unique())
colors = {
    'Japan': '#3498DB', 
    'Taiwan': '#2ECC71', 
    'Singapore': '#F1C40F', 
    'Myanmar': '#E67E22', 
    'Vietnam': '#E74C3C', 
    'Thailand': '#9B59B6', 
    'Indonesia': '#1ABC9C'
}

# Vẽ dot plot - KHÔNG HIỂN THỊ GIÁ TRỊ %
for i, loc in enumerate(locations):
    loc_data = df_plot[df_plot['location'] == loc]
    x_pos = [i + 0.15 * j for j in range(len(loc_data))]
    
    ax1.scatter(x_pos, loc_data['pipeline_recovery_rate'], 
                s=250, color=colors.get(loc, 'gray'), 
                edgecolors='black', linewidth=0.8, alpha=0.85)

# Đường baseline 100%
ax1.axhline(100, color='red', linestyle='--', alpha=0.5, linewidth=1.5)

# Đường mean
mean_val = df_plot['pipeline_recovery_rate'].mean()
ax1.axhline(mean_val, color='blue', linestyle='-', alpha=0.5, linewidth=1.5)

ax1.set_xlabel('Geographic Location', fontsize=13, fontweight='bold')
ax1.set_ylabel('Pipeline Recovery Rate vs RGI (%)', fontsize=13, fontweight='bold')
ax1.set_title('Pipeline Recovery Rate by Geographic Location', fontsize=15, fontweight='bold')
ax1.set_xticks(range(len(locations)))
ax1.set_xticklabels(locations, fontsize=11)
ax1.grid(axis='y', alpha=0.3, linestyle='--')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.set_ylim(90, 105)

plt.tight_layout()
fig1.savefig(FIGURES_DIR / "recovery_rate_dotplot.png", dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Đã lưu: recovery_rate_dotplot.png (không hiển thị giá trị %)")

# ============================================================================
# 6. VẼ ẢNH CÔNG THỨC RECOVERY RATE
# ============================================================================

print("\n📊 Vẽ ảnh công thức Recovery Rate...")

fig2, ax2 = plt.subplots(figsize=(10, 8))
ax2.axis('off')

text_content = """


                    RECOVERY RATE FORMULA


    ┌─────────────────────────────────────────────────────────────┐
    │                                                             │
    │                     common                                  │
    │    Recovery Rate = ──────────────────  × 100%              │
    │                   common + rgi_only                         │
    │                                                             │
    └─────────────────────────────────────────────────────────────┘


    Trong đó:

    • common = Số gene được phát hiện bởi CẢ HAI phương pháp
    • rgi_only = Số gene CHỈ RGI phát hiện (pipeline không phát hiện)


    ─────────────────────────────────────────────────────────────


    GIẢI THÍCH:

    • Recovery Rate là tỷ lệ phần trăm gene của RGI (CARD)
      mà pipeline có thể phát hiện lại được.

    • Recovery Rate = 100%: Pipeline phát hiện được 100%
      gene mà RGI phát hiện.

    • Recovery Rate < 100%: Pipeline bỏ sót một số gene
      mà RGI phát hiện.

    • Recovery Rate càng cao: Pipeline càng tốt trong việc
      phát hiện gene của RGI.


    ─────────────────────────────────────────────────────────────


    VÍ DỤ:

    Nếu common = 31 và rgi_only = 1:
    Recovery Rate = 31/(31+1) × 100% = 96.88%

    → Pipeline phát hiện lại được 96.88% gene của RGI

"""

ax2.text(0.5, 0.5, text_content, 
        transform=ax2.transAxes,
        fontsize=13, 
        verticalalignment='center',
        horizontalalignment='center',
        fontfamily='monospace',
        linespacing=1.8)

plt.tight_layout()
fig2.savefig(FIGURES_DIR / "recovery_rate_formula.png", dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Đã lưu: recovery_rate_formula.png")

# ============================================================================
# 7. IN KẾT QUẢ RA MÀN HÌNH
# ============================================================================

print("\n" + "="*60)
print("✅ HOÀN THÀNH BƯỚC 2: RECOVERY RATE")
print("="*60)

print("\n📁 Các file đã tạo:")
print(f"  tables/recovery_rate.csv              - Bảng dữ liệu recovery rate")
print(f"  figures/recovery_rate_dotplot.png      - Biểu đồ dot plot (không hiển thị %)")
print(f"  figures/recovery_rate_formula.png      - Ảnh công thức")
print(f"  summary/recovery_rate_summary.txt      - Summary report")

print("\n📊 KẾT QUẢ CHÍNH:")
print(f"  - Pipeline Recovery Rate (mean): {pipeline_mean:.2f}%")
print(f"  - Pipeline Recovery Rate (range): {pipeline_min:.2f}% - {pipeline_max:.2f}%")
print(f"  - Common Proportion (mean): {common_mean:.2f}%")
print(f"  - Số mẫu: {len(df_recovery)}")

print("\n📝 GHI CHÚ:")
print(f"  - Pipeline recovers {pipeline_mean:.2f}% of RGI-detected genes on average")
print(f"  - {common_mean:.1f}% of genes are detected by both methods")
print(f"  - Summary report: summary/recovery_rate_summary.txt")
