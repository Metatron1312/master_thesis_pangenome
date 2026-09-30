#!/usr/bin/env python3
"""
Script hoàn chỉnh:
1. Scatter plot (PNG + PDF) — không annotate, in giá trị 2 mẫu ra terminal
2. Bảng Word (.docx) — dạng merge cells, chuẩn tạp chí
3. Bảng ảnh (PNG + PDF) — font Arial cho PowerPoint
"""

import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import wilcoxon, shapiro

from docx import Document
from docx.shared import Pt, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_ALIGN_VERTICAL, WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

# ============================================================
# CẤU HÌNH
# ============================================================
CSV = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/Representative_analysis/base_mapping_rate_v2/base_mapping_rate_identity70.csv"
OUTDIR = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/Representative_analysis/base_mapping_rate_v2/plots"
os.makedirs(OUTDIR, exist_ok=True)

mpl.rcParams.update({
    "font.family"      : "Arial",
    "font.sans-serif"  : ["Arial"],
    "font.size"        : 24,
    "axes.titlesize"   : 24,
    "axes.labelsize"   : 24,
    "xtick.labelsize"  : 20,
    "ytick.labelsize"  : 20,
    "legend.fontsize"  : 20,
    "axes.grid"        : True,
    "grid.alpha"       : 0.25,
    "axes.spines.top"  : False,
    "axes.spines.right": False,
    "figure.dpi"       : 100,
    "savefig.dpi"      : 300,
    "savefig.bbox"     : "tight",
})

COLORS = {
    "linear"    : "#E74C3C",
    "pangenome" : "#2E86C1",
    "win"       : "#27AE60",
    "lose"      : "#E67E22",
    "tie"       : "#95A5A6",
}

# ============================================================
# ĐỌC DỮ LIỆU
# ============================================================
print("Đọc CSV...")
df = pd.read_csv(CSV)
pivot = df.pivot(index="biosample_id", columns="ref_type",
                 values="base_mapping_rate(%)").reset_index()
pivot.columns.name = None
pivot = pivot.rename(columns={"linear": "lin", "pangenome": "pan"})
pivot["delta"] = pivot["pan"] - pivot["lin"]

n_samples = len(pivot)
n_win  = (pivot["delta"] > 0.01).sum()
n_lose = (pivot["delta"] < -0.01).sum()
n_tie  = ((pivot["delta"].abs() <= 0.01)).sum()

# Identity stats
ident_lin_mean = df[df["ref_type"] == "linear"]["mean_identity(%)"].mean()
ident_pan_mean = df[df["ref_type"] == "pangenome"]["mean_identity(%)"].mean()
ident_lin_med  = df[df["ref_type"] == "linear"]["median_identity(%)"].mean()
ident_pan_med  = df[df["ref_type"] == "pangenome"]["median_identity(%)"].mean()

# Kiểm định
w_stat, w_p_two = wilcoxon(pivot["pan"], pivot["lin"], alternative="two-sided")
shapiro_stat, shapiro_p = shapiro(pivot["delta"])

def fmt_p(p):
    if p < 1e-10:
        return "< 0.0000000001"
    elif p < 0.001:
        return f"{p:.10f}".rstrip('0')
    else:
        return f"{p:.6f}".rstrip('0')

p_str = fmt_p(w_p_two)
shapiro_str = fmt_p(shapiro_p)

# ============================================================
# 1. SCATTER PLOT
# ============================================================
print("\nVẽ scatter plot...")

fig, ax = plt.subplots(figsize=(11, 10))

colors_scatter = [
    COLORS["win"]  if d > 0.01 else
    COLORS["lose"] if d < -0.01 else
    COLORS["tie"]
    for d in pivot["delta"]
]
ax.scatter(pivot["lin"], pivot["pan"], c=colors_scatter, s=110,
           alpha=0.75, edgecolor="white", linewidth=1.0, zorder=3)

lim = [60, 102]
ax.plot(lim, lim, "--", color="gray", alpha=0.7, linewidth=2.0,
        label="Pangenome = Linear", zorder=1)
ax.fill_between(lim, lim, [102, 102], color=COLORS["win"],
                alpha=0.06, zorder=0)

legend_elements = [
    Patch(facecolor=COLORS["win"],  alpha=0.75,
          label=f"Pangenome > Linear ({n_win})"),
    Patch(facecolor=COLORS["lose"], alpha=0.75,
          label=f"Pangenome < Linear ({n_lose})"),
    Patch(facecolor=COLORS["tie"],  alpha=0.75,
          label=f"Pangenome = Linear ({n_tie})"),
]
ax.legend(handles=legend_elements,
          loc="upper left", frameon=True,
          fontsize=20, title="%base mapping rate", title_fontsize=20)

ax.set_xlabel("Linear reference\n%base mapping rate", fontsize=24)
ax.set_ylabel("Pangenome\n%base mapping rate", fontsize=24)
ax.set_title(f"Linear vs Pangenome ({n_samples} mẫu)", fontsize=24, pad=15)
ax.set_xlim(lim)
ax.set_ylim(lim)
ax.set_aspect("equal")
ax.grid(True, alpha=0.3)

plt.tight_layout()
fig.savefig(os.path.join(OUTDIR, "fig_scatter_final.png"))
fig.savefig(os.path.join(OUTDIR, "fig_scatter_final.pdf"))
plt.close(fig)
print(f"  ✅ Scatter: {OUTDIR}/fig_scatter_final.png")

# ============================================================
# 2. BẢNG WORD — MERGE CELLS
# ============================================================
print("\nTạo bảng Word (merge cells)...")

# --- Cấu trúc dữ liệu có phân nhóm ---
table_rows = [
    # (Loại, nhãn, [Linear, Pangenome, Delta])
    ("header", "Chỉ số", ["Linear", "Pangenome", "Δ (Pan − Lin)"]),
    ("section", "%base mapping rate", ["", "", ""]),
    ("data", "  Mean (%)",             [f"{pivot['lin'].mean():.2f}",   f"{pivot['pan'].mean():.2f}",   f"{pivot['delta'].mean():+.2f}"]),
    ("data", "  Median (%)",           [f"{pivot['lin'].median():.2f}", f"{pivot['pan'].median():.2f}", f"{pivot['delta'].median():+.2f}"]),
    ("data", "  Std (%)",              [f"{pivot['lin'].std():.2f}",    f"{pivot['pan'].std():.2f}",    f"{pivot['delta'].std():.2f}"]),
    ("data", "  Min (%)",              [f"{pivot['lin'].min():.2f}",    f"{pivot['pan'].min():.2f}",    f"{pivot['delta'].min():+.2f}"]),
    ("data", "  Max (%)",              [f"{pivot['lin'].max():.2f}",    f"{pivot['pan'].max():.2f}",    f"{pivot['delta'].max():+.2f}"]),
    ("section", "So sánh Pangenome vs Linear", ["", "", ""]),
    ("data", "  Pangenome > Linear",   ["—",                            f"{n_win}/{n_samples}",         f"{n_win/n_samples*100:.1f}%"]),
    ("data", "  Pangenome < Linear",   [f"{n_lose}/{n_samples}",        "—",                            f"{n_lose/n_samples*100:.1f}%"]),
    ("data", "  Pangenome = Linear",   ["—",                            "—",                            f"{n_tie}"]),
    ("section", "Kiểm định thống kê", ["", "", ""]),
    ("data", "  Wilcoxon p-value",     ["—", "—", f"{p_str}"]),
    ("data", "  Shapiro-Wilk p-value", ["—", "—", f"{shapiro_str}"]),
    ("section", "%identity (bổ sung)", ["", "", ""]),
    ("data", "  Mean (%)",             [f"{ident_lin_mean:.2f}", f"{ident_pan_mean:.2f}", f"{ident_pan_mean - ident_lin_mean:+.2f}"]),
    ("data", "  Median (%)",           [f"{ident_lin_med:.2f}",  f"{ident_pan_med:.2f}",  f"{ident_pan_med - ident_lin_med:+.2f}"]),
]

doc = Document()

# --- Font mặc định ---
style = doc.styles["Normal"]
style.font.name = "Arial"
style.font.size = Pt(11)

# --- Caption ---
caption = doc.add_paragraph()
run = caption.add_run("Bảng 1. So sánh base mapping rate giữa linear reference và pangenome trên 106 mẫu.")
run.font.name = "Arial"
run.font.size = Pt(11)
run.bold = True

# --- Tạo bảng 4 cột ---
n_rows = len(table_rows)
table = doc.add_table(rows=n_rows, cols=4)
table.style = "Table Grid"
table.alignment = WD_TABLE_ALIGNMENT.CENTER

# --- Điền dữ liệu + style ---
section_rows = []   # lưu vị trí các section để merge
for i, (row_type, label, vals) in enumerate(table_rows):
    row_cells = table.rows[i].cells

    if row_type == "header":
        # Header
        for j, txt in enumerate([label] + vals):
            cell = row_cells[j]
            cell.text = ""
            p = cell.paragraphs[0]
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER if j > 0 else WD_ALIGN_PARAGRAPH.LEFT
            r = p.add_run(txt)
            r.bold = True
            r.font.name = "Arial"
            r.font.size = Pt(11)
            r.font.color.rgb = __import__("docx").shared.RGBColor(0xFF, 0xFF, 0xFF)
        # Nền header xám đậm
        for j in range(4):
            tc_pr = row_cells[j]._tc.get_or_add_tcPr()
            shd = OxmlElement("w:shd")
            shd.set(qn("w:val"), "clear")
            shd.set(qn("w:fill"), "404040")
            tc_pr.append(shd)

    elif row_type == "section":
        # Tiểu mục — merge 4 cột lại, in đậm, nền xám nhạt
        cell = row_cells[0].merge(row_cells[3])
        cell.text = ""
        p = cell.paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.LEFT
        r = p.add_run(label)
        r.bold = True
        r.italic = True
        r.font.name = "Arial"
        r.font.size = Pt(11)
        # Nền xám nhạt
        tc_pr = cell._tc.get_or_add_tcPr()
        shd = OxmlElement("w:shd")
        shd.set(qn("w:val"), "clear")
        shd.set(qn("w:fill"), "E8E8E8")
        tc_pr.append(shd)
        section_rows.append(i)

    else:  # data
        for j, txt in enumerate([label] + vals):
            cell = row_cells[j]
            cell.text = ""
            p = cell.paragraphs[0]
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER if j > 0 else WD_ALIGN_PARAGRAPH.LEFT
            r = p.add_run(txt)
            r.font.name = "Arial"
            r.font.size = Pt(10)
            # p-value in đậm
            if "p-value" in label:
                r.bold = True

# --- Ghi chú ---
note = doc.add_paragraph()
run = note.add_run("Ghi chú: Δ = Pangenome − Linear. Đơn vị: %. "
                   "Kiểm định Wilcoxon signed-rank test cho base mapping rate. "
                   "Identity được báo cáo bổ sung như chỉ số chất lượng alignment.")
run.font.name = "Arial"
run.font.size = Pt(9)
run.italic = True

word_path = os.path.join(OUTDIR, "table1_stats.docx")
doc.save(word_path)
print(f"  ✅ Bảng Word: {word_path}")

# ============================================================
# 3. BẢNG ẢNH — FONT ARIAL
# ============================================================
print("\nTạo bảng ảnh...")

# Dữ liệu cho bảng ảnh — dùng cấu trúc 2 chiều
plot_data = []
plot_data.append(["Chỉ số", "Linear", "Pangenome", "Δ (Pan − Lin)"])
for row_type, label, vals in table_rows[1:]:   # bỏ header
    plot_data.append([label, vals[0], vals[1], vals[2]])

fig, ax = plt.subplots(figsize=(14, 10))
ax.axis("off")

tbl = ax.table(cellText=plot_data,
               cellLoc="center",
               loc="center",
               colWidths=[0.40, 0.17, 0.18, 0.22])
tbl.auto_set_font_size(False)
tbl.set_fontsize(15)
tbl.scale(1.2, 1.9)

# --- Style ---
for i, row in enumerate(plot_data):
    for j in range(4):
        cell = tbl[(i, j)]
        cell.set_text_props(family="Arial", color="black")

    # Header
    if i == 0:
        for j in range(4):
            cell = tbl[(i, j)]
            cell.set_facecolor("#404040")
            cell.set_text_props(color="white", weight="bold",
                                family="Arial", size=15)
    # Tiểu mục (label bắt đầu không có "  ")
    elif not row[0].startswith("  ") and row[0] != "":
        for j in range(4):
            cell = tbl[(i, j)]
            cell.set_facecolor("#E8E8E8")
            cell.set_text_props(weight="bold", style="italic",
                                family="Arial", color="black")
    else:
        # Data row
        if "p-value" in row[0]:
            for j in range(4):
                tbl[(i, j)].set_facecolor("#FFF8E7")
                tbl[(i, j)].set_text_props(weight="bold",
                                           family="Arial", color="black")
        else:
            for j in range(4):
                tbl[(i, j)].set_facecolor("white")

# Left-align cột 1
for i in range(len(plot_data)):
    tbl[(i, 0)].set_text_props(ha="left", family="Arial", color="black")

ax.set_title("Bảng 1. So sánh base mapping rate giữa linear reference và pangenome",
             fontsize=18, pad=20, fontweight="bold",
             color="black", family="Arial")

plt.tight_layout()
fig.savefig(os.path.join(OUTDIR, "fig_stats_table.png"),
            bbox_inches="tight", dpi=300, facecolor="white")
fig.savefig(os.path.join(OUTDIR, "fig_stats_table.pdf"),
            bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"  ✅ Bảng ảnh: {OUTDIR}/fig_stats_table.png")

# ============================================================
# 4. IN GIÁ TRỊ 2 MẪU
# ============================================================
print("\n" + "=" * 70)
print("GIÁ TRỊ ĐỂ TỰ CHÚ THÍCH VÀO SCATTER PLOT (PowerPoint)")
print("=" * 70)

best = pivot.nlargest(1, "delta").iloc[0]
print(f"\n[1] Mẫu pangenome thắng đậm nhất")
print(f"    Sample    : {best['biosample_id']}")
print(f"    Linear    : {best['lin']:.2f}%")
print(f"    Pangenome : {best['pan']:.2f}%")
print(f"    Δ         : {best['delta']:+.2f}%")

worst = pivot.nsmallest(1, "delta").iloc[0]
print(f"\n[2] Mẫu linear thắng (âm nhất)")
print(f"    Sample    : {worst['biosample_id']}")
print(f"    Linear    : {worst['lin']:.2f}%")
print(f"    Pangenome : {worst['pan']:.2f}%")
print(f"    Δ         : {worst['delta']:+.2f}%")

target = "SAMEA8218231"
row = pivot[pivot["biosample_id"] == target]
if len(row) > 0:
    r = row.iloc[0]
    print(f"\n[3] Mẫu {target}")
    print(f"    Linear    : {r['lin']:.2f}%")
    print(f"    Pangenome : {r['pan']:.2f}%")
    print(f"    Δ         : {r['delta']:+.2f}%")

# ============================================================
# 5. TỔNG KẾT
# ============================================================
print("\n" + "=" * 70)
print("✅ HOÀN THÀNH")
print("=" * 70)
print(f"\nOutput: {OUTDIR}")
print(f"  1. fig_scatter_final.png/pdf   — Scatter")
print(f"  2. table1_stats.docx           — Bảng Word (merge cells)")
print(f"  3. fig_stats_table.png/pdf     — Bảng ảnh (Arial)")
print("=" * 70)
