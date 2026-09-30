#!/usr/bin/env python3
"""
Phân tích pangenome: core / unique / accessory
+ phân loại chromosome_specific / plasmid_specific / shared / reference_only

Quy tắc:
- n_samples = TẤT CẢ sample (bao gồm reference)
- Core       = có mặt ở đủ n_samples
- Unique     = có mặt ở đúng 1 sample
- Accessory  = có mặt ở 2 .. n_samples-1
- seq_type   = chỉ mô tả, không ảnh hưởng occupancy
- reference_only = block chỉ xuất hiện trong sample reference
"""

import json
import re
import os
from collections import defaultdict
import pandas as pd

# ============================================================
# CẤU HÌNH
# ============================================================
GRAPH_JSON = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/graph_and_stats/pangenome_106_strains_with_refseq.json"
OUTPUT_DIR = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/graph_and_stats/stats_official"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# BƯỚC 1: ĐỌC JSON
# ============================================================
print("=" * 60)
print("BƯỚC 1: ĐỌC JSON")
print("=" * 60)

with open(GRAPH_JSON, "r") as f:
    data = json.load(f)

paths  = data["paths"]
blocks = data["blocks"]
nodes  = data["nodes"]

print(f"Số path  : {len(paths)}")
print(f"Số block : {len(blocks)}")
print(f"Số node  : {len(nodes)}")

# ============================================================
# BƯỚC 2: HÀM TRÍCH XUẤT
# ============================================================
def get_sample_id(path):
    """Ưu tiên biosample_id trong desc, fallback về name."""
    desc = path.get("desc") or ""
    name = path.get("name") or ""
    m = re.search(r"biosample_id=([^\s,;]+)", desc)
    if m:
        return m.group(1)
    return name if name else "UNKNOWN"

def get_seq_type(path):
    """chromosome / plasmid / unknown — CHỈ để mô tả."""
    text = ((path.get("desc") or "") + " " + (path.get("name") or "")).lower()
    if "plasmid" in text:
        return "plasmid"
    if "chromosome" in text or "chrom" in text:
        return "chromosome"
    return "unknown"

# ============================================================
# BƯỚC 3: NHÓM PATH THEO SAMPLE + PHÁT HIỆN REFERENCE
# ============================================================
print("\n" + "=" * 60)
print("BƯỚC 3: NHÓM PATH THEO SAMPLE")
print("=" * 60)

sample_to_paths  = defaultdict(list)
path_to_sample   = {}
path_to_seqtype  = {}
reference_samples = set()   # sample KHÔNG có biosample_id

for path_id, path in paths.items():
    sid   = get_sample_id(path)
    stype = get_seq_type(path)
    sample_to_paths[sid].append(path_id)
    path_to_sample[path_id]  = sid
    path_to_seqtype[path_id] = stype

    desc = path.get("desc") or ""
    if "biosample_id=" not in desc:
        reference_samples.add(sid)

n_samples = len(sample_to_paths)
print(f"→ n_samples (gồm cả reference) = {n_samples}")
print(f"→ Reference samples            = {sorted(reference_samples)}")

# ============================================================
# BƯỚC 4: XÂY DỰNG BLOCK → SAMPLES + SEQ TYPES
# ============================================================
print("\n" + "=" * 60)
print("BƯỚC 4: XÂY DỰNG BLOCK → SAMPLES")
print("=" * 60)

block_to_samples  = defaultdict(set)
block_to_seqtypes = defaultdict(set)

for node_id, node in nodes.items():
    block_id = node.get("block_id")
    path_id  = node.get("path_id")
    if block_id is None or path_id is None:
        continue

    pid = str(path_id)
    if pid not in path_to_sample:
        continue

    block_to_samples[block_id].add(path_to_sample[pid])

    stype = path_to_seqtype[pid]
    if stype != "unknown":
        block_to_seqtypes[block_id].add(stype)

block_to_samples = {b: s for b, s in block_to_samples.items() if s}
print(f"→ Số block có ít nhất 1 alignment: {len(block_to_samples)}")

# ============================================================
# BƯỚC 5: PHÂN LOẠI
# ============================================================
print("\n" + "=" * 60)
print("BƯỚC 5: PHÂN LOẠI BLOCK")
print("=" * 60)

records = []
for block_id, samples in block_to_samples.items():
    n_present = len(samples)

    # --- Occupancy ---
    if n_present == n_samples:
        occupancy = "core"
    elif n_present == 1:
        occupancy = "unique"
    else:
        occupancy = "accessory"

    # --- Sequence type ---
    stypes  = block_to_seqtypes.get(block_id, set())
    has_chr = "chromosome" in stypes
    has_pls = "plasmid"    in stypes

    if has_chr and has_pls:
        seq_type = "shared"
    elif has_chr:
        seq_type = "chromosome_specific"
    elif has_pls:
        seq_type = "plasmid_specific"
    else:
        # Không xác định được chr/pls → kiểm tra reference_only
        if n_present == 1 and next(iter(samples)) in reference_samples:
            seq_type = "reference_only"
        else:
            seq_type = "unknown"

    block  = blocks.get(str(block_id), {})
    length = len(block.get("consensus", ""))

    records.append({
        "block_id"  : block_id,
        "length"    : length,
        "n_present" : n_present,
        "occupancy" : occupancy,
        "seq_type"  : seq_type,
    })

df = pd.DataFrame(records)
print(f"→ Tổng block phân tích: {len(df)}")

# ============================================================
# BƯỚC 6: THỐNG KÊ
# ============================================================
print("\n" + "=" * 60)
print("BƯỚC 6: THỐNG KÊ")
print("=" * 60)

occ_stats = df.groupby("occupancy").agg(
    n_blocks = ("block_id", "count"),
    total_bp = ("length",   "sum"),
    avg_bp   = ("length",   "mean"),
).round(1)
print("\n[Occupancy — gồm reference trong n_samples]")
print(occ_stats)

seq_stats = df.groupby("seq_type").agg(
    n_blocks = ("block_id", "count"),
    total_bp = ("length",   "sum"),
    avg_bp   = ("length",   "mean"),
).round(1)
print("\n[Sequence type]")
print(seq_stats)

cross    = pd.crosstab(df["occupancy"], df["seq_type"])
cross_bp = df.pivot_table(
    index="occupancy", columns="seq_type",
    values="length", aggfunc="sum", fill_value=0
)
print("\n[Cross-tab block: occupancy × seq_type]")
print(cross)
print("\n[Cross-tab bp: occupancy × seq_type]")
print(cross_bp)

curve = (df.groupby("n_present")
           .agg(n_blocks=("block_id","count"), total_bp=("length","sum"))
           .reset_index().sort_values("n_present"))
curve["cum_blocks"] = curve["n_blocks"].cumsum()
curve["cum_bp"]     = curve["total_bp"].cumsum()

# ============================================================
# BƯỚC 7: XUẤT FILE
# ============================================================
print("\n" + "=" * 60)
print("BƯỚC 7: XUẤT FILE")
print("=" * 60)

df.to_csv(os.path.join(OUTPUT_DIR, "block_classification.tsv"),
          sep="\t", index=False)
occ_stats.to_csv(os.path.join(OUTPUT_DIR, "stats_occupancy.tsv"), sep="\t")
seq_stats.to_csv(os.path.join(OUTPUT_DIR, "stats_seqtype.tsv"), sep="\t")
cross.to_csv(os.path.join(OUTPUT_DIR, "crosstab_occupancy_seqtype.tsv"), sep="\t")
cross_bp.to_csv(os.path.join(OUTPUT_DIR, "crosstab_occupancy_seqtype_bp.tsv"), sep="\t")
curve.to_csv(os.path.join(OUTPUT_DIR, "occupancy_curve.tsv"), sep="\t", index=False)

for label, sub in [("core",           df[df.occupancy=="core"]),
                   ("unique",         df[df.occupancy=="unique"]),
                   ("accessory",      df[df.occupancy=="accessory"]),
                   ("reference_only", df[df.seq_type=="reference_only"])]:
    with open(os.path.join(OUTPUT_DIR, f"{label}_block_ids.txt"), "w") as f:
        f.write("\n".join(map(str, sub["block_id"])))

core_df = df[df["occupancy"] == "core"]
uni_df  = df[df["occupancy"] == "unique"]
acc_df  = df[df["occupancy"] == "accessory"]

summary = pd.DataFrame([
    {"category": "Core",               "n_blocks": len(core_df), "total_bp": int(core_df["length"].sum())},
    {"category": "Accessory (2..n-1)", "n_blocks": len(acc_df),  "total_bp": int(acc_df["length"].sum())},
    {"category": "Unique (1 sample)",  "n_blocks": len(uni_df),  "total_bp": int(uni_df["length"].sum())},
])
summary.to_csv(os.path.join(OUTPUT_DIR, "summary_for_plots.tsv"),
               sep="\t", index=False)

seqtype_block_counts = (df["seq_type"].value_counts()
                          .rename_axis("seq_type").reset_index(name="n_blocks"))
seqtype_block_counts.to_csv(os.path.join(OUTPUT_DIR, "seqtype_block_counts.tsv"),
                            sep="\t", index=False)

print(f"✅ Đã ghi tất cả file vào: {OUTPUT_DIR}")

# ============================================================
# BƯỚC 8: TÓM TẮT
# ============================================================
print("\n" + "=" * 60)
print("TÓM TẮT")
print("=" * 60)
print(f"n_samples (gồm ref)     : {n_samples}")
print(f"Tổng block              : {len(df)}")
print(f"Core                    : {len(core_df)} block  |  {int(core_df['length'].sum()):,} bp")
print(f"Accessory (2..n-1)      : {len(acc_df)} block  |  {int(acc_df['length'].sum()):,} bp")
print(f"Unique (1 sample)       : {len(uni_df)} block  |  {int(uni_df['length'].sum()):,} bp")
print()
print("Theo sequence type:")
for st, n in df["seq_type"].value_counts().items():
    bp = int(df[df.seq_type==st]["length"].sum())
    print(f"  {st:22s}: {n:5d} block  |  {bp:,} bp")
