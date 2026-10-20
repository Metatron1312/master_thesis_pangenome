import json
import re
import sys
from collections import defaultdict

# ============================================================
# CẤU HÌNH
# ============================================================

GRAPH_JSON = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/graph_and_stats/pangenome_106_strains_with_refseq.json"
OUTPUT_DIR = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/iqtree/core_genome"

# ============================================================
# BƯỚC 1: ĐỌC FILE JSON
# ============================================================

print("=" * 60)
print("BƯỚC 1: ĐỌC FILE JSON")
print("=" * 60)

with open(GRAPH_JSON, 'r') as f:
    data = json.load(f)

paths = data['paths']
blocks = data['blocks']
nodes = data['nodes']

print(f"Số path (contig): {len(paths)}")
print(f"Số block: {len(blocks)}")
print(f"Số node: {len(nodes)}")

# ============================================================
# BƯỚC 2: NHÓM PATH THEO SAMPLE ID
# ============================================================

print("\n" + "=" * 60)
print("BƯỚC 2: NHÓM PATH THEO SAMPLE ID")
print("=" * 60)

sample_to_paths = defaultdict(list)

for path_id, path in paths.items():
    desc = path.get('desc', '') or ''
    name = path.get('name', '')
    
    # Trích xuất sample ID từ desc
    match = re.search(r'biosample_id=([^\s]+)', desc)
    if match:
        sample_id = match.group(1)
    else:
        # Nếu không có biosample_id, dùng name (cho reference)
        sample_id = name
    
    sample_to_paths[sample_id].append(path_id)

n_samples = len(sample_to_paths)
print(f"Số sample ID: {n_samples}")
print(f"\n10 sample ID đầu tiên:")
for i, (sample_id, path_ids) in enumerate(list(sample_to_paths.items())[:10]):
    print(f"  {sample_id}: {len(path_ids)} contig")

# ============================================================
# BƯỚC 3: XÁC ĐỊNH CORE BLOCK
# ============================================================

print("\n" + "=" * 60)
print("BƯỚC 3: XÁC ĐỊNH CORE BLOCK")
print("=" * 60)

# Tạo dictionary: block_id -> set of sample_ids
block_to_samples = defaultdict(set)

for node_id, node in nodes.items():
    block_id = node.get('block_id')
    path_id = node.get('path_id')
    
    if block_id is None or path_id is None:
        continue
    
    # Tìm sample ID của path này
    path = paths.get(str(path_id))
    if path is None:
        continue
    
    desc = path.get('desc', '') or ''
    name = path.get('name', '')
    
    match = re.search(r'biosample_id=([^\s]+)', desc)
    if match:
        sample_id = match.group(1)
    else:
        sample_id = name
    
    block_to_samples[block_id].add(sample_id)

print(f"Số block có trong block_to_samples: {len(block_to_samples)}")

# Tìm core block: block có mặt trong TẤT CẢ sample
core_blocks = []
for block_id, samples in block_to_samples.items():
    if len(samples) == n_samples:
        core_blocks.append(block_id)

print(f"\nSố core block (có mặt trong tất cả {n_samples} sample): {len(core_blocks)}")

if len(core_blocks) == 0:
    print("\n⚠️  KHÔNG có core block nào!")
    print("Có thể do:")
    print("  1. Một số sample không có block nào chung")
    print("  2. Lỗi trong cách nhóm sample ID")
    sys.exit(1)

# ============================================================
# BƯỚC 4: THỐNG KÊ CORE BLOCK
# ============================================================

print("\n" + "=" * 60)
print("BƯỚC 4: THỐNG KÊ CORE BLOCK")
print("=" * 60)

# Tính tổng chiều dài core block
total_core_length = 0
for block_id in core_blocks:
    block = blocks.get(str(block_id))
    if block:
        total_core_length += len(block.get('consensus', ''))

print(f"Tổng chiều dài core block: {total_core_length:,} bp")
print(f"Chiều dài trung bình: {total_core_length // len(core_blocks):,} bp")

# In 10 core block đầu tiên
print(f"\n10 core block đầu tiên:")
for i, block_id in enumerate(core_blocks[:10]):
    block = blocks.get(str(block_id))
    if block:
        consensus_len = len(block.get('consensus', ''))
        n_nodes = len(block.get('alignments', {}))
        print(f"  Block {block_id}: {consensus_len} bp, {n_nodes} node")

# ============================================================
# BƯỚC 5: LƯU DANH SÁCH CORE BLOCK
# ============================================================

print("\n" + "=" * 60)
print("BƯỚC 5: LƯU DANH SÁCH CORE BLOCK")
print("=" * 60)

import os
os.makedirs(OUTPUT_DIR, exist_ok=True)

output_file = os.path.join(OUTPUT_DIR, "core_blocks.txt")
with open(output_file, "w") as f:
    for block_id in core_blocks:
        f.write(f"{block_id}\n")

print(f"Đã lưu danh sách core block vào: {output_file}")
print(f"Số core block: {len(core_blocks)}")

# Lưu thêm thông tin chi tiết
output_detail = os.path.join(OUTPUT_DIR, "core_blocks_detail.tsv")
with open(output_detail, "w") as f:
    f.write("block_id\tconsensus_len\tn_nodes\tn_samples\n")
    for block_id in core_blocks:
        block = blocks.get(str(block_id))
        if block:
            consensus_len = len(block.get('consensus', ''))
            n_nodes = len(block.get('alignments', {}))
            n_samples_block = len(block_to_samples[block_id])
            f.write(f"{block_id}\t{consensus_len}\t{n_nodes}\t{n_samples_block}\n")

print(f"Đã lưu thông tin chi tiết vào: {output_detail}")

# ============================================================
# KIỂM TRA CUỐI
# ============================================================

print("\n" + "=" * 60)
print("KIỂM TRA CUỐI")
print("=" * 60)
print(f"Số sample: {n_samples}")
print(f"Số core block: {len(core_blocks)}")
print(f"Tổng chiều dài core: {total_core_length:,} bp")
print(f"File danh sách: {output_file}")
print(f"File chi tiết: {output_detail}")