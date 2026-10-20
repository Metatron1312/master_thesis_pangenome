import os
import re
import json
from collections import defaultdict

# ============================================================
# CẤU HÌNH
# ============================================================

GRAPH_JSON = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/graph_and_stats/pangenome_106_strains_with_refseq.json"
PROCESSED_DIR = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/vg_graph/unaligned_block_seqs_processed"
CORE_BLOCKS = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/iqtree/core_genome/core_blocks.txt"
OUTPUT_DIR = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/iqtree/core_genome/core_block_seqs"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# BƯỚC 1: ĐỌC GRAPH JSON
# ============================================================

print("=" * 60)
print("BƯỚC 1: ĐỌC GRAPH JSON")
print("=" * 60)

with open(GRAPH_JSON) as f:
    data = json.load(f)

paths = data['paths']
blocks = data['blocks']
nodes = data['nodes']

print(f"Số path: {len(paths)}")
print(f"Số block: {len(blocks)}")
print(f"Số node: {len(nodes)}")

# ============================================================
# BƯỚC 2: TẠO ÁNH XẠ PATH_ID → SAMPLE_ID
# ============================================================

print("\n" + "=" * 60)
print("BƯỚC 2: TẠO ÁNH XẠ PATH_ID → SAMPLE_ID")
print("=" * 60)

path_to_sample = {}

for path_id, path in paths.items():
    desc = path.get('desc', '') or ''
    name = path.get('name', '')
    
    # Trích xuất biosample_id từ desc
    match = re.search(r'biosample_id=([^\s]+)', desc)
    if match:
        sample_id = match.group(1)
    else:
        # Nếu không có biosample_id, dùng name
        sample_id = name
    
    path_to_sample[path_id] = sample_id

# Thống kê
sample_counts = defaultdict(int)
for sample_id in path_to_sample.values():
    sample_counts[sample_id] += 1

print(f"Số sample ID unique: {len(sample_counts)}")
print(f"\n10 sample ID đầu tiên:")
for i, (sample_id, count) in enumerate(list(sample_counts.items())[:10]):
    print(f"  {sample_id}: {count} path")

# ============================================================
# BƯỚC 3: TẠO ÁNH XẠ NODE_ID → SAMPLE_ID
# ============================================================

print("\n" + "=" * 60)
print("BƯỚC 3: TẠO ÁNH XẠ NODE_ID → SAMPLE_ID")
print("=" * 60)

node_to_sample = {}

for node_id, node in nodes.items():
    path_id = node.get('path_id')
    if path_id is None:
        continue
    
    path_id_str = str(path_id)
    if path_id_str in path_to_sample:
        node_to_sample[node_id] = path_to_sample[path_id_str]

print(f"Số node có sample ID: {len(node_to_sample)}")

# ============================================================
# BƯỚC 4: ĐỌC DANH SÁCH CORE BLOCK
# ============================================================

print("\n" + "=" * 60)
print("BƯỚC 4: ĐỌC DANH SÁCH CORE BLOCK")
print("=" * 60)

with open(CORE_BLOCKS) as f:
    core_blocks = [line.strip() for line in f if line.strip()]

print(f"Số core block: {len(core_blocks)}")

# ============================================================
# BƯỚC 5: XỬ LÝ TỪNG CORE BLOCK
# ============================================================

print("\n" + "=" * 60)
print("BƯỚC 5: XỬ LÝ TỪNG CORE BLOCK")
print("=" * 60)

stats = {
    'total_blocks': 0,
    'blocks_with_107': 0,
    'blocks_with_less': 0,
    'blocks_with_more': 0,
    'skipped': 0,
}

for i, block_id in enumerate(core_blocks):
    input_file = os.path.join(PROCESSED_DIR, f"block_{block_id}.fa")
    
    if not os.path.exists(input_file):
        stats['skipped'] += 1
        continue
    
    # Đọc tất cả sequence
    sequences = []
    current_header = None
    current_seq = ""
    
    with open(input_file) as f:
        for line in f:
            if line.startswith('>'):
                if current_header and current_seq:
                    sequences.append((current_header, current_seq))
                current_header = line.strip()
                current_seq = ""
            else:
                current_seq += line.strip()
        if current_header and current_seq:
            sequences.append((current_header, current_seq))
    
    # Nhóm theo sample ID, chọn sequence đầu tiên cho mỗi sample
    sample_to_seq = {}
    for header, seq in sequences:
        # Trích xuất node_id từ header
        # Header có dạng: >14811260016675454272 {"path_name":...}
        match = re.match(r'>(\d+)', header)
        if not match:
            continue
        node_id = match.group(1)
        
        # Tra cứu sample ID
        if node_id in node_to_sample:
            sample_id = node_to_sample[node_id]
        else:
            # Fallback: dùng path_name
            match2 = re.search(r'"path_name":"([^"]+)"', header)
            if match2:
                path_name = match2.group(1)
                if '.' in path_name:
                    sample_id = path_name.split('.')[0]
                else:
                    sample_id = path_name
            else:
                continue
        
        # Chỉ lấy sequence đầu tiên cho mỗi sample
        if sample_id not in sample_to_seq:
            sample_to_seq[sample_id] = (header, seq)
    
    n_samples = len(sample_to_seq)
    
    # Ghi file output
    output_file = os.path.join(OUTPUT_DIR, f"block_{block_id}.fa")
    with open(output_file, "w") as f:
        for sample_id, (header, seq) in sample_to_seq.items():
            f.write(f">{sample_id}\n{seq}\n")
    
    stats['total_blocks'] += 1
    if n_samples == 107:
        stats['blocks_with_107'] += 1
    elif n_samples < 107:
        stats['blocks_with_less'] += 1
    else:
        stats['blocks_with_more'] += 1
    
    if (i + 1) % 100 == 0:
        print(f"  Đã xử lý {i + 1}/{len(core_blocks)} block")

print(f"\nHoàn tất xử lý!")
print(f"  Tổng block xử lý: {stats['total_blocks']}")
print(f"  Block có đúng 107 sample: {stats['blocks_with_107']}")
print(f"  Block có < 107 sample: {stats['blocks_with_less']}")
print(f"  Block có > 107 sample: {stats['blocks_with_more']}")
print(f"  Block bị bỏ qua: {stats['skipped']}")

# ============================================================
# BƯỚC 6: KIỂM TRA KẾT QUẢ
# ============================================================

print("\n" + "=" * 60)
print("BƯỚC 6: KIỂM TRA KẾT QUẢ")
print("=" * 60)

sample_block = core_blocks[0]
sample_output = os.path.join(OUTPUT_DIR, f"block_{sample_block}.fa")

print(f"\nFile mẫu: {sample_output}")
if os.path.exists(sample_output):
    with open(sample_output) as f:
        headers = [line.strip() for line in f if line.startswith('>')]
    print(f"Số sequence: {len(headers)}")
    print(f"5 header đầu tiên:")
    for h in headers[:5]:
        print(f"  {h}")

print(f"\nOutput dir: {OUTPUT_DIR}")
n_files = len([f for f in os.listdir(OUTPUT_DIR) if f.endswith('.fa')])
print(f"Số file output: {n_files}")
