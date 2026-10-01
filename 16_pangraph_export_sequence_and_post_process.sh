#!/bin/bash

# Script: pangraph_export_and_process_fixed.sh
# Description: Export sequences với pangraph, sau đó process để NC_021726.1 đầu tiên + thêm coordinates
# Usage: ./pangraph_export_and_process_fixed.sh

# ===== CẤU HÌNH =====
GRAPH_JSON="/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/graph_and_stats/pangenome_106_strains_with_refseq.json"
EXPORT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/vg_graph/block_seqs_unaligned"
FINAL_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/vg_graph/block_seqs_processed"
REF_NAME="NC_021726.1"

echo "🚀 PANGRAHP EXPORT + POST-PROCESSING"
echo "Graph JSON: $GRAPH_JSON"
echo "Export dir: $EXPORT_DIR"
echo "Final dir: $FINAL_DIR"
echo "Reference: $REF_NAME"
echo

# ===== BƯỚC 1: EXPORT VỚI PANGRAHP =====
echo "1. 📤 Export sequences với pangraph..."
pangraph export block-sequences "$GRAPH_JSON" --unaligned -o "$EXPORT_DIR"

if [ $? -ne 0 ]; then
    echo "❌ Pangraph export failed!"
    exit 1
fi

echo "✅ Pangraph export completed"

# ===== BƯỚC 2: EXTRACT COORDINATES TỪ GRAPH.JSON =====
echo "2. 🧮 Extract coordinates từ graph.json..."

python3 << EOF
import json
import os
import csv

# Cấu hình
graph_json_path = "$GRAPH_JSON"
coords_file = "$FINAL_DIR/block_coordinates.csv"
ref_name = "$REF_NAME"

# Đọc graph.json
with open(graph_json_path, 'r') as f:
    pangraph = json.load(f)

# Tìm reference path và extract coordinates
ref_coordinates = {}
for path_id, path_info in pangraph['paths'].items():
    if ref_name in path_info.get('name', ''):
        print(f"✅ Found reference: {path_info['name']}")
        
        # Extract coordinates từ reference path
        for node_id in path_info['nodes']:
            node_str = str(node_id)
            if node_str in pangraph['nodes']:
                node_info = pangraph['nodes'][node_str]
                block_id = node_info['block_id']
                block_str = str(block_id)
                
                if block_str in pangraph['blocks']:
                    position = node_info.get('position', [])
                    if len(position) == 2:
                        ref_start, ref_end = position
                        strand = node_info.get('strand', '+')
                        
                        ref_coordinates[block_str] = {
                            'ref_start': ref_start,
                            'ref_end': ref_end,
                            'strand': strand,
                            'length': ref_end - ref_start,
                            'node_id': node_id
                        }
        break

# Tạo thư mục output
os.makedirs("$FINAL_DIR", exist_ok=True)

# Ghi coordinates file
with open(coords_file, 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['block_id', 'ref_start', 'ref_end', 'strand', 'length', 'node_id'])
    for block_id, coords in ref_coordinates.items():
        writer.writerow([block_id, coords['ref_start'], coords['ref_end'], coords['strand'], coords['length'], coords['node_id']])

print(f"✅ Extracted coordinates for {len(ref_coordinates)} blocks")
print(f"📊 Coverage: {sum(c['length'] for c in ref_coordinates.values())} bp")
EOF

# ===== BƯỚC 3: PROCESS FILES BẰNG PYTHON =====
echo "3. 🔄 Processing files - Đảm bảo NC_021726.1 đầu tiên..."

python3 << EOF
import os
import glob

EXPORT_DIR = "$EXPORT_DIR"
FINAL_DIR = "$FINAL_DIR"
REF_NAME = "$REF_NAME"

# Đọc coordinates
coordinates = {}
coords_file = "$FINAL_DIR/block_coordinates.csv"
if os.path.exists(coords_file):
    with open(coords_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split(',')
            if len(parts) >= 5:
                block_id, ref_start, ref_end, strand, length = parts[0], parts[1], parts[2], parts[3], parts[4]
                coordinates[block_id] = (ref_start, ref_end, strand, length)

# Process từng file
processed = 0
with_ref = 0
without_ref = 0

for input_file in glob.glob(os.path.join(EXPORT_DIR, "block_*.fa")):
    filename = os.path.basename(input_file)
    block_id = filename[6:-3]  # Remove 'block_' and '.fa'
    
    # Đọc sequences từ file
    sequences = []
    current_header = ""
    current_sequence = ""
    
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Save previous sequence
                if current_header and current_sequence:
                    sequences.append((current_header, current_sequence))
                # Start new sequence
                current_header = line.strip()
                current_sequence = ""
            else:
                current_sequence += line.strip()
        
        # Add last sequence
        if current_header and current_sequence:
            sequences.append((current_header, current_sequence))
    
    # Phân loại: NC_021726.1 vs others
    nc_sequence = None
    other_sequences = []
    
    for header, sequence in sequences:
        if REF_NAME in header:
            # Thêm coordinates vào header nếu có
            if block_id in coordinates:
                ref_start, ref_end, strand, length = coordinates[block_id]
                header = f"{header}|ref_coords={ref_start}-{ref_end}|strand={strand}|length={length}"
            nc_sequence = (header, sequence)
        else:
            other_sequences.append((header, sequence))
    
    # Ghi file mới với NC_021726.1 đầu tiên nếu có
    output_file = os.path.join(FINAL_DIR, filename)
    with open(output_file, 'w') as f:
        if nc_sequence:
            f.write(f"{nc_sequence[0]}\n")
            f.write(f"{nc_sequence[1]}\n")
            with_ref += 1
        else:
            without_ref += 1
        
        for header, sequence in other_sequences:
            f.write(f"{header}\n")
            f.write(f"{sequence}\n")
    
    processed += 1
    if processed % 100 == 0:
        print(f"   📦 Processed {processed} files...")

print(f"✅ Processing completed: {processed} files")
print(f"   - With {REF_NAME}: {with_ref}")
print(f"   - Without {REF_NAME}: {without_ref}")
EOF

# ===== BƯỚC 4: THỐNG KÊ =====
echo "4. 📊 Thống kê kết quả:"
total_files=$(find "$FINAL_DIR" -name "block_*.fa" | wc -l)
echo "   Tổng files trong final dir: $total_files"

# Kiểm tra kết quả
echo
echo "🔍 Kiểm tra final files:"
find "$FINAL_DIR" -name "block_*.fa" | head -n 2 | while read file; do
    echo "   $(basename "$file"):"
    head -n 2 "$file"
    echo
done

echo "🎉 HOÀN THÀNH! Files đã sẵn sàng cho MAFFT:"
echo "   - NC_021726.1 luôn đầu tiên khi có"
echo "   - Có reference coordinates trong header"
echo "   - Strain-specific sequences được giữ nguyên"
echo "   - Final directory: $FINAL_DIR"
