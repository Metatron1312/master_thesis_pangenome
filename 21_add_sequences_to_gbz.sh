#!/bin/bash

# Cấu hình
EXISTING_GBZ="/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/vg_graph/graph/vg_autoindex/pangenome_106_strains_with_refseq.giraffe.gbz"
FASTA_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/vg_graph/unaligned_block_seqs_processed"
VCF_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/vg_graph/msa2vcf_v19/processed"
OUTPUT_PREFIX="complete_106_strains_with_refseq_pangenome"
BATCH_SIZE=5  # Giảm batch size xuống 5

echo "🏗️ SCALING UP: Building complete pangenome graph (Memory-optimized)..."

# Bước 1: Extract base graph từ GBZ
echo "1. 🔄 Extracting base graph from GBZ..."
vg convert "$EXISTING_GBZ" > ${OUTPUT_PREFIX}.vg
ORIGINAL_NODES=$(vg stats -l ${OUTPUT_PREFIX}.vg | grep -o '[0-9]\+' | head -1)
echo "   📏 Base graph: $ORIGINAL_NODES nodes"

# Bước 2: Xác định sequences độc nhất
echo "2. 🔍 Identifying unique sequences..."
find "$VCF_DIR" -name "block_*.vcf.gz" | sed 's/.*block_//; s/.vcf.gz//' > existing_blocks.txt

find "$FASTA_DIR" -name "block_*.fa" | while read fa_file; do
    block_id=$(basename "$fa_file" .fa | sed 's/block_//')
    if ! grep -q "^$block_id$" existing_blocks.txt; then
        echo "$fa_file"
    fi
done > valid_unique_fasta.txt

TOTAL_FILES=$(wc -l < valid_unique_fasta.txt)
echo "   📊 Total unique sequences: $TOTAL_FILES"

if [[ $TOTAL_FILES -eq 0 ]]; then
    echo "❌ No unique sequences found"
    exit 1
fi

# Bước 3: Xử lý từng file một (tránh memory issues)
echo "3. 🏗️ Processing sequences one by one..."
PROCESSED=0
ADDED_COUNT=0
FAILED_COUNT=0

while read fa_file; do
    ((PROCESSED++))
    echo "   ➕ [$PROCESSED/$TOTAL_FILES] Processing: $(basename "$fa_file")"
    
    # Tạo graph từ file đơn lẻ
    vg construct -r "$fa_file" -m 32 > single_graph.vg 2>/dev/null
    
    if [[ -s single_graph.vg ]]; then
        SINGLE_NODES=$(vg stats -l single_graph.vg | grep -o '[0-9]\+' | head -1)
        echo "     📏 Single graph: $SINGLE_NODES nodes"
        
        # Combine với graph chính
        vg combine ${OUTPUT_PREFIX}.vg single_graph.vg > combined_graph.vg 2>/dev/null
        
        if [[ -s combined_graph.vg ]]; then
            mv combined_graph.vg ${OUTPUT_PREFIX}.vg
            ((ADDED_COUNT++))
            echo "     ✅ Added successfully"
        else
            echo "     ❌ Combine failed"
            ((FAILED_COUNT++))
        fi
    else
        echo "     ❌ Failed to create graph from file"
        ((FAILED_COUNT++))
    fi
    
    # Hiển thị progress mỗi 50 files
    if [[ $((PROCESSED % 50)) -eq 0 ]]; then
        CURRENT_NODES=$(vg stats -l ${OUTPUT_PREFIX}.vg | grep -o '[0-9]\+' | head -1)
        echo "   📈 Progress: $PROCESSED/$TOTAL_FILES files, $ADDED_COUNT added, $CURRENT_NODES total nodes"
    fi
    
    # Dọn dẹp
    rm -f single_graph.vg combined_graph.vg
    
done < valid_unique_fasta.txt

# Bước 4: Finalize
echo "4. 📊 Finalizing complete graph..."
FINAL_NODES=$(vg stats -l ${OUTPUT_PREFIX}.vg | grep -o '[0-9]\+' | head -1)
FINAL_EDGES=$(vg stats -z ${OUTPUT_PREFIX}.vg | grep edges | grep -o '[0-9]\+')

echo "🎉 PROCESSING COMPLETE!"
echo "📊 Results:"
echo "   - Original graph: $ORIGINAL_NODES nodes"
echo "   - Final graph: $FINAL_NODES nodes, $FINAL_EDGES edges"
echo "   - Processed files: $PROCESSED/$TOTAL_FILES"
echo "   - Successfully added: $ADDED_COUNT"
echo "   - Failed: $FAILED_COUNT"

# Bước 5: Index graph cuối cùng
if [[ $ADDED_COUNT -gt 0 ]]; then
    echo "5. 🔨 Indexing final graph..."
    vg index ${OUTPUT_PREFIX}.vg -x ${OUTPUT_PREFIX}.xg
    
    if [[ -f "${OUTPUT_PREFIX}.xg" ]]; then
        echo "   ✅ XG index created: ${OUTPUT_PREFIX}.xg"
        
        # Tạo GBWT từ graph hoàn chỉnh
        echo "6. 📊 Building GBWT from complete graph..."
        vg gbwt -x ${OUTPUT_PREFIX}.xg -E -o ${OUTPUT_PREFIX}.gbwt -p
        if [[ -f "${OUTPUT_PREFIX}.gbwt" ]]; then
            echo "   ✅ GBWT created: ${OUTPUT_PREFIX}.gbwt"
            echo "   📊 GBWT info:"
            vg gbwt -M ${OUTPUT_PREFIX}.gbwt
        fi
    fi
fi

# Dọn dẹp
rm -f existing_blocks.txt valid_unique_fasta.txt

echo "🧹 Cleanup completed"
echo "📁 Final outputs:"
ls -lh ${OUTPUT_PREFIX}.*
