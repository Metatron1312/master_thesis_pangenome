#!/bin/bash

# Script: align_blocks_with_coordinates.sh
# Description: Align blocks với reference coordinates information
# Usage: ./align_blocks_with_coordinates.sh

# ===== CẤU HÌNH =====
INPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/vg_graph/unaligned_block_seqs_processed"
OUTPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/vg_graph/aligned_blocks_with_coords"
COORDS_FILE="$INPUT_DIR/block_coordinates.csv"
DOCKER_IMAGE="staphb/mafft:latest"
MAFFT_THREADS=8

# ===== TẠO THƯ MỤC OUTPUT =====
echo "Tạo thư mục output: $OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

# Copy coordinates file sang output dir
cp "$COORDS_FILE" "$OUTPUT_DIR/"
echo "📋 Copied coordinates file to output directory"

LOG_FILE="$OUTPUT_DIR/alignment_log.txt"
STATS_FILE="$OUTPUT_DIR/alignment_stats.csv"

echo "=== MAFFT ALIGNMENT WITH COORDINATES - $(date) ===" > "$LOG_FILE"
echo "block_id,seq_count,has_ref_coords,alignment_status,output_file" > "$STATS_FILE"

# ===== HÀM KIỂM TRA =====
count_sequences() {
    local file="$1"
    grep -c "^>" "$file"
}

has_reference_coords() {
    local file="$1"
    grep -q "ref_coords=" "$file"
}

convert_to_uppercase() {
    local input_file="$1"
    local output_file="$2"
    
    awk '
    /^>/ {
        if (seq != "") {
            print toupper(seq)
            seq = ""
        }
        print
        next
    }
    {
        seq = seq $0
    }
    END {
        if (seq != "") {
            print toupper(seq)
        }
    }' "$input_file" > "$output_file"
}

# ===== HÀM XỬ LÝ FILE =====
process_file() {
    local input_file="$1"
    local file_name=$(basename "$input_file")
    local base_name="${file_name%.fa}"
    local block_id="${base_name#block_}"
    local temp_error=$(mktemp)
    
    echo "⚙️  Xử lý: $file_name"
    
    # Kiểm tra số sequence
    seq_count=$(count_sequences "$input_file")
    has_coords=$(has_reference_coords "$input_file" && echo "yes" || echo "no")
    
    echo "   Số sequence: $seq_count, Có ref coords: $has_coords"
    
    # TRƯỜNG HỢP 1: File chỉ có 1 sequence
    if [ "$seq_count" -eq 1 ]; then
        output_file="$OUTPUT_DIR/${base_name}_single.fasta"
        cp "$input_file" "$output_file"
        echo "ℹ️  Single sequence: $file_name"
        echo "$block_id,$seq_count,$has_coords,single,$output_file" >> "$STATS_FILE"
        return 0
    fi
    
    # TRƯỜNG HỢP 2: Nhiều sequences - Align với MAFFT
    echo "   Aligning $seq_count sequences..."
    
    docker run --platform linux/amd64 --rm \
        -v "$INPUT_DIR":/input \
        -v "$OUTPUT_DIR":/output \
        $DOCKER_IMAGE \
        mafft --auto --quiet --thread $MAFFT_THREADS "/input/$file_name" > "$OUTPUT_DIR/${base_name}_aligned.fasta" 2>"$temp_error"
    
    # Kiểm tra kết quả
    if [ $? -eq 0 ] && [ -s "$OUTPUT_DIR/${base_name}_aligned.fasta" ]; then
        # Chuyển sang chữ hoa
        convert_to_uppercase "$OUTPUT_DIR/${base_name}_aligned.fasta" "${OUTPUT_DIR}/${base_name}_aligned.fasta.tmp"
        mv "${OUTPUT_DIR}/${base_name}_aligned.fasta.tmp" "$OUTPUT_DIR/${base_name}_aligned.fasta"
        
        final_count=$(count_sequences "$OUTPUT_DIR/${base_name}_aligned.fasta")
        
        if [ "$final_count" -eq "$seq_count" ]; then
            echo "✅ Aligned: $file_name ($seq_count sequences)"
            echo "$block_id,$seq_count,$has_coords,aligned,${base_name}_aligned.fasta" >> "$STATS_FILE"
        else
            echo "⚠️  Sequence count changed: $file_name ($seq_count → $final_count)"
            echo "$block_id,$seq_count,$has_coords,warning_count_changed,${base_name}_aligned.fasta" >> "$STATS_FILE"
        fi
        
        # Giữ nguyên reference coordinates trong header
        if [ "$has_coords" = "yes" ]; then
            echo "   ✓ Giữ nguyên reference coordinates"
        fi
        
        rm -f "$temp_error"
        return 0
    else
        echo "❌ MAFFT failed: $file_name"
        echo "FAILED: $file_name" >> "$LOG_FILE"
        cat "$temp_error" >> "$LOG_FILE"
        
        # Fallback: copy file gốc
        output_file="$OUTPUT_DIR/${base_name}_unaligned.fasta"
        cp "$input_file" "$output_file"
        echo "$block_id,$seq_count,$has_coords,failed,${base_name}_unaligned.fasta" >> "$STATS_FILE"
        
        rm -f "$temp_error" "$OUTPUT_DIR/${base_name}_aligned.fasta" 2>/dev/null
        return 1
    fi
}

# ===== XỬ LÝ CHÍNH =====
echo "🚀 BẮT ĐẦU ALIGNMENT VỚI REFERENCE COORDINATES"
echo "Input: $INPUT_DIR"
echo "Output: $OUTPUT_DIR"
echo "MAFFT threads: $MAFFT_THREADS"
echo

START_TIME=$(date +%s)

processed=0
aligned=0
single=0
failed=0
with_coords=0
no_coords=0

total_files=$(find "$INPUT_DIR" -name "block_*.fa" | wc -l)
echo "📊 Tổng số block files: $total_files"

# Xử lý từng file
for input_file in "$INPUT_DIR"/block_*.fa; do
    if [ -f "$input_file" ]; then
        process_file "$input_file"
        result=$?
        ((processed++))
        
        # Thống kê
        if [ $result -eq 0 ]; then
            if grep -q "_single.fasta" <<< "$(ls "$OUTPUT_DIR"/${input_file##*/}* 2>/dev/null)"; then
                ((single++))
            else
                ((aligned++))
            fi
            
            # Kiểm tra có coordinates không
            if has_reference_coords "$input_file"; then
                ((with_coords++))
            else
                ((no_coords++))
            fi
        else
            ((failed++))
        fi
        
        # Hiển thị progress
        if [ $((processed % 50)) -eq 0 ]; then
            elapsed=$(( $(date +%s) - START_TIME ))
            progress_percent=$(( processed * 100 / total_files ))
            echo "📊 Progress: $processed/$total_files ($progress_percent%)"
            echo "   Aligned: $aligned | Single: $single | Failed: $failed"
            echo "   With coords: $with_coords | No coords: $no_coords"
        fi
    fi
done

END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))

# ===== THỐNG KÊ KẾT QUẢ =====
echo
echo "📊 KẾT QUẢ ALIGNMENT:"
echo "   Tổng blocks: $total_files"
echo "   Aligned thành công: $aligned"
echo "   Single sequences: $single"
echo "   Failed: $failed"
echo "   Blocks có reference coordinates: $with_coords"
echo "   Blocks không có coordinates: $no_coords"
echo "   Thời gian: $((TOTAL_TIME / 3600))h $(( (TOTAL_TIME % 3600) / 60 ))m $((TOTAL_TIME % 60))s"

# Ghi summary
echo "" >> "$LOG_FILE"
echo "=== SUMMARY ===" >> "$LOG_FILE"
echo "Total blocks: $total_files" >> "$LOG_FILE"
echo "Successfully aligned: $aligned" >> "$LOG_FILE"
echo "Single sequences: $single" >> "$LOG_FILE"
echo "Failed: $failed" >> "$LOG_FILE"
echo "With reference coordinates: $with_coords" >> "$LOG_FILE"
echo "Without coordinates: $no_coords" >> "$LOG_FILE"
echo "Total time: ${TOTAL_TIME}s" >> "$LOG_FILE"

echo
echo "🎉 ALIGNMENT HOÀN THÀNH!"
echo "📁 Output: $OUTPUT_DIR"
echo "📊 Stats: $STATS_FILE"
echo "📋 Log: $LOG_FILE"
echo
echo "NEXT STEP: Chạy MSA-to-VCF script với reference coordinates"
