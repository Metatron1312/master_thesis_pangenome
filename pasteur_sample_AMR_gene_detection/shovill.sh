#!/bin/bash

# =============================================
# CẤU HÌNH ĐƯỜNG DẪN
INPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/Assess_SNP_SV/data/scale_up/filtered_reads"
OUTPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/Assess_SNP_SV/data/scale_up/assembled"
REPORT_FILE="$OUTPUT_DIR/assembly_report_$(date +%Y%m%d).txt"
# =============================================

# Hàm kiểm tra RAM thực tế (không cần Alpine)
get_safe_ram() {
    # Sử dụng giá trị mặc định 4GB để đảm bảo ổn định
    echo 4
}

# Hàm tạo báo cáo
init_report() {
    echo "SHOVILL ASSEMBLY REPORT - $(date)" > "$REPORT_FILE"
    echo "=================================" >> "$REPORT_FILE"
    printf "%-30s | %-8s | %-10s | %s\n" "Sample" "Contigs" "MaxLen" "Status" >> "$REPORT_FILE"
    echo "------------------------------------------------------------" >> "$REPORT_FILE"
}

add_to_report() {
    local sample=$1
    local contigs=$2
    local maxlen=$3
    local status=$4
    printf "%-30s | %-8s | %-10s | %s\n" "$sample" "$contigs" "$maxlen" "$status" >> "$REPORT_FILE"
}

# Kiểm tra Docker
if ! command -v docker &> /dev/null; then
    echo "Error: Docker not found!" | tee -a "$REPORT_FILE"
    exit 1
fi

# Khởi tạo báo cáo
mkdir -p "$OUTPUT_DIR" || { echo "Error: Cannot create output directory"; exit 1; }
init_report

# Hàm phân tích kết quả assembly
analyze_assembly() {
    local contigs_file="$1"
    if [ ! -f "$contigs_file" ]; then
        echo "0 0"
        return
    fi
    awk '/^>/ {if (seqlen){print seqlen}; seqlen=0; next} {seqlen+=length($0)} END {print seqlen}' "$contigs_file" | 
    sort -nr | 
    awk 'BEGIN {max=0; count=0} 
         {if($1>max) max=$1; count++} 
         END {print count, max}'
}

# Hàm chạy Shovill
run_shovill() {
    local R1_FILE="$1"
    local R2_FILE="$2"
    local SAMPLE_NAME="$3"
    local SAMPLE_OUTDIR="$OUTPUT_DIR/$SAMPLE_NAME"
    
    echo "🔍 Processing: $SAMPLE_NAME"
    
    # Xóa thư mục cũ nếu tồn tại
    if [ -d "$SAMPLE_OUTDIR" ]; then
        echo "   - Removing previous output directory..."
        rm -rf "$SAMPLE_OUTDIR"
    fi
    
    mkdir -p "$SAMPLE_OUTDIR" || { echo "❌ Cannot create output directory"; return 1; }

    local AVAILABLE_RAM=$(get_safe_ram)
    local CPU_CORES=4  # Sử dụng 4 cores để ổn định
    
    echo "   - Using $CPU_CORES CPU cores"
    echo "   - Using $AVAILABLE_RAM GB RAM"

    # Tham số tối ưu cho A. baumannii
    local SHOVILL_OPTS=(
        --assembler spades
        --trim
        --cpus "$CPU_CORES"
        --ram "$AVAILABLE_RAM"
        --force  # Thêm --force để ghi đè thư mục
    )

    echo "   - Starting assembly..."
    
    # Chạy Shovill trong Docker
    docker run --rm -it \
        --platform linux/amd64 \
        --memory "${AVAILABLE_RAM}g" \
        --cpus "$CPU_CORES" \
        -v "$INPUT_DIR":/input \
        -v "$SAMPLE_OUTDIR":/output \
        staphb/shovill \
        shovill \
            --R1 "/input/$(basename "$R1_FILE")" \
            --R2 "/input/$(basename "$R2_FILE")" \
            --outdir /output \
            "${SHOVILL_OPTS[@]}"
    
    # Phân tích kết quả
    if [ -f "$SAMPLE_OUTDIR/contigs.fa" ]; then
        local stats=($(analyze_assembly "$SAMPLE_OUTDIR/contigs.fa"))
        add_to_report "$SAMPLE_NAME" "${stats[0]}" "${stats[1]}" "Success"
        echo "✅ Success: ${stats[0]} contigs (maxlen=${stats[1]} bp)"
        return 0
    else
        add_to_report "$SAMPLE_NAME" "0" "0" "Failed"
        echo "❌ Assembly failed"
        return 1
    fi
}

# Main pipeline
echo "Starting Shovill Assembly Pipeline..."
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Report file: $REPORT_FILE"
echo "----------------------------------------"

TOTAL=0
SUCCESS=0

for R1_FILE in "$INPUT_DIR"/*_R1.fastq.gz; do
    SAMPLE_NAME=$(basename "$R1_FILE" | sed 's/_R1.fastq.gz//')
    R2_FILE="${R1_FILE/_R1.fastq.gz/_R2.fastq.gz}"
    
    ((TOTAL++))
    
    if [ ! -f "$R2_FILE" ]; then
        add_to_report "$SAMPLE_NAME" "0" "0" "Missing R2"
        echo "⚠️ Missing R2 for: $SAMPLE_NAME"
        continue
    fi
    
    echo "------------------------------"
    if run_shovill "$R1_FILE" "$R2_FILE" "$SAMPLE_NAME"; then
        ((SUCCESS++))
    fi
    echo "------------------------------"
done

# Final report
echo "================================="
echo "SUMMARY: $SUCCESS/$TOTAL successful"
echo "Full report saved to: $REPORT_FILE"

# Hiển thị báo cáo cuối cùng
echo ""
echo "FINAL REPORT:"
cat "$REPORT_FILE"

exit 0