#!/bin/bash

# ================= CẤU HÌNH ĐƯỜNG DẪN =================
BASE_INPUT_DIR="/Users/shinra/gene_detection_v2/validate_database/SAMN49741230/raw_reads"
OUTPUT_DIR="/Users/shinra/gene_detection_v2/validate_database/SAMN49741230/fastp_trimmed"
REPORT_DIR="/Users/shinra/gene_detection_v2/validate_database/SAMN49741230/fastp_trimmed/fastp_reports"

# Danh sách mẫu
SAMPLES=("SAMN49741230")
# =======================================================

# Tạo thư mục output nếu chưa tồn tại
mkdir -p "$OUTPUT_DIR" "$REPORT_DIR"

# Số luồng CPU (cho macOS)
THREADS=$(sysctl -n hw.ncpu 2>/dev/null || echo 4)

echo "🔍 Bắt đầu xử lý..."

# Kiểm tra thư mục input tồn tại
if [ ! -d "$BASE_INPUT_DIR" ]; then
    echo "❌ Không tìm thấy thư mục: $BASE_INPUT_DIR"
    echo "📂 Các thư mục có sẵn:"
    ls -la "/Users/shinra/SAMN49741230/"
    exit 1
fi

echo "▶️ Đang xử lý từ thư mục: $BASE_INPUT_DIR"

# Tìm tất cả file R1 (dạng *_1.fastq.gz)
shopt -s nullglob
R1_FILES=("$BASE_INPUT_DIR"/*_1.fastq.gz)
shopt -u nullglob

if [ ${#R1_FILES[@]} -eq 0 ]; then
    echo "❌ Không tìm thấy file R1 (*_1.fastq.gz) trong $BASE_INPUT_DIR!"
    echo "📂 Các file trong thư mục:"
    ls -la "$BASE_INPUT_DIR"
    exit 1
fi

echo "✅ Đã tìm thấy ${#R1_FILES[@]} file R1"

for SAMPLE in "${SAMPLES[@]}"; do
    for R1_FILE in "${R1_FILES[@]}"; do
        # Xác định file R2 tương ứng (thay _1 bằng _2)
        R2_FILE="${R1_FILE/_1.fastq.gz/_2.fastq.gz}"
        
        # Kiểm tra file R2 có tồn tại không
        if [ ! -f "$R2_FILE" ]; then
            echo "⚠️ Không tìm thấy file R2 tương ứng: $(basename "$R2_FILE")"
            continue
        fi
        
        echo "   - R1: $(basename "$R1_FILE")"
        echo "   - R2: $(basename "$R2_FILE")"
        echo "   - Sample: $SAMPLE"
        
        # Chạy fastp qua Docker - sử dụng tên sample cho output
        docker run --rm \
            -v "$BASE_INPUT_DIR":/input \
            -v "$OUTPUT_DIR":/output \
            -v "$REPORT_DIR":/reports \
            quay.io/biocontainers/fastp:1.0.1--heae3180_0 \
            fastp \
            -i "/input/$(basename "$R1_FILE")" \
            -I "/input/$(basename "$R2_FILE")" \
            -o "/output/${SAMPLE}_1_trimmed.fastq.gz" \
            -O "/output/${SAMPLE}_2_trimmed.fastq.gz" \
            --cut_front \
            --cut_tail \
            --cut_window_size 4 \
            --cut_mean_quality 20 \
            --length_required 50 \
            --correction \
            --overrepresentation_analysis \
            --thread "$THREADS" \
            --html "/reports/${SAMPLE}_fastp_report.html" \
            --json "/reports/${SAMPLE}_fastp_report.json"
        
        if [ $? -eq 0 ]; then
            echo "   ✅ Fastp thành công cho $SAMPLE"
            echo "   📊 Report: ${REPORT_DIR}/${SAMPLE}_fastp_report.html"
        else
            echo "   ❌ Fastp thất bại cho $SAMPLE"
        fi
        echo "   ---"
    done
done

echo "🎉 Đã hoàn thành!"
echo "📂 Reads đã trim: $OUTPUT_DIR"
echo "📊 Báo cáo: $REPORT_DIR"

# Kiểm tra kết quả
echo "🔍 Kiểm tra kết quả:"
ls -la "$OUTPUT_DIR"
