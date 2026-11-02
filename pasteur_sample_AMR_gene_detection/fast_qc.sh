#!/bin/bash

# Script: fast_qc.sh
# Phiên bản tích hợp sẵn đường dẫn input/output

# ================= CẤU HÌNH ĐƯỜNG DẪN =================
INPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/Assess_SNP_SV/data/raw/SAMEA8218231"    # Thay bằng thư mục chứa file reads của bạn
OUTPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/Assess_SNP_SV/data/QC_checked"   # Thay bằng thư mục bạn muốn lưu kết quả
# =======================================================

# Kiểm tra FastQC
if ! command -v fastqc &> /dev/null; then
    echo "Lỗi: Cần cài đặt FastQC trước!"
    echo "Cài đặt: sudo apt-get install fastqc hoặc conda install -c bioconda fastqc"
    exit 1
fi

# Tạo thư mục output nếu chưa tồn tại
mkdir -p "$OUTPUT_DIR"

# Số luồng CPU (cho macOS và Linux)
THREADS=$(sysctl -n hw.ncpu 2>/dev/null || echo 2)

# Tìm file R1
R1_FILES=("$INPUT_DIR"/*_R1.fastq "$INPUT_DIR"/*_R1.fastq.gz)

if [ ${#R1_FILES[@]} -eq 0 ]; then
    echo "Không tìm thấy file R1 nào trong $INPUT_DIR!"
    echo "Đảm bảo file có tên dạng *_R1.fastq hoặc *_R1.fastq.gz"
    exit 1
fi

echo "🔍 Đã tìm thấy ${#R1_FILES[@]} cặp reads trong $INPUT_DIR:"

for R1_FILE in "${R1_FILES[@]}"; do
    # Xác định file R2 (thay R1 bằng R2)
    R2_FILE="${R1_FILE/_R1./_R2.}"
    R2_FILE="${R2_FILE/_R1_/_R2_}"
    
    if [ ! -f "$R2_FILE" ]; then
        echo "⚠️ Không tìm thấy file R2 tương ứng cho: $(basename "$R1_FILE")"
        continue
    fi
    
    echo "▶️ Đang xử lý cặp:"
    echo "   - R1: $(basename "$R1_FILE")"
    echo "   - R2: $(basename "$R2_FILE")"
    
    fastqc "$R1_FILE" "$R2_FILE" -o "$OUTPUT_DIR" -t $THREADS
    
    if [ $? -eq 0 ]; then
        echo "   ✅ Thành công"
    else
        echo "   ❌ Thất bại"
    fi
done

echo "✅ Đã hoàn thành FastQC!"
echo "📂 Kết quả lưu tại: $OUTPUT_DIR"

# Tự động chạy MultiQC nếu đã cài đặt
if command -v multiqc &> /dev/null; then
    echo "🔄 Đang tạo báo cáo tổng hợp..."
    multiqc "$OUTPUT_DIR" -o "${OUTPUT_DIR}_multiqc"
    echo "📊 Báo cáo đa hợp: ${OUTPUT_DIR}_multiqc/multiqc_report.html"
else
    echo "ℹ️ Cài MultiQC để tổng hợp kết quả: conda install -c bioconda multiqc"
fi
