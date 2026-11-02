#!/bin/bash

# ================= CẤU HÌNH =================
CONTIGS_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur/combined_contigs"
OUTPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur/quast_results_combined"
THREADS=6
# ============================================

# Kiểm tra thư mục tồn tại
if [ ! -d "$CONTIGS_DIR" ]; then
    echo "❌ Lỗi: Thư mục không tồn tại - $CONTIGS_DIR"
    exit 1
fi

# Tạo output directory
mkdir -p "$OUTPUT_DIR"

# Tạo labels tự động từ tên file (lấy phần trước dấu .)
# Ví dụ: AB_Precipitated_Purple_20.final.complete.fasta → AB_Precipitated_Purple_20
labels=$(ls "$CONTIGS_DIR"/*.complete.fasta 2>/dev/null | \
         xargs -n1 basename | \
         sed 's/\..*//' | \
         tr '\n' ',' | \
         sed 's/,$//')

# Kiểm tra có file không
if [ -z "$labels" ]; then
    echo "❌ Lỗi: Không tìm thấy file .complete.fasta trong $CONTIGS_DIR"
    exit 1
fi

echo "🔍 Đang chạy QUAST cho các file:"
echo "$labels" | tr ',' '\n' | sed 's/^/   - /'
echo ""

# Chạy QUAST
quast.py "$CONTIGS_DIR"/*.complete.fasta \
    -o "$OUTPUT_DIR" \
    --labels "$labels" \
    --threads "$THREADS" \
    --gene-finding \
    --est-ref-size 4000000

# Kiểm tra kết quả
if [ $? -eq 0 ]; then
    echo ""
    echo "✅ QUAST đã hoàn thành thành công!"
    echo "📂 Kết quả được lưu tại: $(pwd)/$OUTPUT_DIR"
    echo ""
    echo "📊 Các file báo cáo chính:"
    echo "   - $(pwd)/$OUTPUT_DIR/report.html (Báo cáo HTML)"
    echo "   - $(pwd)/$OUTPUT_DIR/report.tsv (Báo cáo TSV)"
    echo "   - $(pwd)/$OUTPUT_DIR/report.txt (Báo cáo text)"
else
    echo "❌ QUAST thất bại. Vui lòng kiểm tra lỗi."
    exit 1
fi