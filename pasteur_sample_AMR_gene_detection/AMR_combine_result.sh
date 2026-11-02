#!/bin/bash

# --------------------------
# CẤU HÌNH ĐƯỜNG DẪN
# --------------------------
INPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur/AMR_Finder_results"          # Thư mục chứa file CSV đã có
COMBINED_RESULT="/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur/AMR_Finder_results/combined_amr_results.csv"  # File tổng hợp đầu ra

# --------------------------
# TẠO FILE TỔNG HỢP
# --------------------------
echo "📊 Bắt đầu tổng hợp kết quả AMRFinderPlus..."
echo "📁 Input directory: $INPUT_DIR"
echo "📄 Output file: $COMBINED_RESULT"
echo ""

# Kiểm tra thư mục input
if [ ! -d "$INPUT_DIR" ]; then
    echo "❌ Lỗi: Thư mục input không tồn tại: $INPUT_DIR"
    exit 1
fi

# Tạo header cho file tổng hợp
{
    echo "Type,Sample,Protein id,Contig id,Start,Stop,Strand,Element symbol,Element name,Scope,Type,Subtype,Class,Subclass,Method,Target length,Reference sequence length,% Coverage of reference,% Identity to reference,Alignment length,Closest reference accession,Closest reference name,HMM accession,HMM description"
} > "$COMBINED_RESULT"

# --------------------------
# GHÉP KẾT QUẢ TỪ CÁC FILE CSV
# --------------------------
TOTAL_FILES=0
TOTAL_ROWS=0

for csv_file in "$INPUT_DIR"/*.amrfinder.csv; do
    [ -f "$csv_file" ] || continue

    # Xác định loại (chromosome/plasmid) và tên mẫu từ tên file
    if [[ "$csv_file" == *"chromosome"* ]]; then
        type="chromosome"
    elif [[ "$csv_file" == *"plasmid"* ]]; then
        type="plasmid"
    else
        type="unknown"
    fi

    # Lấy tên sample (bỏ phần đuôi)
    sample_name=$(basename "$csv_file" | sed 's/\.amrfinder\.csv$//' | sed 's/\.contigs$//')
    
    TOTAL_FILES=$((TOTAL_FILES + 1))
    
    # Đếm số dòng trong file (không tính header)
    ROW_COUNT=$(tail -n +2 "$csv_file" 2>/dev/null | wc -l | tr -d ' ')
    TOTAL_ROWS=$((TOTAL_ROWS + ROW_COUNT))
    
    echo "   📄 Xử lý: $(basename "$csv_file") ($ROW_COUNT rows, type: $type, sample: $sample_name)"

    # Thêm dữ liệu vào file tổng hợp
    awk -v type="$type" -v sample="$sample_name" 'NR>1 {print type "," sample "," $0}' "$csv_file" >> "$COMBINED_RESULT"
done

# --------------------------
# HIỂN THỊ KẾT QUẢ
# --------------------------
echo ""
echo "============================================="
echo "✅ Đã tổng hợp hoàn tất!"
echo "============================================="
echo "📊 Thống kê:"
echo "   - Số file đã xử lý: $TOTAL_FILES"
echo "   - Tổng số dòng dữ liệu: $TOTAL_ROWS"
echo "   - File output: $COMBINED_RESULT"
echo "   - Kích thước file: $(du -h "$COMBINED_RESULT" 2>/dev/null | cut -f1 || echo "N/A")"
echo ""

# Hiển thị preview kết quả
if [ -f "$COMBINED_RESULT" ] && [ "$TOTAL_ROWS" -gt 0 ]; then
    echo "📋 Preview 5 dòng đầu tiên:"
    echo "----------------------------------------"
    head -n 6 "$COMBINED_RESULT" | column -t -s ','
    echo "----------------------------------------"
    echo ""
    
    # Thống kê theo loại
    echo "📈 Phân bố theo loại:"
    echo "   Chromosome: $(grep -c '^chromosome,' "$COMBINED_RESULT" 2>/dev/null || echo 0) genes"
    echo "   Plasmid: $(grep -c '^plasmid,' "$COMBINED_RESULT" 2>/dev/null || echo 0) genes"
    echo ""
    
    # Thống kê theo gene
    echo "🧬 Top 10 gene kháng tìm thấy:"
    tail -n +2 "$COMBINED_RESULT" 2>/dev/null | cut -d',' -f8 | sort | uniq -c | sort -rn | head -10 | while read count gene; do
        if [ -n "$gene" ] && [ "$gene" != "Element" ]; then
            echo "   - $gene: $count"
        fi
    done
else
    echo "⚠️ Không có dữ liệu nào được tổng hợp (có thể các file CSV trống)"
fi

echo "============================================="

# --------------------------
# KIỂM TRA CHI TIẾT
# --------------------------
echo ""
echo "🔍 Kiểm tra chi tiết từng file:"
for csv_file in "$INPUT_DIR"/*.amrfinder.csv; do
    [ -f "$csv_file" ] || continue
    ROW_COUNT=$(tail -n +2 "$csv_file" 2>/dev/null | wc -l | tr -d ' ')
    echo "   - $(basename "$csv_file"): $ROW_COUNT genes"
done

echo ""
echo "🎯 Hoàn thành!"