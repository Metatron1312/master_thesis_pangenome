#!/bin/bash

# ================= CẤU HÌNH ĐƯỜNG DẪN =================
INPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur/filtered_contigs/filtered_AMR"
OUTPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur/combined_contigs"
# ======================================================

# Tạo thư mục output nếu chưa tồn tại
mkdir -p "$OUTPUT_DIR"

echo "🔍 Đang xử lý ghép contig..."
echo "📂 Input: $INPUT_DIR"
echo "📂 Output: $OUTPUT_DIR"
echo ""

# Định nghĩa biosample (chỉ có 1 sample)
biosample="AB_Precipitated_Purple_20"

# Xác định các file cần ghép
chrom_file="$INPUT_DIR/${biosample}_contigs.with_metadata.chromosome.filtered.fasta"
plasmid_file="$INPUT_DIR/${biosample}_contigs.with_metadata.plasmid.filtered.fasta"

# File output
output_file="$OUTPUT_DIR/${biosample}.final.complete.fasta"

# Đếm số contig
chrom_count=0
plasmid_count=0

# Xóa file output cũ
rm -f "$output_file"

echo "🔬 Đang xử lý: $biosample"

# Xử lý chromosome
if [[ -f "$chrom_file" ]]; then
    echo "🟦 Đang thêm chromosome contigs từ: $(basename "$chrom_file")"
    chrom_count=$(grep -c '^>' "$chrom_file")
    cat "$chrom_file" >> "$output_file"
else
    echo "⚠️ Không tìm thấy file chromosome: $chrom_file"
fi

# Xử lý plasmid
if [[ -f "$plasmid_file" ]]; then
    echo "🟪 Đang thêm plasmid contigs từ: $(basename "$plasmid_file")"
    plasmid_count=$(grep -c '^>' "$plasmid_file")
    cat "$plasmid_file" >> "$output_file"
else
    echo "⚠️ Không tìm thấy file plasmid: $plasmid_file"
fi

# Thống kê
echo ""
echo "✅ Hoàn thành:"
echo "   - Chromosome contigs: $chrom_count"
echo "   - Plasmid contigs: $plasmid_count"
echo "   - Tổng contigs: $((chrom_count + plasmid_count))"
echo "   - File output: $(basename "$output_file")"
echo "   - Đường dẫn: $output_file"

# Kiểm tra file output
if [[ -f "$output_file" ]]; then
    echo ""
    echo "📊 Kiểm tra file kết quả:"
    echo "   - Kích thước: $(du -h "$output_file" 2>/dev/null | cut -f1)"
    echo "   - Tổng số dòng: $(wc -l < "$output_file")"
    echo "   - Tổng số contig: $(grep -c '^>' "$output_file")"
else
    echo ""
    echo "❌ Không tạo được file output!"
fi

echo -e "\n✨ Đã hoàn thành!"