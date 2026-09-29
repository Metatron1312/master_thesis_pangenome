#!/bin/bash

# ================= CẤU HÌNH ĐƯỜNG DẪN =================
INPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/filter_contigs/filtered_AMR"  # Thư mục chứa cả chromosome và plasmid
OUTPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/combined_contigs"            # Thư mục kết quả
# ======================================================

# Tạo thư mục output nếu chưa tồn tại
mkdir -p "$OUTPUT_DIR"

# Lấy danh sách các biosample duy nhất
biosamples=$(ls "$INPUT_DIR" | grep -o '^[^.]*' | sort -u)

echo "🔍 Đang xử lý ghép contig cho $(echo "$biosamples" | wc -l) biosamples..."

# Duyệt qua từng biosample
for biosample in $biosamples; do
    echo -e "\n🔬 Đang xử lý: $biosample"
    
    # Xác định các file cần ghép
    chrom_file=$(ls "$INPUT_DIR"/"$biosample".chromosome.*.fasta 2>/dev/null | head -n1)
    plasmid_file=$(ls "$INPUT_DIR"/"$biosample".plasmid.*.fasta 2>/dev/null | head -n1)
    
    # Tạo file output
    output_file="$OUTPUT_DIR/${biosample}.filtered.complete.fasta"
    
    # Đếm số contig
    chrom_count=0
    plasmid_count=0
    
    # Xử lý chromosome
    if [[ -f "$chrom_file" ]]; then
        echo "🟦 Đang thêm chromosome contigs từ: $(basename "$chrom_file")"
        chrom_count=$(grep -c '^>' "$chrom_file")
        cat "$chrom_file" >> "$output_file"
    else
        echo "⚠️ Không tìm thấy file chromosome cho $biosample"
    fi
    
    # Xử lý plasmid
    if [[ -f "$plasmid_file" ]]; then
        echo "🟪 Đang thêm plasmid contigs từ: $(basename "$plasmid_file")"
        plasmid_count=$(grep -c '^>' "$plasmid_file")
        cat "$plasmid_file" >> "$output_file"
    else
        echo "⚠️ Không tìm thấy file plasmid cho $biosample"
    fi
    
    # Thống kê
    echo "✅ Hoàn thành:"
    echo "   - Chromosome contigs: $chrom_count"
    echo "   - Plasmid contigs: $plasmid_count"
    echo "   - Tổng contigs: $((chrom_count + plasmid_count))"
    echo "   - File output: $(basename "$output_file")"
done

echo -e "\n✨ Tất cả đã hoàn thành!"
echo "📂 Kết quả được lưu tại: $OUTPUT_DIR"
echo "📊 Tổng số file đã tạo: $(ls "$OUTPUT_DIR"/*.complete.fasta 2>/dev/null | wc -l)"
