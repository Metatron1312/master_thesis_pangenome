#!/bin/bash

# ================= CẤU HÌNH ĐƯỜNG DẪN =================
FILTERED_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/filter_contigs/filtered"
CHROMOSOME_CONTAM_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/filter_contigs/contaminated_chromosome"
PLASMID_CONTAM_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/filter_contigs/contaminated_plasmid"
AMR_TSV="/Users/shinra/Bioinformatics/WGS_AMR/data_public/AMR_Finder_results/final_amr_results_v6.tsv"
OUTPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/filter_contigs/filtered_AMR"
THREADS=$(sysctl -n hw.ncpu 2>/dev/null || echo 2) # Tương thích macOS
# ======================================================

# Kiểm tra các thư mục tồn tại
for dir in "$FILTERED_DIR" "$CHROMOSOME_CONTAM_DIR" "$PLASMID_CONTAM_DIR"; do
    if [ ! -d "$dir" ]; then
        echo "❌ Lỗi: Thư mục không tồn tại - $dir"
        exit 1
    fi
done

# Tạo thư mục output
mkdir -p "$OUTPUT_DIR"

# Bước 1: Đọc file TSV (phiên bản tương thích Bash <4)
CONTAM_WITH_AMR=""
while IFS=$'\t' read -r -a fields; do
    # Bỏ qua header và dòng không hợp lệ
    [[ "${fields[0]}" == "Biosample_ID" ]] && continue
    [[ -z "${fields[1]}" ]] && continue
    
    # Lưu dạng chuỗi thay vì mảng kết hợp
    CONTAM_WITH_AMR+="${fields[0]}.${fields[2]}.${fields[1]}"$'\n'
done < <(tr -d '\r' < "$AMR_TSV")

total_amr_contigs=$(echo "$CONTAM_WITH_AMR" | grep -c '^')
if [ "$total_amr_contigs" -eq 0 ]; then
    echo "⚠️ Cảnh báo: Không tìm thấy contig nào có gene kháng sinh trong file TSV!"
    exit 1
fi

echo "🔍 Đã tìm thấy $total_amr_contigs contig nhiễm có gene kháng sinh"

# Bước 2: Xử lý từng file filtered.fasta
processed_files=0
for filtered_file in "$FILTERED_DIR"/*.{chromosome,plasmid}.annotated.filtered.fasta; do
    [ -e "$filtered_file" ] || continue
    
    filename=$(basename "$filtered_file")
    biosample=$(echo "$filename" | cut -d'.' -f1)
    file_type=$(echo "$filename" | cut -d'.' -f2)
    
    echo -e "\n🔬 Đang xử lý: $biosample ($file_type)"

    # Xác định file contaminated
    if [[ "$file_type" == "plasmid" ]]; then
        contam_file="$PLASMID_CONTAM_DIR/${biosample}.plasmid.annotated.contaminated.fasta"
    else
        contam_file="$CHROMOSOME_CONTAM_DIR/${biosample}.chromosome.annotated.contaminated.fasta"
    fi
    
    if [[ ! -f "$contam_file" ]]; then
        echo "⚠️ Không tìm thấy file contaminated, copy nguyên filtered"
        cp "$filtered_file" "$OUTPUT_DIR/$filename"
        ((processed_files++))
        continue
    fi
    
    # Tạo danh sách contig nhiễm
    CONTAM_CONTIGS=$(grep '^>' "$contam_file" | cut -d' ' -f1 | tr -d '>')
    
    # Xử lý file filtered
    output_file="$OUTPUT_DIR/$filename"
    temp_file="${output_file}.tmp"
    rm -f "$temp_file"
    
    current_contig=""
    keep_contig=0
    removed_contigs=0
    kept_contigs=0
    
    while read -r line; do
        if [[ "$line" =~ ^\> ]]; then
            contig_id=$(echo "$line" | cut -d' ' -f1 | tr -d '>')
            
            # Kiểm tra contig nhiễm
            if echo "$CONTAM_CONTIGS" | grep -q "^${contig_id}$"; then
                # Kiểm tra có gene kháng không
                if echo "$CONTAM_WITH_AMR" | grep -q "^${biosample}.${file_type}.${contig_id}$"; then
                    keep_contig=1
                    ((kept_contigs++))
                    echo "💊 Giữ lại contig có gene kháng: $contig_id"
                else
                    keep_contig=0
                    ((removed_contigs++))
                    echo "🗑️ Loại bỏ contig không có gene kháng: $contig_id"
                fi
            else
                keep_contig=1
            fi
        fi
        
        if [[ "$keep_contig" -eq 1 ]]; then
            echo "$line" >> "$temp_file"
        fi
    done < "$filtered_file"
    
    mv "$temp_file" "$output_file"
    echo "✅ Đã xử lý:"
    echo "   - Contig giữ lại: $kept_contigs"
    echo "   - Contig loại bỏ: $removed_contigs"
    ((processed_files++))
done

echo -e "\n✨ Hoàn thành! Đã xử lý $processed_files files"
echo "📂 Kết quả lưu tại: $OUTPUT_DIR"