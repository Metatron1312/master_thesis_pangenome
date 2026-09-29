#!/bin/bash

# --------------------------
# CẤU HÌNH ĐƯỜNG DẪN
# --------------------------
INPUT_CHROMOSOME_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/filter_contigs/contaminated_chromosome"
INPUT_PLASMID_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/filter_contigs/contaminated_plasmid"
OUTPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/AMR_Finder_results"

# --------------------------
# THAM SỐ AMRFinderPlus
# --------------------------
AMR_IMAGE="ncbi/amr:latest"  # Sửa lại tên image theo đúng tên bạn có
THREADS=8

# Kiểm tra Docker image
check_docker_image() {
    if ! docker image inspect "$AMR_IMAGE" &> /dev/null; then
        echo "❌ Lỗi: Image Docker '$AMR_IMAGE' không tồn tại"
        echo "👉 Vui lòng cài đặt bằng lệnh:"
        echo "   docker pull ncbi/amr"
        exit 1
    fi
}

# Hàm chuyển TSV sang CSV
convert_to_csv() {
    local tsv_file="$1"
    [ ! -f "$tsv_file" ] && return 1
    
    local csv_file="${tsv_file%.tsv}.csv"
    awk 'BEGIN {FS="\t"; OFS=","} {
        gsub(/"/, "\"\"", $0)
        for (i=1; i<=NF; i++) {
            if ($i ~ /,/ || $i ~ /"/ || $i ~ /\s/) $i = "\"" $i "\""
        }
        print
    }' "$tsv_file" > "$csv_file" && rm "$tsv_file"
    
    echo "   → Đã chuyển thành CSV: $(basename "$csv_file")"
}

# Hàm xử lý từng file
process_file() {
    local input_path="$1"
    local file_type="$2"
    local sample_name=$(basename "$input_path" | cut -d '.' -f 1)
    local output_tsv="${OUTPUT_DIR}/${sample_name}.${file_type}.amrfinder.tsv"
    
    echo "🔍 Đang xử lý: $(basename "$input_path")..."
    
    # Mount volumes và chạy AMRFinderPlus
    if docker run --rm \
        -v "$(realpath "$INPUT_CHROMOSOME_DIR"):/chromosome" \
        -v "$(realpath "$INPUT_PLASMID_DIR"):/plasmid" \
        -v "$(realpath "$OUTPUT_DIR"):/output" \
        "$AMR_IMAGE" \
        amrfinder -n "/${file_type}/$(basename "$input_path")" \
        --plus \
        --threads "$THREADS" \
        -o "/output/$(basename "$output_tsv")"; then
        
        convert_to_csv "$output_tsv" || echo "   ⚠️ Không thể chuyển đổi sang CSV"
    else
        echo "   ❌ Xảy ra lỗi khi xử lý file"
    fi
    echo "--------------------------------------"
}

# --------------------------
# XỬ LÝ CHÍNH
# --------------------------
echo "🛠️ Bắt đầu quét AMR genes..."
check_docker_image
mkdir -p "$OUTPUT_DIR"

# Xử lý chromosome
for chrom_file in "$INPUT_CHROMOSOME_DIR"/*.chromosome.*.fasta; do
    [ -f "$chrom_file" ] || continue
    process_file "$chrom_file" "chromosome"
done

# Xử lý plasmid
for plasmid_file in "$INPUT_PLASMID_DIR"/*.plasmid.*.fasta; do
    [ -f "$plasmid_file" ] || continue
    process_file "$plasmid_file" "plasmid"
done

echo "🎉 Hoàn thành! Kết quả được lưu tại: $(realpath "$OUTPUT_DIR")"