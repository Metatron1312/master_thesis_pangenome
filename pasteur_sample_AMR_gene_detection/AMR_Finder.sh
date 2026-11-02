#!/bin/bash

# --------------------------
# CÀI ĐẶT AMRFINDERPLUS
# --------------------------
install_amrfinder() {
    echo "📦 Đang cài đặt AMRFinderPlus Docker image..."
    echo "   (Quá trình này có thể mất vài phút)"
    
    docker pull ncbi/amr:latest
    
    if [ $? -eq 0 ]; then
        echo "✅ Cài đặt AMRFinderPlus thành công!"
        return 0
    else
        echo "❌ Cài đặt AMRFinderPlus thất bại!"
        echo "👉 Thử pull từ nguồn khác: ncbi/amrfinderplus"
        docker pull ncbi/amrfinderplus:latest
        if [ $? -eq 0 ]; then
            echo "✅ Cài đặt AMRFinderPlus thành công (alternative image)!"
            AMR_IMAGE="ncbi/amrfinderplus:latest"
            return 0
        else
            echo "❌ Cài đặt AMRFinderPlus thất bại hoàn toàn!"
            return 1
        fi
    fi
}

# --------------------------
# CẤU HÌNH ĐƯỜNG DẪN
# --------------------------
INPUT_CHROMOSOME_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur/filtered_contigs/contaminated_chromosome"
INPUT_PLASMID_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur/filtered_contigs/contaminated_plasmid"
OUTPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur/AMR_Finder_results"

# --------------------------
# THAM SỐ AMRFinderPlus
# --------------------------
AMR_IMAGE="ncbi/amr:latest"  # Sẽ được cập nhật nếu cần
THREADS=8

# --------------------------
# KIỂM TRA VÀ CÀI ĐẶT
# --------------------------
check_and_install_amrfinder() {
    # Kiểm tra Docker image
    if ! docker image inspect "$AMR_IMAGE" &> /dev/null; then
        echo "⚠️ Image Docker '$AMR_IMAGE' không tồn tại"
        echo "👉 Bắt đầu cài đặt..."
        echo ""
        
        if install_amrfinder; then
            echo ""
            echo "✅ Đã cài đặt AMRFinderPlus thành công!"
            echo ""
            return 0
        else
            echo ""
            echo "❌ Không thể cài đặt AMRFinderPlus"
            echo "👉 Vui lòng thử cài đặt thủ công với lệnh:"
            echo "   docker pull ncbi/amr:latest"
            echo "   hoặc"
            echo "   docker pull ncbi/amrfinderplus:latest"
            return 1
        fi
    else
        echo "✅ AMRFinderPlus image đã tồn tại"
        return 0
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
echo "============================================="

# Tạo thư mục output
mkdir -p "$OUTPUT_DIR"

# Kiểm tra và cài đặt AMRFinderPlus
echo "📌 Bước 1: Kiểm tra AMRFinderPlus..."
if ! check_and_install_amrfinder; then
    exit 1
fi

echo ""
echo "📌 Bước 2: Kiểm tra thư mục input..."

# Kiểm tra thư mục chromosome
if [ ! -d "$INPUT_CHROMOSOME_DIR" ]; then
    echo "⚠️ Thư mục chromosome không tồn tại: $INPUT_CHROMOSOME_DIR"
    echo "   (Bỏ qua bước xử lý chromosome)"
else
    CHROM_FILES=$(find "$INPUT_CHROMOSOME_DIR" -name "*.chromosome.*.fasta" 2>/dev/null | wc -l)
    echo "   Chromosome files: $CHROM_FILES files"
fi

# Kiểm tra thư mục plasmid
if [ ! -d "$INPUT_PLASMID_DIR" ]; then
    echo "⚠️ Thư mục plasmid không tồn tại: $INPUT_PLASMID_DIR"
    echo "   (Bỏ qua bước xử lý plasmid)"
else
    PLASMID_FILES=$(find "$INPUT_PLASMID_DIR" -name "*.plasmid.*.fasta" 2>/dev/null | wc -l)
    echo "   Plasmid files: $PLASMID_FILES files"
fi

echo ""
echo "📌 Bước 3: Xử lý AMR detection..."

# Xử lý chromosome
if [ -d "$INPUT_CHROMOSOME_DIR" ]; then
    echo ""
    echo "🔬 Xử lý CHROMOSOME files:"
    echo "--------------------------------------"
    for chrom_file in "$INPUT_CHROMOSOME_DIR"/*.chromosome.*.fasta; do
        [ -f "$chrom_file" ] || continue
        process_file "$chrom_file" "chromosome"
    done
fi

# Xử lý plasmid
if [ -d "$INPUT_PLASMID_DIR" ]; then
    echo ""
    echo "🧬 Xử lý PLASMID files:"
    echo "--------------------------------------"
    for plasmid_file in "$INPUT_PLASMID_DIR"/*.plasmid.*.fasta; do
        [ -f "$plasmid_file" ] || continue
        process_file "$plasmid_file" "plasmid"
    done
fi

# Tổng kết
echo ""
echo "============================================="
echo "🎉 Hoàn thành! Kết quả được lưu tại:"
echo "   $(realpath "$OUTPUT_DIR")"

# Hiển thị kết quả
if [ -d "$OUTPUT_DIR" ]; then
    echo ""
    echo "📊 Kết quả tìm thấy:"
    ls -la "$OUTPUT_DIR" | grep -E "\.(tsv|csv)$" | awk '{print "   - " $9 " (" $5 " bytes)"}'
    
    # Tổng số gene kháng tìm thấy
    echo ""
    echo "📈 Thống kê nhanh:"
    for tsv_file in "$OUTPUT_DIR"/*.tsv; do
        [ -f "$tsv_file" ] || continue
        # Bỏ header (dòng 1) và đếm số dòng
        COUNT=$(tail -n +2 "$tsv_file" 2>/dev/null | wc -l | tr -d ' ')
        if [ "$COUNT" -gt 0 ]; then
            echo "   - $(basename "$tsv_file"): $COUNT AMR genes detected"
        fi
    done
fi

echo "============================================="