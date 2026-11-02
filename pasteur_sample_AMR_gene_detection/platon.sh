#!/bin/bash

# Script chạy Platon trong Docker (phiên bản 1.0.0) cho các file genome dạng *_contigs.fa

# Cấu hình
INPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur/assembled/AB_Precipitated_Purple_20"
DB_PATH="/Users/shinra/platon/db"
OUTPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur/platon_results"
THREADS=8
VERBOSE=true
DOCKER_IMAGE="quay.io/biocontainers/platon:1.7--pyhdfd78af_0"

# Kiểm tra hệ thống
if [ ! -d "$INPUT_DIR" ]; then
    echo "Lỗi: Thư mục đầu vào không tồn tại!" >&2
    exit 1
fi

if ! command -v docker &> /dev/null; then
    echo "Lỗi: Docker không khả dụng!" >&2
    exit 1
fi

# Tạo thư mục output
mkdir -p "$OUTPUT_DIR"

# Xử lý từng file
find "$INPUT_DIR" -type f -name "*_contigs.fa" | while read -r file; do
    # Lấy biosample ID (bỏ _contigs.fa ở cuối)
    biosample_id=$(basename "$file" | sed 's/_contigs\.fa$//')
    output_subdir="${OUTPUT_DIR}/${biosample_id}"
    
    # Tạo thư mục kết quả
    mkdir -p "$output_subdir"
    
    # Lệnh Docker
    docker_cmd=(
        docker run --rm --platform linux/amd64
        -v "${DB_PATH}:/db"
        -v "$(dirname "$file"):/input"
        -v "${OUTPUT_DIR}:/output"
        "$DOCKER_IMAGE"
        platon --db /db --output "/output/${biosample_id}" --threads "$THREADS"
    )
    
    [ "$VERBOSE" = true ] && docker_cmd+=(--verbose)
    docker_cmd+=("/input/$(basename "$file")")
    
    # Thực thi
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Chạy Platon cho ${biosample_id}..."
    if "${docker_cmd[@]}"; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Hoàn thành ${biosample_id}"
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Lỗi xử lý ${biosample_id}" >&2
    fi
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Kết thúc quá trình!"