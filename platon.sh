#!/bin/bash

# Script chạy Platon trong Docker (phiên bản 1.0.0) cho các file genome

# Cấu hình
INPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/genomes"
DB_PATH="/Users/shinra/platon/db"
OUTPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/Platon_results"
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
find "$INPUT_DIR" -type f \( -name "*.fa.gz" -o -name "*.fna.gz" \) | while read -r file; do
    # Lấy biosample ID (bỏ _genomic.fna.gz, .contigs.fa.gz)
    biosample_id=$(basename "$file" | sed -E 's/\.(fa|fna)\.gz$//' | sed -E 's/(_genomic|\.contigs)//')
    temp_fasta="${OUTPUT_DIR}/${biosample_id}.temp.fasta"
    
    # Giải nén
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Giải nén $file..."
    if ! gunzip -ck "$file" > "$temp_fasta"; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Lỗi giải nén!" >&2
        continue
    fi
    
    # Tạo thư mục kết quả
    mkdir -p "${OUTPUT_DIR}/${biosample_id}"
    
    # Lệnh Docker (phiên bản 1.0.0 không có --mode)
    docker_cmd=(
        docker run --rm --platform linux/amd64
        -v "${DB_PATH}:/db"
        -v "${OUTPUT_DIR}:/data"
        "$DOCKER_IMAGE"
        platon --db /db --output "/data/${biosample_id}" --threads "$THREADS"
    )
    
    [ "$VERBOSE" = true ] && docker_cmd+=(--verbose)
    docker_cmd+=("/data/${biosample_id}.temp.fasta")
    
    # Thực thi
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Chạy Platon cho ${biosample_id}..."
    if "${docker_cmd[@]}"; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Hoàn thành ${biosample_id}"
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Lỗi xử lý ${biosample_id}" >&2
    fi
    
    # Dọn dẹp
    rm -f "$temp_fasta"
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Kết thúc quá trình!"
