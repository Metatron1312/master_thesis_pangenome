#!/bin/bash

# Cấu hình đường dẫn
DB_PATH="/Users/shinra/Kraken2/db_standard_8GB_2025"
OUTPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/Kraken2_results_v2"
DATA_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/Platon_results_annotated"

# Tạo thư mục output
mkdir -p "$OUTPUT_DIR"

classify_contigs() {
    local file="$1"
    local biosample_id=$(basename "$file" | cut -d'.' -f1)
    local file_type=$(basename "$file" | cut -d'.' -f2)
    local output_prefix="$OUTPUT_DIR/${biosample_id}_${file_type}"

    echo "🔬 Processing: $file"

    # Chạy Docker với đường dẫn đơn giản hóa
    docker run --platform linux/amd64 --rm \
        -v "$DB_PATH:/db" \
        -v "$DATA_DIR:/data" \
        -v "$OUTPUT_DIR:/output" \
        quay.io/biocontainers/kraken2:2.1.6--pl5321h077b44d_0 \
        kraken2 \
            --db /db \
            --threads 8 \
            --report "/output/${biosample_id}_${file_type}.report" \
            --output "/output/${biosample_id}_${file_type}.output" \
            "/data/${biosample_id}/$(basename "$file")"
}

# Tìm và xử lý file
find "$DATA_DIR" -type f \( -name "*.chromosome.annotated.fasta" -o -name "*.plasmid.annotated.fasta" \) | while read file; do
    classify_contigs "$file"
done