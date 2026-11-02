#!/bin/bash

# Cấu hình đường dẫn
DB_PATH="/Users/shinra/Kraken2/db_standard_8GB_2025"
OUTPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur/kraken2_results"
ROOT_INPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur/platon_results"  # Thư mục gốc chứa các thư mục biosample

# Tạo thư mục output
mkdir -p "$OUTPUT_DIR"

classify_contigs() {
    local file="$1"
    local biosample_id=$(basename $(dirname "$file"))  # Lấy tên thư mục cha (biosample_id)
    local file_type=$(basename "$file" | rev | cut -d'.' -f2 | rev)  # Lấy chromosome hoặc plasmid từ tên file
    local output_prefix="$OUTPUT_DIR/${biosample_id}_${file_type}"

    echo "🔬 Processing: $file"

    # Chạy Docker với đường dẫn đơn giản hóa
    docker run --platform linux/amd64 --rm \
        -v "$DB_PATH:/db" \
        -v "$(dirname "$file"):/data" \
        -v "$OUTPUT_DIR:/output" \
        quay.io/biocontainers/kraken2:2.1.6--pl5321h077b44d_0 \
        kraken2 \
            --db /db \
            --threads 8 \
            --report "/output/${biosample_id}_${file_type}.report" \
            --output "/output/${biosample_id}_${file_type}.output" \
            "/data/$(basename "$file")"
}

# Tìm và xử lý file trong tất cả thư mục con
find "$ROOT_INPUT_DIR" -type f \( -name "*.with_metadata.chromosome.fasta" -o -name "*.with_metadata.plasmid.fasta" \) | while read file; do
    classify_contigs "$file"
done

echo "✅ Completed processing all samples"