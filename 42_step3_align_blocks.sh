#!/bin/bash
set -e

INPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/iqtree/core_genome/core_block_seqs"
OUTPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/iqtree/core_genome/core_block_aligned"
CORE_BLOCKS="/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/iqtree/core_genome/core_blocks.txt"
MAFFT_IMAGE="staphb/mafft:latest"

mkdir -p "${OUTPUT_DIR}"

echo "============================================================"
echo "ALIGN TẤT CẢ CORE BLOCK BẰNG MAFFT (DOCKER)"
echo "============================================================"
echo "Thời gian bắt đầu: $(date)"
echo ""

# Copy danh sách core block vào thư mục input để container đọc được
cp "${CORE_BLOCKS}" "${INPUT_DIR}/core_blocks.txt"

# Chạy một container duy nhất, thực hiện vòng lặp bên trong
docker run --rm \
    -v "${INPUT_DIR}:/data/input:ro" \
    -v "${OUTPUT_DIR}:/data/output" \
    "${MAFFT_IMAGE}" \
    bash -c '
        count=0
        total=$(wc -l < /data/input/core_blocks.txt)
        while read block_id; do
            input_file="/data/input/block_${block_id}.fa"
            output_file="/data/output/block_${block_id}.aln.fa"
            
            if [ ! -f "${input_file}" ]; then
                continue
            fi
            
            mafft --auto --quiet "${input_file}" > "${output_file}"
            
            count=$((count + 1))
            if [ $((count % 50)) -eq 0 ]; then
                echo "  Đã align ${count}/${total} block"
            fi
        done < /data/input/core_blocks.txt
        
        echo "Hoàn tất: ${count} block"
    '

echo ""
echo "Thời gian kết thúc: $(date)"
echo ""
echo "Output dir: ${OUTPUT_DIR}"
echo "Số file output: $(ls "${OUTPUT_DIR}"/*.aln.fa 2>/dev/null | wc -l)"
