#!/bin/bash
set -e

# ============================================================
# CẤU HÌNH
# ============================================================

# File alignment đầu vào (core genome)
INPUT_ALN="/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/iqtree/core_genome/core_genome_aln.fa"

# Thư mục output cho core_genome tree
OUTPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/iqtree/core_genome"

# Docker image
IQTREE_IMAGE="quay.io/biocontainers/iqtree:3.1.3--h8471819_0"

# Số threads (AUTO = tự động dùng tất cả CPU)
THREADS="AUTO"

# Số lần bootstrap
BOOTSTRAP=1000

# ============================================================
# KIỂM TRA INPUT
# ============================================================

echo "============================================================"
echo "KIỂM TRA INPUT"
echo "============================================================"

if [ ! -f "${INPUT_ALN}" ]; then
    echo "LỖI: Không tìm thấy file alignment: ${INPUT_ALN}"
    exit 1
fi

echo "File alignment: ${INPUT_ALN}"
echo "Kích thước: $(ls -lh "${INPUT_ALN}" | awk '{print $5}')"
echo "Số mẫu: $(grep -c ">" "${INPUT_ALN}")"
echo ""

# ============================================================
# TẠO THƯ MỤC OUTPUT
# ============================================================

mkdir -p "${OUTPUT_DIR}"
echo "Thư mục output: ${OUTPUT_DIR}"
echo ""

# ============================================================
# CHẠY IQ-TREE QUA DOCKER
# ============================================================

echo "============================================================"
echo "BẮT ĐẦU CHẠY IQ-TREE CHO CORE GENOME"
echo "============================================================"
echo "Thời gian bắt đầu: $(date)"
echo ""

docker run --rm \
    -v "$(dirname "${INPUT_ALN}"):/data/input:ro" \
    -v "${OUTPUT_DIR}:/data/output" \
    "${IQTREE_IMAGE}" \
    iqtree \
    -s "/data/input/$(basename "${INPUT_ALN}")" \
    -m MFP \
    -B "${BOOTSTRAP}" \
    -T "${THREADS}" \
    --prefix /data/output/core_genome_tree

echo ""
echo "============================================================"
echo "IQ-TREE ĐÃ CHẠY XONG"
echo "============================================================"
echo "Thời gian kết thúc: $(date)"
echo ""
echo "Kiểm tra output tại: ${OUTPUT_DIR}"
echo ""

# ============================================================
# LIỆT KÊ OUTPUT
# ============================================================

echo "Các file output:"
ls -lh "${OUTPUT_DIR}"
