#!/bin/bash
set -e

# ============================================================
# CẤU HÌNH
# ============================================================

# Thư mục chứa kết quả Prokka
PROKKA_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/Prokka_results"

# Thư mục output cho Panaroo
PANAROO_OUT="/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/panaroo_output"

# Số threads (giảm để tránh lỗi KeyError)
THREADS=4

# Docker image
PANAROO_IMAGE="quay.io/biocontainers/panaroo:1.8.0--pyhdfd78af_0"

# ============================================================
# TẠO THƯ MỤC OUTPUT
# ============================================================

mkdir -p "${PANAROO_OUT}/output"

# ============================================================
# TẠO DANH SÁCH GFF3 (đường dẫn container)
# ============================================================

GFF_LIST="${PANAROO_OUT}/gff_list.txt"
> "${GFF_LIST}"

for sample_dir in "${PROKKA_DIR}"/*/; do
    # Lấy tên thư mục, ví dụ: BJAB07104_prokka
    dir_name=$(basename "${sample_dir}")
    
    # Bỏ hậu tố "_prokka" để lấy sample name, ví dụ: BJAB07104
    sample_name="${dir_name%_prokka}"
    
    # File GFF nằm trong thư mục con, tên là <sample_name>.gff
    gff_file="${sample_dir}/${sample_name}.gff"
    
    if [ -f "${gff_file}" ]; then
        # Đường dẫn bên trong container
        echo "/data/gff/${dir_name}/${sample_name}.gff" >> "${GFF_LIST}"
        echo "Đã thêm: /data/gff/${dir_name}/${sample_name}.gff"
    else
        echo "CẢNH BÁO: Không tìm thấy ${gff_file}"
    fi
done

echo ""
echo "Tổng số file GFF3: $(wc -l < "${GFF_LIST}")"
echo ""

# ============================================================
# CHẠY PANAROO QUA DOCKER
# ============================================================

echo "Bắt đầu chạy Panaroo..."
echo "Thời gian bắt đầu: $(date)"
echo ""

docker run --rm \
    -v "${PROKKA_DIR}:/data/gff:ro" \
    -v "${PANAROO_OUT}/output:/data/output" \
    -v "${GFF_LIST}:/data/gff_list.txt:ro" \
    "${PANAROO_IMAGE}" \
    panaroo \
    -i /data/gff_list.txt \
    -o /data/output \
    --clean-mode strict \
    --remove-invalid-genes \
    --alignment core \
    --aligner mafft \
    --core_threshold 1.0 \
    --threads "${THREADS}" \
    --quiet

echo ""
echo "Panaroo đã chạy xong!"
echo "Thời gian kết thúc: $(date)"
echo ""
echo "Kiểm tra output tại: ${PANAROO_OUT}/output"