#!/bin/bash
# =============================================================================
# bwa_for_read_mapping_comparion.sh
# Map 1 mẫu lên linear reference bằng BWA-MEM.
# SAMPLE_ID được điền trực tiếp trong script.
# =============================================================================

set -euo pipefail

# --- SAMPLE ID (điền trực tiếp ở đây) ---
SAMPLE="SAMN48847114"
# ========================================

# --- Cấu hình ---
READS_BASE_DIR="/Users/shinra/gene_detection_v2/validate_database"
REFERENCE="/Users/shinra/Bioinformatics/WGS_AMR/data_public/combined_contigs_and_refseq/BJAB07104.filtered.complete.fasta"
OUTPUT_BASE_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/Mapping_and_variant_assessment/linear_reference/read_mapping_comparison"
THREADS=8
# =============================================================================

# --- Suy đường dẫn R1/R2 từ SAMPLE_ID ---
READS_DIR="${READS_BASE_DIR}/${SAMPLE}/fastp_trimmed"
R1="${READS_DIR}/${SAMPLE}_1_trimmed.fastq.gz"
R2="${READS_DIR}/${SAMPLE}_2_trimmed.fastq.gz"

echo "=== BWA-MEM MAPPING: ${SAMPLE} ==="
echo "Bắt đầu: $(date)"
echo ""

# --- Kiểm tra tools ---
command -v bwa      >/dev/null 2>&1 || { echo "❌ Không tìm thấy bwa"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "❌ Không tìm thấy samtools"; exit 1; }

# --- Kiểm tra input ---
[ -f "${R1}" ]        || { echo "❌ Không tìm thấy R1: ${R1}"; exit 1; }
[ -f "${R2}" ]        || { echo "❌ Không tìm thấy R2: ${R2}"; exit 1; }
[ -f "${REFERENCE}" ] || { echo "❌ Không tìm thấy reference: ${REFERENCE}"; exit 1; }

echo "📁 R1:  ${R1}"
echo "📁 R2:  ${R2}"
echo "📁 Ref: ${REFERENCE}"
echo ""

# --- Tạo BWA index nếu cần ---
if [ ! -f "${REFERENCE}.bwt" ]; then
    echo "🔄 Tạo BWA index cho reference..."
    bwa index "${REFERENCE}"
    echo "✅ Xong index"
fi

# --- Chuẩn bị thư mục ---
SAMPLE_OUT="${OUTPUT_BASE_DIR}/${SAMPLE}"
mkdir -p "${SAMPLE_OUT}"
TEMP_DIR="${SAMPLE_OUT}/tmp"
mkdir -p "${TEMP_DIR}"

SORTED_BAM="${SAMPLE_OUT}/${SAMPLE}.sorted.bam"
FLAGSTAT="${SAMPLE_OUT}/${SAMPLE}.flagstat.txt"

# --- Map ---
echo "🔄 BWA-MEM + sort..."
READ_GROUP="@RG\\tID:${SAMPLE}\\tSM:${SAMPLE}\\tPL:ILLUMINA\\tLB:lib1\\tPU:unit1"

bwa mem -t ${THREADS} -Y -K 100000000 -R "${READ_GROUP}" \
    "${REFERENCE}" "${R1}" "${R2}" | \
samtools sort -@ ${THREADS} -m 2G -T "${TEMP_DIR}" -o "${SORTED_BAM}" -

echo "🔄 Index BAM..."
samtools index "${SORTED_BAM}"

echo "📊 Flagstat..."
samtools flagstat "${SORTED_BAM}" > "${FLAGSTAT}"

# Hiển thị nhanh
TOTAL=$(grep "in total" "${FLAGSTAT}" | awk '{print $1}')
MAPPED=$(grep "mapped (" "${FLAGSTAT}" | awk '{print $1}')
if [ "${TOTAL}" -gt 0 ]; then
    RATE=$(awk "BEGIN {printf \"%.2f\", ${MAPPED}*100/${TOTAL}}")
    echo "🎯 Mapping rate (linear): ${RATE}%"
fi

# --- Dọn temp ---
rm -rf "${TEMP_DIR}"

echo ""
echo "🎉 HOÀN THÀNH ${SAMPLE}: $(date)"
echo "📂 BAM:      ${SORTED_BAM}"
echo "📂 Flagstat: ${FLAGSTAT}"
echo ""
echo "👉 Bước tiếp theo: chạy script 2 để thu thập stats."
