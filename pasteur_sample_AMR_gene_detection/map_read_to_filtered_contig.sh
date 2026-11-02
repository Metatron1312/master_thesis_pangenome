#!/bin/bash
# =============================================================================
# MAP READS VÀO CONTIG VÀ LỌC (SỬ DỤNG MINIMAP2)
# =============================================================================

# ======================== CẤU HÌNH ========================
BASE_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur"
SAMPLE_ID="AB_Precipitated_Purple_20"

# Đường dẫn
CONTIG_FILE="/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur/combined_contigs/AB_Precipitated_Purple_20.final.complete.fasta"
READ1="/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur/fastp_trimmed/AB_Precipitated_Purple_20_1_trimmed.fastq.gz"
READ2="/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur/fastp_trimmed/AB_Precipitated_Purple_20_2_trimmed.fastq.gz"
OUTPUT_DIR="${BASE_DIR}/mapped_reads_to_contig"

mkdir -p ${OUTPUT_DIR}

echo "=========================================================================="
echo "MAP READS VÀO CONTIG (MINIMAP2) VÀ LỌC"
echo "=========================================================================="
echo "Sample: ${SAMPLE_ID}"
echo "Contig: ${CONTIG_FILE}"
echo "Read1: ${READ1}"
echo "Read2: ${READ2}"
echo "Output: ${OUTPUT_DIR}"
echo "=========================================================================="

# ======================== KIỂM TRA FILE ========================
if [ ! -f "${CONTIG_FILE}" ]; then
    echo "❌ Không tìm thấy contig file: ${CONTIG_FILE}"
    exit 1
fi
echo "✅ Contig file: $(basename ${CONTIG_FILE})"

if [ ! -f "${READ1}" ]; then
    echo "❌ Không tìm thấy read1: ${READ1}"
    exit 1
fi
echo "✅ Read1: $(basename ${READ1})"

if [ ! -f "${READ2}" ]; then
    echo "❌ Không tìm thấy read2: ${READ2}"
    exit 1
fi
echo "✅ Read2: $(basename ${READ2})"

# ======================== MAP READS VÀO CONTIG (MINIMAP2) ========================
echo ""
echo "[Bước 1] Map reads vào contig với minimap2..."

# Tạo SAM file
SAM_FILE="${OUTPUT_DIR}/${SAMPLE_ID}.sam"
BAM_FILE="${OUTPUT_DIR}/${SAMPLE_ID}.bam"
SORTED_BAM="${OUTPUT_DIR}/${SAMPLE_ID}.sorted.bam"
FILTERED_BAM="${OUTPUT_DIR}/${SAMPLE_ID}.filtered.bam"

# Chạy minimap2 với tham số giống script của bạn
minimap2 -ax sr -t 8 --secondary=no -Y --eqx \
    ${CONTIG_FILE} \
    ${READ1} \
    ${READ2} \
    -o ${SAM_FILE} 2>/dev/null

if [ ! -s "${SAM_FILE}" ]; then
    echo "❌ Minimap2 mapping thất bại!"
    exit 1
fi
echo "✅ Minimap2 mapping thành công"

# ======================== CHUYỂN SANG BAM VÀ FILTER ========================
echo ""
echo "[Bước 2] Chuyển SAM -> BAM và filter..."

# SAM -> BAM
samtools view -bS ${SAM_FILE} > ${BAM_FILE} 2>/dev/null

# Sort BAM
samtools sort ${BAM_FILE} -o ${SORTED_BAM} 2>/dev/null

# Index BAM
samtools index ${SORTED_BAM} 2>/dev/null

# Filter: chỉ giữ reads có mapping quality >= 20 và không phải secondary alignment
# (Giống tham số min_mapping_quality: 20 trong script của bạn)
samtools view -b -q 20 -F 256 ${SORTED_BAM} > ${FILTERED_BAM} 2>/dev/null

# Index filtered BAM
samtools index ${FILTERED_BAM} 2>/dev/null

echo "✅ Filtered BAM created"

# ======================== LẤY READS ĐÃ FILTER ========================
echo ""
echo "[Bước 3] Lấy reads đã filter từ BAM..."

# Chuyển BAM thành FASTQ
samtools fastq -1 ${OUTPUT_DIR}/${SAMPLE_ID}_filtered_1.fastq.gz \
                -2 ${OUTPUT_DIR}/${SAMPLE_ID}_filtered_2.fastq.gz \
                ${FILTERED_BAM} 2>/dev/null

echo "✅ Filtered FASTQ files created"

# ======================== THỐNG KÊ ========================
echo ""
echo "[Bước 4] Thống kê..."

# Đếm số reads ban đầu
TOTAL_READS=$(zcat ${READ1} | echo $((`wc -l`/4)) 2>/dev/null)
echo "📊 Total reads ban đầu: ${TOTAL_READS} pairs"

# Đếm số reads sau filter
FILTERED_READS=$(zcat ${OUTPUT_DIR}/${SAMPLE_ID}_filtered_1.fastq.gz | echo $((`wc -l`/4)) 2>/dev/null)
echo "📊 Reads sau filter: ${FILTERED_READS} pairs"

# Tính tỷ lệ
RATIO=$(echo "scale=2; ${FILTERED_READS}/${TOTAL_READS}*100" | bc 2>/dev/null || echo "0")
echo "📊 Tỷ lệ giữ lại: ${RATIO}%"

# ======================== LƯU BÁO CÁO ========================
cat > ${OUTPUT_DIR}/filter_report.txt << EOF
READ FILTERING REPORT
=====================
Sample: ${SAMPLE_ID}
Date: $(date)
Contig file: ${CONTIG_FILE}
Mapping tool: minimap2
Parameters: -ax sr -t 8 --secondary=no -Y --eqx
Total reads: ${TOTAL_READS} pairs
Filtered reads: ${FILTERED_READS} pairs
Retention rate: ${RATIO}%
Filter criteria:
  - Mapping quality >= 20
  - Not secondary alignment
  - Read mapped to contig
Output files:
  - Filtered BAM: ${FILTERED_BAM}
  - Filtered FASTQ: ${OUTPUT_DIR}/${SAMPLE_ID}_filtered_1.fastq.gz
  - Filtered FASTQ: ${OUTPUT_DIR}/${SAMPLE_ID}_filtered_2.fastq.gz
EOF

echo ""
echo "📁 Filter report: ${OUTPUT_DIR}/filter_report.txt"

# ======================== CLEANUP ========================
echo ""
echo "[Bước 5] Dọn dẹp file tạm..."

rm -f ${SAM_FILE}
rm -f ${BAM_FILE}
rm -f ${SORTED_BAM}
rm -f ${SORTED_BAM}.bai

echo "✅ Dọn dẹp hoàn tất"

echo ""
echo "=========================================================================="
echo "✅ MAP READS VÀO CONTIG COMPLETE"
echo "=========================================================================="
echo ""
echo "📁 Output files:"
echo "   - ${FILTERED_BAM}"
echo "   - ${FILTERED_BAM}.bai"
echo "   - ${OUTPUT_DIR}/${SAMPLE_ID}_filtered_1.fastq.gz"
echo "   - ${OUTPUT_DIR}/${SAMPLE_ID}_filtered_2.fastq.gz"
echo "   - ${OUTPUT_DIR}/filter_report.txt"
echo "=========================================================================="
