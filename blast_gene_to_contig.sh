#!/bin/bash
# =============================================================================
# SCRIPT: BLAST GENE LÊN CONTIG (FIXED PARSE)
# =============================================================================
# HƯỚNG DẪN: Chỉ cần sửa 3 biến dưới đây cho mỗi gene
# =============================================================================

# ======================== CẤU HÌNH (SỬA 3 DÒNG NÀY) ========================
# 1. Tên gene (ví dụ: adeB, adeB, tetB, OXA-23, ...)
GENE_NAME="adeB"

# 2. Đường dẫn đến file gene (có thể là .txt hoặc .fasta)
GENE_FILE="/Users/shinra/Reference_AMR_gene/adeB.txt"

# 3. Sample ID
SAMPLE_ID="SAMN48847114"
# ============================================================================

# ======================== CẤU HÌNH CỐ ĐỊNH ========================
BASE_DIR="/Users/shinra/gene_detection_v2/validate_database"
CONTIG_FILE="${BASE_DIR}/${SAMPLE_ID}/assembled_contigs/${SAMPLE_ID}/contigs.fa"
OUTPUT_DIR="${BASE_DIR}/${SAMPLE_ID}/${GENE_NAME}_blast"
mkdir -p ${OUTPUT_DIR}

echo "=========================================================================="
echo "BLAST ${GENE_NAME} LÊN CONTIG"
echo "=========================================================================="
echo "Sample: ${SAMPLE_ID}"
echo "Contig file: ${CONTIG_FILE}"
echo "${GENE_NAME} file: ${GENE_FILE}"
echo "Output: ${OUTPUT_DIR}"
echo "=========================================================================="

# ======================== KIỂM TRA FILE ========================
echo ""
echo "[Bước 1] Kiểm tra file..."

if [ ! -f "${CONTIG_FILE}" ]; then
    echo "❌ Không tìm thấy contig file: ${CONTIG_FILE}"
    exit 1
fi
echo "✅ Contig file: $(basename ${CONTIG_FILE})"

if [ ! -f "${GENE_FILE}" ]; then
    echo "❌ Không tìm thấy ${GENE_NAME} file: ${GENE_FILE}"
    exit 1
fi
echo "✅ ${GENE_NAME} file: $(basename ${GENE_FILE})"

# ======================== CHUYỂN ĐỔI SANG FASTA ========================
echo ""
echo "[Bước 2] Chuyển đổi file .txt sang FASTA format..."

GENE_FASTA="${OUTPUT_DIR}/${GENE_NAME}_reference.fa"

if grep -q "^>" "${GENE_FILE}" 2>/dev/null; then
    echo "   File đã có header FASTA, copy trực tiếp..."
    cp "${GENE_FILE}" "${GENE_FASTA}"
else
    echo "   File chưa có header FASTA, thêm header..."
    echo ">${GENE_NAME}_reference" > "${GENE_FASTA}"
    cat "${GENE_FILE}" >> "${GENE_FASTA}"
fi

GENE_LEN=$(grep -v "^>" "${GENE_FASTA}" | tr -d '\n' | wc -c | tr -d ' ')
echo "   Độ dài ${GENE_NAME}: ${GENE_LEN} bp"

# ======================== TẠO BLAST DATABASE ========================
echo ""
echo "[Bước 3] Tạo BLAST database từ contig..."

makeblastdb -in ${CONTIG_FILE} -dbtype nucl -out ${OUTPUT_DIR}/contig_db 2>/dev/null

if [ ! -f "${OUTPUT_DIR}/contig_db.nhr" ]; then
    echo "❌ Tạo BLAST database thất bại!"
    exit 1
fi
echo "✅ Đã tạo BLAST database"

# ======================== BLAST GENE ========================
echo ""
echo "[Bước 4] BLAST ${GENE_NAME} lên contig..."

BLAST_RESULT="${OUTPUT_DIR}/${GENE_NAME}_blast_results.txt"

blastn -query ${GENE_FASTA} \
       -db ${OUTPUT_DIR}/contig_db \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
       -out ${BLAST_RESULT} \
       -evalue 1e-5 \
       -num_threads 4 2>/dev/null

if [ ! -s ${BLAST_RESULT} ]; then
    echo ""
    echo "❌ KHÔNG TÌM THẤY ${GENE_NAME} TRONG CONTIG!"
    echo ""
    echo "💡 KẾT LUẬN:"
    echo "   - Sample KHÔNG có ${GENE_NAME} trong contig lắp ráp"
    exit 0
fi

# ======================== PHÂN TÍCH KẾT QUẢ ========================
echo ""
echo "[Bước 5] Phân tích kết quả BLAST..."

BLAST_COUNT=$(wc -l < ${BLAST_RESULT})
echo "✅ Tìm thấy ${BLAST_COUNT} hit(s)"

# Lấy hit tốt nhất và parse đúng các cột
BEST_HIT=$(head -1 ${BLAST_RESULT})

# Sử dụng awk để lấy từng cột
QSEQID=$(echo "$BEST_HIT" | awk -F'\t' '{print $1}')
SSEQID=$(echo "$BEST_HIT" | awk -F'\t' '{print $2}')
IDENTITY=$(echo "$BEST_HIT" | awk -F'\t' '{print $3}')
LENGTH=$(echo "$BEST_HIT" | awk -F'\t' '{print $4}')
QSTART=$(echo "$BEST_HIT" | awk -F'\t' '{print $7}')
QEND=$(echo "$BEST_HIT" | awk -F'\t' '{print $8}')
SSTART=$(echo "$BEST_HIT" | awk -F'\t' '{print $9}')
SEND=$(echo "$BEST_HIT" | awk -F'\t' '{print $10}')
EVALUE=$(echo "$BEST_HIT" | awk -F'\t' '{print $11}')

echo ""
echo "DEBUG: raw values"
echo "  qstart=$QSTART, qend=$QEND, sstart=$SSTART, send=$SEND"

# Kiểm tra nếu sstart > send (alignment trên mạch ngược)
if [ ${SSTART} -gt ${SEND} ]; then
    REAL_START=${SEND}
    REAL_END=${SSTART}
    STRAND="Reverse (-)"
else
    REAL_START=${SSTART}
    REAL_END=${SEND}
    STRAND="Forward (+)"
fi

# Tính coverage
COVERAGE=$(echo "scale=2; ${LENGTH}/${GENE_LEN}*100" | bc)

echo ""
echo "📊 HIT TỐT NHẤT:"
echo "---------------"
echo "  Query (${GENE_NAME}): ${QSEQID}"
echo "  Contig: ${SSEQID}"
echo "  Strand: ${STRAND}"
echo "  Position trên contig: ${REAL_START} - ${REAL_END}"
echo "  Position trên query: ${QSTART} - ${QEND}"
echo "  Alignment length: ${LENGTH} bp"
echo "  Identity: ${IDENTITY}%"
echo "  E-value: ${EVALUE}"
echo "  Coverage: ${COVERAGE}%"

# ======================== HIỂN THỊ TẤT CẢ HITS ========================
echo ""
echo "📋 TẤT CẢ HITS:"
echo "---------------"
echo -e "#\tContig\tIdentity\tLength\tCoverage\tPosition\tStrand"
echo "--------------------------------------------------------"

awk -F'\t' -v gene_len="${GENE_LEN}" '{
    qstart = $7; qend = $8; sstart = $9; send = $10
    if (sstart < send) {
        start = sstart; end = send; strand = "+"
    } else {
        start = send; end = sstart; strand = "-"
    }
    cov = ($4 / gene_len) * 100
    printf "%d\t%s\t%.1f%%\t%d\t%.1f%%\t%s-%s\t%s\n", NR, $2, $3, $4, cov, start, end, strand
}' ${BLAST_RESULT}

# ======================== KẾT LUẬN ========================
echo ""
echo "[Bước 6] Kết luận..."
echo "-------------------"

if (( $(echo "$IDENTITY > 95" | bc -l 2>/dev/null || echo 0) )) && (( $(echo "$COVERAGE > 90" | bc -l 2>/dev/null || echo 0) )); then
    echo "✅ ${GENE_NAME} ĐƯỢC PHÁT HIỆN với identity ${IDENTITY}% và coverage ${COVERAGE}%"
    echo "   → RGI strict mode NÊN phát hiện được gene này"
    echo "   → Nếu RGI không phát hiện, kiểm tra database CARD"
    STATUS="FOUND_HIGH"
elif (( $(echo "$IDENTITY > 80" | bc -l 2>/dev/null || echo 0) )); then
    echo "⚠️ ${GENE_NAME} được phát hiện với identity ${IDENTITY}% và coverage ${COVERAGE}%"
    echo "   → Cần dùng RGI loose mode để phát hiện"
    STATUS="FOUND_MEDIUM"
else
    echo "❌ ${GENE_NAME} KHÔNG được phát hiện (identity ${IDENTITY}%)"
    echo "   → Sample KHÔNG có ${GENE_NAME}"
    STATUS="NOT_FOUND"
fi

# ======================== LƯU BÁO CÁO ========================
echo ""
echo "[Bước 7] Lưu báo cáo..."

cat > ${OUTPUT_DIR}/${GENE_NAME}_report.txt << EOF
${GENE_NAME} BLAST REPORT
=================
Analysis Date: $(date)
Sample: ${SAMPLE_ID}
Contig file: ${CONTIG_FILE}
${GENE_NAME} file: ${GENE_FILE}
${GENE_NAME} length: ${GENE_LEN} bp

KẾT QUẢ BLAST:
--------------
Tổng số hits: ${BLAST_COUNT}

HIT TỐT NHẤT:
-------------
Query: ${QSEQID}
Contig: ${SSEQID}
Strand: ${STRAND}
Position trên contig: ${REAL_START} - ${REAL_END}
Position trên query: ${QSTART} - ${QEND}
Alignment length: ${LENGTH} bp
Identity: ${IDENTITY}%
Coverage: ${COVERAGE}%
E-value: ${EVALUE}

TẤT CẢ HITS:
------------
EOF

awk -F'\t' -v gene_len="${GENE_LEN}" '{
    qstart = $7; qend = $8; sstart = $9; send = $10
    if (sstart < send) {
        start = sstart; end = send; strand = "+"
    } else {
        start = send; end = sstart; strand = "-"
    }
    cov = ($4 / gene_len) * 100
    printf "%d. %s: %.1f%% identity, %.1f%% coverage, position %s-%s (%s)\n", NR, $2, $3, cov, start, end, strand
}' ${BLAST_RESULT} >> ${OUTPUT_DIR}/${GENE_NAME}_report.txt

cat >> ${OUTPUT_DIR}/${GENE_NAME}_report.txt << EOF

KẾT LUẬN:
---------
EOF

if [ "$STATUS" == "FOUND_HIGH" ]; then
    echo "✅ ${GENE_NAME} được phát hiện với identity ${IDENTITY}% và coverage ${COVERAGE}%" >> ${OUTPUT_DIR}/${GENE_NAME}_report.txt
    echo "   → RGI strict mode nên phát hiện được gene này" >> ${OUTPUT_DIR}/${GENE_NAME}_report.txt
elif [ "$STATUS" == "FOUND_MEDIUM" ]; then
    echo "⚠️ ${GENE_NAME} được phát hiện với identity ${IDENTITY}% và coverage ${COVERAGE}%" >> ${OUTPUT_DIR}/${GENE_NAME}_report.txt
    echo "   → Cần dùng RGI loose mode để phát hiện" >> ${OUTPUT_DIR}/${GENE_NAME}_report.txt
else
    echo "❌ ${GENE_NAME} KHÔNG được phát hiện (identity ${IDENTITY}%)" >> ${OUTPUT_DIR}/${GENE_NAME}_report.txt
    echo "   → Sample KHÔNG có ${GENE_NAME}" >> ${OUTPUT_DIR}/${GENE_NAME}_report.txt
fi

echo "✅ Đã lưu báo cáo: ${OUTPUT_DIR}/${GENE_NAME}_report.txt"

# ======================== TỔNG KẾT ========================
echo ""
echo "=========================================================================="
echo "✅ BLAST ${GENE_NAME} COMPLETE"
echo "=========================================================================="
echo ""
echo "📁 Output directory: ${OUTPUT_DIR}"
echo "   - ${GENE_NAME}_reference.fa        : File FASTA ${GENE_NAME}"
echo "   - ${GENE_NAME}_blast_results.txt   : Kết quả BLAST raw"
echo "   - ${GENE_NAME}_report.txt          : Báo cáo chi tiết"
echo ""
echo "📋 Kết quả tóm tắt:"
echo "   Contig: ${SSEQID}"
echo "   Position: ${REAL_START} - ${REAL_END}"
echo "   Strand: ${STRAND}"
echo "   Identity: ${IDENTITY}%"
echo "   Coverage: ${COVERAGE}%"
echo "   Status: ${STATUS}"
echo "=========================================================================="