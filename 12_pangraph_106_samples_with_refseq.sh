#!/bin/bash

# ==============================================
# CÀI ĐẶT THAM SỐ
# ==============================================
INPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/combined_contigs_and_refseq_v2"  # Thay bằng đường dẫn thực tế
OUTPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq"  # Thư mục output mới
OUTPUT_JSON="${OUTPUT_DIR}/pangenome_106_strains_with_refseq.json"
LOG_FILE="${OUTPUT_DIR}/pangenome_106_strains_with_refseq.log"
THREADS=8

# Tạo thư mục output nếu chưa tồn tại
mkdir -p "$OUTPUT_DIR"

# Bắt đầu tính thời gian
START_TIME=$(date +%s)

# ==============================================
# KIỂM TRA ĐẦU VÀO
# ==============================================
echo "[$(date +'%H:%M:%S')] 🚀 Bắt đầu chạy Pangraph"
echo "[$(date +'%H:%M:%S')] 📂 Kiểm tra thư mục input..."

# Kiểm tra thư mục input
if [ ! -d "$INPUT_DIR" ]; then
  echo "❌ Lỗi: Thư mục input $INPUT_DIR không tồn tại!"
  exit 1
fi

# Tìm các file fasta theo mẫu
FASTA_FILES=($(ls "${INPUT_DIR}"/*.filtered.complete.fasta 2>/dev/null))
TOTAL_GENOMES=${#FASTA_FILES[@]}

if [ "$TOTAL_GENOMES" -eq 0 ]; then
  echo "❌ Lỗi: Không tìm thấy file *.filtered.complete.fasta trong thư mục input!"
  exit 1
fi

echo "[$(date +'%H:%M:%S')] 🔍 Tìm thấy $TOTAL_GENOMES genome trong input"

# ==============================================
# CHẠY PANGRAPH VỚI THAM SỐ TỐI ƯU
# ==============================================
{
  echo "[$(date +'%H:%M:%S')] ⚙️ Tham số:"
  echo " - Aligner: minimap2 (-s 20)"
  echo " - Threads: $THREADS"
  echo " - Block min length: 300 (tăng từ 100 để giảm blocks)"
  echo " - Band width: 10 (tăng từ 5 để tránh band boundary warnings)"
  echo "[$(date +'%H:%M:%S')] 🏃 Đang chạy..."

  # Chạy Pangraph với tham số tối ưu
  pangraph build \
    -k minimap2 -s 20 -b 10 -a 300 -l 300 -x 300 -j "$THREADS" \
    -o "$OUTPUT_JSON" \
    "${FASTA_FILES[@]}" 2>&1
} | tee "$LOG_FILE"

# ==============================================
# KIỂM TRA KẾT QUẢ
# ==============================================
echo -e "\n[$(date +'%H:%M:%S')] ✅ Hoàn thành!"
echo "[$(date +'%H:%M:%S')] 📊 Kiểm tra kết quả:"

if [ -f "$OUTPUT_JSON" ]; then
  # Sử dụng jq -e để phát hiện lỗi cấu trúc JSON
  if BLOCKS=$(jq -e '.blocks | length' "$OUTPUT_JSON" 2>/dev/null); then
    STRAINS=$(jq -e '.strains | length' "$OUTPUT_JSON" 2>/dev/null || echo "Không đọc được strains")
    echo " - Số blocks: $BLOCKS"
    echo " - Số strains: $STRAINS"
    
    # Phân tích thêm về phân bố độ dài block
    echo -e "\n[$(date +'%H:%M:%S')] 📈 Phân bố độ dài block:"
    jq '.blocks[].alignment_length' "$OUTPUT_JSON" | sort -n | awk '
      NR==1 {min=$1}
      {sum+=$1; if($1>max)max=$1; count++}
      END {
        print " - Block ngắn nhất: " min " bp";
        print " - Block dài nhất: " max " bp";
        print " - Độ dài trung bình: " sum/count " bp";
      }'
  else
    echo "❌ Lỗi: File JSON không hợp lệ hoặc không có khóa 'blocks'"
    echo "Kiểm tra log file: $LOG_FILE"
  fi
else
  echo "❌ Lỗi: Không tìm thấy file output!"
fi

# Tính thời gian chạy chính xác
END_TIME=$(date +%s)
RUNTIME=$((END_TIME-START_TIME))
echo "[$(date +'%H:%M:%S')] ⏱ Thời gian chạy: $((RUNTIME/3600)) giờ $(( (RUNTIME%3600)/60 )) phút $((RUNTIME%60)) giây"
