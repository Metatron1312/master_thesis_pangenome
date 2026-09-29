#!/bin/bash

# ================= CẤU HÌNH ĐƯỜNG DẪN =================
INPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/combined_contigs"  # Thư mục chứa file .fasta
OUTPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/busco_results"     # Thư mục kết quả
DATASET_DIR="/Users/shinra/busco_dataset"   # Thư mục dataset (tạo trước)
THREADS=8                               # Số thread
LINEAGE="acinetobacter_odb12"           # Dataset cho A. baumannii
# ======================================================

# 0. Tạo các thư mục cần thiết với quyền user hiện tại
mkdir -p "$DATASET_DIR" "$OUTPUT_DIR"
chmod -R 755 "$DATASET_DIR" "$OUTPUT_DIR"

# 1. Tải dataset nếu chưa có
echo -e "\n🔍 Kiểm tra dataset $LINEAGE..."
if [ ! -d "$DATASET_DIR/$LINEAGE" ]; then
    echo "📥 Dataset chưa có, đang tải về..."
    docker run --rm --platform linux/amd64 \
      -v "$DATASET_DIR":/busco_downloads \
      -u $(id -u):$(id -g) \
      ezlabgva/busco:v6.0.0_cv1 \
      busco --download "$LINEAGE" --download_path /busco_downloads --force
else
    echo "✅ Dataset đã tồn tại tại $DATASET_DIR/$LINEAGE"
fi

# 2. Tạo config file
CONFIG_FILE="$DATASET_DIR/busco_config.ini"
echo -e "[busco_run]\ndownload_path = /busco_downloads" > "$CONFIG_FILE"

# 3. Chạy BUSCO cho từng mẫu
for fasta_file in "$INPUT_DIR"/*.filtered.complete.fasta; do
    biosample=$(basename "$fasta_file" | cut -d'.' -f1)
    
    echo -e "\n🔬 Đang xử lý: $biosample"
    
    docker run --rm --platform linux/amd64 \
      -v "$INPUT_DIR":/input \
      -v "$OUTPUT_DIR":/output \
      -v "$DATASET_DIR":/busco_downloads \
      -v "$CONFIG_FILE":/config.ini \
      -u $(id -u):$(id -g) \
      -w /output \
      ezlabgva/busco:v6.0.0_cv1 \
      busco \
      -i "/input/$(basename "$fasta_file")" \
      -o "busco_${biosample}" \
      -m genome \
      -l "$LINEAGE" \
      -c "$THREADS" \
      --config /config.ini \
      --offline
    
    if [ $? -ne 0 ]; then
        echo "⚠️ Lỗi khi xử lý $biosample"
        echo "$(date) - Lỗi $biosample" >> "$OUTPUT_DIR/error.log"
    else
        echo "✅ Hoàn thành $biosample"
    fi
done

# 4. Tạo plot nếu có kết quả
if [ -n "$(ls -A "$OUTPUT_DIR"/busco_*/short_summary.*.txt 2>/dev/null)" ]; then
    echo -e "\n📊 Đang tạo plot tổng hợp..."
    docker run --rm --platform linux/amd64 \
      -v "$OUTPUT_DIR":/output \
      -u $(id -u):$(id -g) \
      -w /output \
      ezlabgva/busco:v6.0.0_cv1 \
      busco --plot /output
    
    mkdir -p "$OUTPUT_DIR/plots"
    mv "$OUTPUT_DIR"/busco_results*.png "$OUTPUT_DIR/plots/" 2>/dev/null
fi

# 5. Tổng hợp kết quả
echo -e "\n📝 Tổng hợp kết quả..."
grep -h "C:" "$OUTPUT_DIR"/busco_*/short_summary.*.txt > "$OUTPUT_DIR/summary_report.txt" 2>/dev/null || echo "Không có kết quả nào được tạo" > "$OUTPUT_DIR/summary_report.txt"

echo -e "\n✨ Hoàn thành! Kết quả:"
echo " - Dataset: $DATASET_DIR/$LINEAGE"
echo " - Kết quả BUSCO: $OUTPUT_DIR"
echo " - Plot: $OUTPUT_DIR/plots/"
echo " - Báo cáo: $OUTPUT_DIR/summary_report.txt"