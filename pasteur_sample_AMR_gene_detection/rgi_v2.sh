#!/bin/bash
# RGI ANALYSIS - FIXED DOCKER VOLUME MAPPING

# ================================
# CONFIGURATION
# ================================
INPUT_FASTA="/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur/combined_contigs/AB_Precipitated_Purple_20.final.complete.fasta"
OUTPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur/rgi_result"
SAMPLE_NAME="AB_Precipitated_Purple_20"

# ================================
# SCRIPT
# ================================

echo "============================================="
echo "RGI ANALYSIS - FINAL VERSION"
echo "Sample: $SAMPLE_NAME"
echo "============================================="
echo "📁 Input: $INPUT_FASTA"
echo "📁 Output: $OUTPUT_DIR"
echo "============================================="

# Kiểm tra input files
if [ ! -f "$INPUT_FASTA" ]; then
    echo "❌ ERROR: Input file not found: $INPUT_FASTA"
    echo "🔍 Đang tìm file trong thư mục: $(dirname "$INPUT_FASTA")"
    ls -la "$(dirname "$INPUT_FASTA")" 2>/dev/null || echo "Thư mục không tồn tại"
    exit 1
fi

# Tạo output directory
mkdir -p "$OUTPUT_DIR"

# FIX: Map cả input directory VÀ output directory vào container
INPUT_DIR=$(dirname "$INPUT_FASTA")
INPUT_FILE=$(basename "$INPUT_FASTA")

echo ""
echo "🚀 Running RGI with proper volume mapping..."

# SỬA LỖI: Map cả input và output directories
docker run --platform linux/amd64 \
  -v "$INPUT_DIR":/input \
  -v "$OUTPUT_DIR":/output \
  -w /input \
  quay.io/biocontainers/rgi:6.0.5--pyh05cac1d_0 \
  rgi main \
  -i "$INPUT_FILE" \
  -o "/output/rgi_strict" \
  -t contig \
  --clean

# Kiểm tra kết quả
if [ $? -eq 0 ]; then
    echo ""
    echo "✅ RGI analysis completed successfully!"
    echo ""
    echo "📊 OUTPUT FILES:"
    ls -la "$OUTPUT_DIR/rgi_strict"*
else
    echo "❌ RGI analysis failed!"
    exit 1
fi

# Summary statistics
if [ -f "$OUTPUT_DIR/rgi_strict.txt" ]; then
    echo ""
    echo "📈 SUMMARY STATISTICS:"
    TOTAL_HITS=$(tail -n +2 "$OUTPUT_DIR/rgi_strict.txt" | wc -l)
    echo "   Total AMR hits: $TOTAL_HITS"
    
    # Phân loại theo drug class
    echo ""
    echo "💊 DRUG CLASS DISTRIBUTION:"
    tail -n +2 "$OUTPUT_DIR/rgi_strict.txt" | cut -f9 | sort | uniq -c | sort -nr | head -10
    
    # Phân loại theo model (Perfect/Strict/Loose)
    echo ""
    echo "🏷️ MODEL DISTRIBUTION:"
    tail -n +2 "$OUTPUT_DIR/rgi_strict.txt" | cut -f5 | sort | uniq -c | sort -nr | head -5
fi

echo ""
echo "🎯 RGI ANALYSIS COMPLETED!"
echo "📂 Results saved in: $OUTPUT_DIR"