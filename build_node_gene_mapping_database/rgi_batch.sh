#!/bin/bash
# BATCH RGI ANALYSIS - IMPROVED VERSION

# ================================
# CONFIGURATION
# ================================
INPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/combined_contigs_and_refseq_v2"
OUTPUT_BASE_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/Rgi_results"
THREADS_PER_SAMPLE=2

echo "============================================="
echo "BATCH RGI ANALYSIS - 107 SAMPLES"
echo "Started: $(date)"
echo "============================================="

# Kiểm tra Docker
if ! docker info > /dev/null 2>&1; then
    echo "❌ Docker is not running. Please start Docker Desktop first!"
    exit 1
fi

echo "✅ Docker is running"

mkdir -p "$OUTPUT_BASE_DIR"

# Counter
TOTAL_SAMPLES=0
SUCCESS_COUNT=0
FAIL_COUNT=0

echo "Scanning for assembly files..."
ASSEMBLY_FILES=("$INPUT_DIR"/*.filtered.complete.fasta)
TOTAL_FILES=${#ASSEMBLY_FILES[@]}

if [ $TOTAL_FILES -eq 0 ]; then
    echo "❌ No assembly files found in $INPUT_DIR"
    exit 1
fi

echo "Found $TOTAL_FILES assembly files"
echo ""

# ================================
# MAIN PROCESSING LOOP
# ================================
for assembly_file in "${ASSEMBLY_FILES[@]}"; do
    # Extract sample ID từ filename
    filename=$(basename "$assembly_file")
    SAMPLE_ID="${filename%%.filtered.complete.fasta}"
    
    ((TOTAL_SAMPLES++))
    
    echo "============================================="
    echo "[$TOTAL_SAMPLES/$TOTAL_FILES] Processing: $SAMPLE_ID"
    echo "Start time: $(date '+%H:%M:%S')"
    echo "============================================="
    
    # Tạo output directory cho sample
    SAMPLE_OUTPUT_DIR="$OUTPUT_BASE_DIR/$SAMPLE_ID"
    mkdir -p "$SAMPLE_OUTPUT_DIR"
    
    echo "📁 Input: $filename"
    echo "📁 Output: $SAMPLE_OUTPUT_DIR"
    
    # Chạy RGI với Docker - THÊM error handling
    echo "🚀 Running RGI analysis..."
    
    # Chạy RGI và capture output
    RGI_OUTPUT=$(docker run --rm --platform linux/amd64 \
      -v "$INPUT_DIR":/input \
      -v "$SAMPLE_OUTPUT_DIR":/output \
      -w /input \
      quay.io/biocontainers/rgi:6.0.5--pyh05cac1d_0 \
      rgi main \
      -i "$filename" \
      -o "/output/rgi_results" \
      -t contig \
      --clean \
      --low_quality \
      --num_threads "$THREADS_PER_SAMPLE" 2>&1)
    
    RGI_EXIT_CODE=$?
    
    # Kiểm tra kết quả
    if [ $RGI_EXIT_CODE -eq 0 ] && [ -f "$SAMPLE_OUTPUT_DIR/rgi_results.txt" ]; then
        # Thống kê kết quả
        TOTAL_HITS=$(tail -n +2 "$SAMPLE_OUTPUT_DIR/rgi_results.txt" 2>/dev/null | wc -l | tr -d ' ' || echo "0")
        
        echo "✅ RGI SUCCESS: $SAMPLE_ID"
        echo "   Total AMR hits: $TOTAL_HITS"
        
        # Hiển thị thông tin chi tiết
        if [ $TOTAL_HITS -gt 0 ]; then
            # Hiển thị drug classes với counts
            echo "   Drug class distribution:"
            tail -n +2 "$SAMPLE_OUTPUT_DIR/rgi_results.txt" | cut -f9 2>/dev/null | sort | uniq -c | sort -nr | head -5 | while read count class; do
                echo "     - $class: $count"
            done
            
            # Hiển thị AMR genes với drug class
            echo "   Top AMR genes:"
            tail -n +2 "$SAMPLE_OUTPUT_DIR/rgi_results.txt" | awk -F'\t' '{print $8 "\t" $9}' | sort | uniq | head -5 | while IFS=$'\t' read gene class; do
                echo "     - $gene ($class)"
            done
            
            # Hiển thị resistance mechanisms
            echo "   Resistance mechanisms:"
            tail -n +2 "$SAMPLE_OUTPUT_DIR/rgi_results.txt" | cut -f16 2>/dev/null | sort | uniq -c | sort -nr | head -3 | while read count mechanism; do
                echo "     - $mechanism: $count"
            done
        else
            echo "   ⚠️  No AMR hits found"
        fi
        
        ((SUCCESS_COUNT++))
        
        # Log warnings nếu có
        if echo "$RGI_OUTPUT" | grep -q "WARNING\|ERROR"; then
            echo "   ⚠️  RGI warnings detected (check detailed logs)"
        fi
        
    else
        echo "❌ RGI FAILED: $SAMPLE_ID (exit code: $RGI_EXIT_CODE)"
        echo "   Error output:"
        echo "$RGI_OUTPUT" | grep -i "error\|failed" | head -5
        ((FAIL_COUNT++))
    fi
    
    echo "Finish time: $(date '+%H:%M:%S')"
    echo "Elapsed time: ~1 minute"
    
    # Progress update
    echo ""
    echo "📈 Progress: $TOTAL_SAMPLES/$TOTAL_FILES ($SUCCESS_COUNT success, $FAIL_COUNT failed)"
    
    # Nghỉ giữa các samples để tránh quá tải
    if [ $TOTAL_SAMPLES -lt $TOTAL_FILES ]; then
        echo "💤 Taking a 10-second break..."
        for i in {10..1}; do
            echo -ne "   Resuming in $i seconds...\r"
            sleep 1
        done
        echo -ne "\n"
    fi
    echo ""
done

# ================================
# GENERATE SUMMARY REPORT
# ================================
echo ""
echo "============================================="
echo "GENERATING COMPREHENSIVE SUMMARY REPORT"
echo "============================================="

SUMMARY_FILE="$OUTPUT_BASE_DIR/rgi_batch_summary.tsv"
DETAILED_REPORT="$OUTPUT_BASE_DIR/rgi_comprehensive_report.txt"

# Tạo header cho file summary
echo -e "Sample_ID\tTotal_Hits\tDrug_Classes_Count\tTop_Drug_Class\tTop_Count\tAMR_Genes_Count\tResistance_Mechanisms" > "$SUMMARY_FILE"

{
    echo "RGI BATCH ANALYSIS - COMPREHENSIVE REPORT"
    echo "=========================================="
    echo "Generated: $(date)"
    echo "Total samples processed: $TOTAL_SAMPLES"
    echo "Successful analyses: $SUCCESS_COUNT"
    echo "Failed analyses: $FAIL_COUNT"
    echo "Success rate: $(echo "scale=1; $SUCCESS_COUNT * 100 / $TOTAL_SAMPLES" | bc)%"
    echo ""
    
    TOTAL_HITS_ALL=0
    SAMPLES_WITH_HITS=0
    
    echo "DETAILED SAMPLE RESULTS:"
    echo "========================"
    
    # Process each sample for summary
    for sample_dir in "$OUTPUT_BASE_DIR"/*/; do
        if [ -d "$sample_dir" ]; then
            SAMPLE_ID=$(basename "$sample_dir")
            RGI_FILE="$sample_dir/rgi_results.txt"
            
            if [ -f "$RGI_FILE" ] && [ -s "$RGI_FILE" ]; then
                TOTAL_HITS=$(tail -n +2 "$RGI_FILE" 2>/dev/null | wc -l | tr -d ' ' || echo "0")
                TOTAL_HITS_ALL=$((TOTAL_HITS_ALL + TOTAL_HITS))
                
                if [ $TOTAL_HITS -gt 0 ]; then
                    ((SAMPLES_WITH_HITS++))
                    
                    # Phân tích chi tiết
                    DRUG_CLASSES_COUNT=$(tail -n +2 "$RGI_FILE" | cut -f9 2>/dev/null | sort -u | wc -l)
                    TOP_CLASS_INFO=$(tail -n +2 "$RGI_FILE" | cut -f9 2>/dev/null | sort | uniq -c | sort -nr | head -1)
                    TOP_CLASS=$(echo "$TOP_CLASS_INFO" | awk '{$1=""; print $0}' | sed 's/^ //')
                    TOP_COUNT=$(echo "$TOP_CLASS_INFO" | awk '{print $1}')
                    
                    AMR_GENES_COUNT=$(tail -n +2 "$RGI_FILE" | cut -f8 2>/dev/null | sort -u | wc -l)
                    RESISTANCE_MECHANISMS=$(tail -n +2 "$RGI_FILE" | cut -f16 2>/dev/null | sort -u | tr '\n' ';' | sed 's/;$//')
                    
                    echo -e "$SAMPLE_ID\t$TOTAL_HITS\t$DRUG_CLASSES_COUNT\t$TOP_CLASS\t$TOP_COUNT\t$AMR_GENES_COUNT\t$RESISTANCE_MECHANISMS" >> "$SUMMARY_FILE"
                    
                    echo "🔬 $SAMPLE_ID"
                    echo "   Hits: $TOTAL_HITS | Genes: $AMR_GENES_COUNT | Drug Classes: $DRUG_CLASSES_COUNT"
                    echo "   Top: $TOP_CLASS ($TOP_COUNT)"
                    echo ""
                else
                    echo -e "$SAMPLE_ID\t0\t0\tNone\t0\t0\tNone" >> "$SUMMARY_FILE"
                    echo "🔬 $SAMPLE_ID: No AMR hits"
                    echo ""
                fi
            else
                echo -e "$SAMPLE_ID\t-1\t-1\tFailed\t-1\t-1\tFailed" >> "$SUMMARY_FILE"
                echo "❌ $SAMPLE_ID: Analysis failed or no results"
                echo ""
            fi
        fi
    done
    
    # Overall statistics
    echo "=========================================="
    echo "OVERALL STATISTICS"
    echo "=========================================="
    echo "Total samples with AMR hits: $SAMPLES_WITH_HITS"
    echo "Total AMR hits across all samples: $TOTAL_HITS_ALL"
    if [ $SUCCESS_COUNT -gt 0 ]; then
        AVG_HITS=$(echo "scale=2; $TOTAL_HITS_ALL / $SUCCESS_COUNT" | bc)
        echo "Average hits per sample: $AVG_HITS"
        HIT_RATE=$(echo "scale=1; $SAMPLES_WITH_HITS * 100 / $SUCCESS_COUNT" | bc)
        echo "Sample hit rate: $HIT_RATE%"
    fi
    echo ""
    
    # Top drug classes overall
    echo "TOP 10 DRUG CLASSES (Overall):"
    echo "------------------------------------------"
    find "$OUTPUT_BASE_DIR" -name "rgi_results.txt" -exec tail -n +2 {} \; 2>/dev/null | cut -f9 | sort | uniq -c | sort -nr | head -10
    echo ""
    
    # Top AMR genes overall
    echo "TOP 10 AMR GENES (Overall):"
    echo "------------------------------------------"
    find "$OUTPUT_BASE_DIR" -name "rgi_results.txt" -exec tail -n +2 {} \; 2>/dev/null | cut -f8 | sort | uniq -c | sort -nr | head -10
    echo ""
    
    # Resistance mechanisms
    echo "RESISTANCE MECHANISMS DISTRIBUTION:"
    echo "------------------------------------------"
    find "$OUTPUT_BASE_DIR" -name "rgi_results.txt" -exec tail -n +2 {} \; 2>/dev/null | cut -f16 | sort | uniq -c | sort -nr | head -10
    
} > "$DETAILED_REPORT"

# ================================
# FINAL SUMMARY
# ================================
echo ""
echo "============================================="
echo "🎉 BATCH RGI ANALYSIS COMPLETED!"
echo "Finished: $(date)"
echo "============================================="
echo "📊 FINAL RESULTS:"
echo "   Total samples processed: $TOTAL_SAMPLES"
echo "   ✅ Successful analyses: $SUCCESS_COUNT"
echo "   ❌ Failed analyses: $FAIL_COUNT"
echo "   🔬 Samples with AMR hits: $SAMPLES_WITH_HITS"
echo "   💊 Total AMR hits found: $TOTAL_HITS_ALL"
echo ""
echo "📁 OUTPUT FILES:"
echo "   📋 Summary TSV: $SUMMARY_FILE"
echo "   📄 Detailed report: $DETAILED_REPORT"
echo "   🔍 Individual results: $OUTPUT_BASE_DIR"
echo ""
echo "⏰ Total execution time: ~$(echo "scale=1; $TOTAL_SAMPLES * 1.2" | bc) minutes"
echo ""
echo "📍 Next step: Run Bakta-RGI integration when ready"