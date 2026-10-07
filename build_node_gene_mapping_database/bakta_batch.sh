#!/bin/bash
# GENTLE BAKTA BATCH PROCESSING - FIXED WITH ENTRYPOINT

INPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/combined_contigs_and_refseq_v2"
OUTPUT_BASE_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/Bakta_results"
DB_DIR="/Users/shinra/bakta_db/db-light"

echo "============================================="
echo "GENTLE BAKTA BATCH PROCESSING"
echo "Started: $(date)"
echo "Running SLOWLY to prevent overheating"
echo "============================================="

# Kiểm tra Docker
if ! docker info > /dev/null 2>&1; then
    echo "❌ Docker is not running. Please start Docker Desktop first!"
    exit 1
fi

echo "✅ Docker is running"

mkdir -p "$OUTPUT_BASE_DIR"

count=0
success_count=0
total_files=$(find "$INPUT_DIR" -name "*.filtered.complete.fasta" | wc -l | tr -d ' ')

echo "Found $total_files assembly files"
echo "This will take a while... please be patient"
echo ""

for assembly_file in "$INPUT_DIR"/*.filtered.complete.fasta; do
    if [ ! -f "$assembly_file" ]; then
        echo "No assembly files found!"
        exit 1
    fi
    
    filename=$(basename "$assembly_file")
    biosample_id="${filename%%.filtered.complete.fasta}"
    output_dir="$OUTPUT_BASE_DIR/$biosample_id"
    
    ((count++))
    
    echo "============================================="
    echo "[$count/$total_files] Processing: $biosample_id"
    echo "Start time: $(date '+%H:%M:%S')"
    echo "============================================="
    
    # Tạo thư mục output
    mkdir -p "$output_dir"
    
    # SỬ DỤNG ĐÚNG CÚ PHÁP NHƯ SCRIPT CŨ - với --entrypoint bash -c
    docker run --rm --platform linux/amd64 \
      -v "$INPUT_DIR":/input \
      -v "$output_dir":/output \
      -v "/Users/shinra/bakta_db":/db \
      --entrypoint bash \
      oschwengers/bakta:latest \
      -c "bakta --db /db/db-light /input/$filename --output /output --prefix $biosample_id --threads 2 --compliant --force"
    
    # Kiểm tra kết quả
    if [ $? -eq 0 ] && [ -f "$output_dir/${biosample_id}.gff3" ]; then
        # Đếm features để xác nhận
        cds_count=$(grep -c "CDS" "$output_dir/${biosample_id}.gff3" 2>/dev/null || echo "0")
        amr_count=$(grep -c -i "antibiotic.resistance" "$output_dir/${biosample_id}.gff3" 2>/dev/null || echo "0")
        
        echo "✅ SUCCESS: $biosample_id"
        echo "   CDS features: $cds_count"
        echo "   AMR genes: $amr_count"
        ((success_count++))
    else
        echo "❌ FAILED: $biosample_id"
        echo "   Check the output in: $output_dir"
    fi
    
    echo "Finish time: $(date '+%H:%M:%S')"
    
    # NGHỈ GIẢI LAO GIỮA CÁC SAMPLES
    if [ $count -lt $total_files ]; then
        echo ""
        echo "💤 Taking a 30-second break before next sample..."
        echo "   (Letting the Mac cool down)"
        for i in {30..1}; do
            echo -ne "   Resuming in $i seconds...\r"
            sleep 1
        done
        echo -ne "\n🔄 Resuming processing..."
        echo ""
    fi
done

# Summary
echo ""
echo "============================================="
echo "GENTLE PROCESSING COMPLETED! 🎉"
echo "Finished: $(date)"
echo "============================================="
echo "📊 RESULTS:"
echo "   Total processed: $count"
echo "   Successful: $success_count"
echo "   Failed: $((count - success_count))"
if [ $count -gt 0 ]; then
    echo "   Success rate: $((success_count * 100 / count))%"
fi
echo ""
echo "📍 Output directory: $OUTPUT_BASE_DIR"
echo ""
echo "😊 All done! Your Mac should be happy and cool."