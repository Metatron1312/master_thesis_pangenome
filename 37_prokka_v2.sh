#!/bin/bash

# =============================================
# PROKKA SCALE UP - ALL SAMPLES WITH CARD DATABASE
# =============================================

# Configurations
INPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/combined_contigs_and_refseq_v2"
OUTPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/Prokka_results"
CARD_DB="/Users/shinra/CARD_db/protein_fasta_protein_homolog_model.fasta"
THREADS=8
DOCKER_IMAGE="staphb/prokka:1.14.6"

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "============================================="
echo "PROKKA SCALE UP - ALL SAMPLES WITH CARD DATABASE"
echo "============================================="
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "CARD database: $(basename "$CARD_DB")"
echo "Threads: $THREADS"
echo "Start time: $(date)"
echo "============================================="

# Check if CARD database exists
if [ ! -f "$CARD_DB" ]; then
    echo "❌ ERROR: CARD database not found: $CARD_DB"
    echo "Available CARD files:"
    ls -la /Users/shinra/CARD_db/*.fasta
    exit 1
fi

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "❌ ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

# Get all assembly files
ASSEMBLY_FILES=("$INPUT_DIR"/*.filtered.complete.fasta)
TOTAL_SAMPLES=${#ASSEMBLY_FILES[@]}

if [ $TOTAL_SAMPLES -eq 0 ]; then
    echo "❌ ERROR: No assembly files found in $INPUT_DIR"
    exit 1
fi

echo "📁 Found $TOTAL_SAMPLES assembly files"
echo ""

# Counter for tracking progress
SUCCESS_COUNT=0
FAILED_COUNT=0
FAILED_SAMPLES=()

# Process each assembly file
for ((i=0; i<${#ASSEMBLY_FILES[@]}; i++)); do
    ASSEMBLY="${ASSEMBLY_FILES[i]}"
    FILENAME=$(basename "$ASSEMBLY")
    SAMPLE_NAME=$(echo "$FILENAME" | sed 's/\.filtered\.complete\.fasta//g')
    
    echo "============================================="
    echo "🔬 [$((i+1))/$TOTAL_SAMPLES] Processing: $SAMPLE_NAME"
    echo "============================================="
    
    # Create sample-specific output directory
    SAMPLE_OUTDIR="$OUTPUT_DIR/${SAMPLE_NAME}_prokka"
    
    # Skip if already processed successfully
    if [ -f "$SAMPLE_OUTDIR/${SAMPLE_NAME}.gff" ]; then
        echo "⚠️  Already processed, skipping: $SAMPLE_NAME"
        ((SUCCESS_COUNT++))
        echo ""
        continue
    fi
    
    # Remove if partial results exist
    if [ -d "$SAMPLE_OUTDIR" ]; then
        echo "⚠️  Partial results found, removing: $SAMPLE_OUTDIR"
        rm -rf "$SAMPLE_OUTDIR"
    fi
    
    # Run Prokka with Docker
    echo "🚀 Running Prokka with CARD database..."
    echo "⏰ Start: $(date)"
    
    docker run --rm --platform linux/amd64 \
        -v "$INPUT_DIR":/input \
        -v "$OUTPUT_DIR":/output \
        -v "/Users/shinra/CARD_db":/database \
        "$DOCKER_IMAGE" \
        prokka \
        --outdir "/output/${SAMPLE_NAME}_prokka" \
        --prefix "$SAMPLE_NAME" \
        --cpus "$THREADS" \
        --locustag "$SAMPLE_NAME" \
        --kingdom Bacteria \
        --proteins "/database/$(basename "$CARD_DB")" \
        --compliant \
        "/input/$(basename "$ASSEMBLY")"
    
    # Check if Prokka completed successfully
    if [ -f "$SAMPLE_OUTDIR/${SAMPLE_NAME}.gff" ]; then
        echo "✅ SUCCESS: $SAMPLE_NAME"
        
        # Count resistance genes
        RESISTANCE_COUNT=$(grep -c -i "resistance\|antibiotic" "$SAMPLE_OUTDIR/${SAMPLE_NAME}.gff" 2>/dev/null || echo "0")
        echo "🦠 Resistance genes: $RESISTANCE_COUNT"
        
        ((SUCCESS_COUNT++))
    else
        echo "❌ FAILED: $SAMPLE_NAME"
        FAILED_SAMPLES+=("$SAMPLE_NAME")
        ((FAILED_COUNT++))
        
        # Show error if exists
        if [ -f "$SAMPLE_OUTDIR/${SAMPLE_NAME}.err" ]; then
            echo "📋 Error log (first 5 lines):"
            head -5 "$SAMPLE_OUTDIR/${SAMPLE_NAME}.err"
        fi
    fi
    
    echo "⏰ End: $(date)"
    echo "📊 Progress: $SUCCESS_COUNT/$TOTAL_SAMPLES completed"
    echo ""
    
    # Optional: Add delay between samples to avoid overloading system
    # sleep 10
done

# Generate comprehensive summary report
echo "============================================="
echo "📈 FINAL SUMMARY REPORT"
echo "============================================="
echo "Start time: $(date -r "$OUTPUT_DIR" 2>/dev/null || echo "N/A")"
echo "End time: $(date)"
echo "Total samples: $TOTAL_SAMPLES"
echo "✅ Successfully annotated: $SUCCESS_COUNT"
echo "❌ Failed: $FAILED_COUNT"
echo ""

# Show resistance gene summary
if [ $SUCCESS_COUNT -gt 0 ]; then
    echo "🦠 RESISTANCE GENE SUMMARY"
    echo "=========================="
    
    # Collect resistance gene counts
    declare -A resistance_totals
    TOTAL_RESISTANCE=0
    
    for SAMPLE_DIR in "$OUTPUT_DIR"/*_prokka; do
        if [ -f "$SAMPLE_DIR"/*.gff ]; then
            SAMPLE=$(basename "$SAMPLE_DIR" _prokka)
            COUNT=$(grep -c -i "resistance\|antibiotic" "$SAMPLE_DIR"/*.gff 2>/dev/null || echo "0")
            TOTAL_RESISTANCE=$((TOTAL_RESISTANCE + COUNT))
            
            # Store for top samples
            if [ "$COUNT" -gt 0 ]; then
                resistance_totals["$SAMPLE"]=$COUNT
            fi
        fi
    done
    
    echo "Total resistance genes found: $TOTAL_RESISTANCE"
    echo "Average per sample: $((TOTAL_RESISTANCE / SUCCESS_COUNT))"
    echo ""
    
    # Show top samples with most resistance genes
    if [ ${#resistance_totals[@]} -gt 0 ]; then
        echo "🏆 TOP SAMPLES BY RESISTANCE GENES:"
        for SAMPLE in "${!resistance_totals[@]}"; do
            echo "  $SAMPLE: ${resistance_totals[$SAMPLE]}"
        done | sort -nr -k2 | head -10
    fi
fi

# Show failed samples
if [ $FAILED_COUNT -gt 0 ]; then
    echo ""
    echo "❌ FAILED SAMPLES:"
    for SAMPLE in "${FAILED_SAMPLES[@]}"; do
        echo "  - $SAMPLE"
    done
fi

# Next steps
echo ""
echo "============================================="
echo "🎉 PROKKA PIPELINE COMPLETED!"
echo "============================================="
echo "📁 Output directory: $OUTPUT_DIR"
echo ""
echo "🚀 NEXT STEP - Annotate your pangenome graph:"
echo "vg annotate -x your_graph.xg -f $OUTPUT_DIR/*_prokka/*.gff -g -t $THREADS > pangenome_annotated.ggff"
echo ""
echo "💡 TIP: You can run failed samples individually later"
