#!/bin/bash
# =============================================================================
# EXTRACT GENE SEQUENCES PER SAMPLE - NHẬN THAM SỐ BẮT BUỘC
# =============================================================================
# Script: extract_sequence_from_node_gene_database_v3.sh
# Usage: ./extract_sequence_from_node_gene_database_v3.sh <SAMPLE_ID>
# Example: ./extract_sequence_from_node_gene_database_v3.sh SAMEA8218227
# =============================================================================

# NHẬN SAMPLE ID TỪ THAM SỐ DÒNG LỆNH
if [ $# -eq 0 ]; then
    echo "❌ ERROR: Cần cung cấp SAMPLE_ID"
    echo "Usage: $0 <SAMPLE_ID>"
    echo "Example: $0 SAMEA8218227"
    exit 1
fi

SAMPLE_ID="$1"
BASE_DIR="/Users/shinra"

GRAPH_GBZ="${BASE_DIR}/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/vg_graph/graph/106_strains_with_refseq_complete_vg_autoindex/complete_106_strains_with_refseq_pangenome.giraffe.gbz"
TSV_FILE="${BASE_DIR}/Bioinformatics/WGS_AMR/data_public/Graph_annotation/amr_gene_nodes_database_v10.tsv"
OUTPUT_DIR="${BASE_DIR}/sequence_extraction"
SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/${SAMPLE_ID}"
OUTPUT_FILE="${SAMPLE_OUTPUT_DIR}/${SAMPLE_ID}_AMR_gene_sequences.fasta"
QC_REPORT="${SAMPLE_OUTPUT_DIR}/${SAMPLE_ID}_qc_summary.txt"
EXTRACT_LOG="${SAMPLE_OUTPUT_DIR}/${SAMPLE_ID}_extraction_log.tsv"

mkdir -p "$SAMPLE_OUTPUT_DIR"

echo "============================================="
echo "EXTRACT GENE SEQUENCES - THAM SỐ VERSION"
echo "============================================="
echo "📋 Sample: $SAMPLE_ID"
echo "📋 Output: $OUTPUT_FILE"
echo "============================================="

# Kiểm tra files
[ ! -f "$GRAPH_GBZ" ] && echo "❌ Graph not found: $GRAPH_GBZ" && exit 1
[ ! -f "$TSV_FILE" ] && echo "❌ Database not found: $TSV_FILE" && exit 1

# =============================================================================
# FUNCTION: DETECT GENE PREFIX - CHỈ DỰA VÀO SAMPLE ID
# =============================================================================
detect_gene_prefix_from_database() {
    local sample="$1"
    
    echo "🔍 Đang kiểm tra sample $sample trong database..." >&2
    
    # Đếm số genes cho sample này
    local gene_count=$(awk -F'\t' -v s="$sample" '$3 == s {count++} END {print count}' "$TSV_FILE")
    
    if [ "$gene_count" -eq 0 ]; then
        echo "❌ ERROR: KHÔNG TÌM THẤY SAMPLE $sample TRONG DATABASE!" >&2
        echo "" >&2
        echo "📋 CÁC SAMPLES CÓ TRONG DATABASE:" >&2
        echo "----------------------------------------" >&2
        awk -F'\t' 'NR>1 {print $3}' "$TSV_FILE" | sort -u | head -10 | while read s; do
            c=$(awk -F'\t' -v samp="$s" '$3 == samp {count++} END {print count}' "$TSV_FILE")
            printf "  • %-15s %3d genes\n" "$s" "$c" >&2
        done
        exit 1
    fi
    
    echo "✅ Tìm thấy $gene_count genes cho sample $sample" >&2
    
    # Lấy gene IDs để xác định prefix
    local gene_ids=$(awk -F'\t' -v s="$sample" '$3 == s {print $1; if (NR>=5) exit}' "$TSV_FILE")
    
    echo "📋 Một số genes:" >&2
    echo "$gene_ids" | head -3 | while read gene; do
        echo "  • $gene" >&2
    done
    
    # Lấy prefix từ gene đầu tiên
    local first_gene=$(awk -F'\t' -v s="$sample" '$3 == s {print $1; exit}' "$TSV_FILE")
    local prefix=$(echo "$first_gene" | cut -d'_' -f1)
    
    echo "$prefix"
}

# =============================================================================
# DETECT PREFIX
# =============================================================================
echo ""
echo "🔍 Auto-detecting gene prefix from database..."
GENE_PREFIX=$(detect_gene_prefix_from_database "$SAMPLE_ID")
echo "✅ Detected prefix: ${GENE_PREFIX}_"

# Lấy vài gene samples để hiển thị
SAMPLE_GENES=$(awk -F'\t' -v s="$SAMPLE_ID" '$3 == s {print $1; if (NR>=3) exit}' "$TSV_FILE" | tr '\n' ' ')
echo "   Sample genes: $SAMPLE_GENES"
echo ""

# Kiểm tra đã extract chưa
if [ -f "$OUTPUT_FILE" ] && [ -s "$OUTPUT_FILE" ]; then
    EXISTING_COUNT=$(grep -c "^>" "$OUTPUT_FILE")
    echo "⚠️  FASTA already exists with $EXISTING_COUNT sequences"
    read -p "   Overwrite? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "❌ Extraction cancelled"
        exit 0
    fi
    echo "✅ Will overwrite existing file"
fi

# Tạo file mới
> "$OUTPUT_FILE"
> "$QC_REPORT"
> "$EXTRACT_LOG"

# Tạo extraction log header
echo -e "gene_id\tgene_name\tsample_id\textraction_date\textraction_status\tsequence_length\tnode_count\toriginal_identity\tvalidation_identity\tqc_status" > "$EXTRACT_LOG"

# Parse database để lấy genes cho sample này
echo "📋 Reading database for sample $SAMPLE_ID with prefix ${GENE_PREFIX}_..."

TEMP_FILE="${SAMPLE_OUTPUT_DIR}/temp_genes_$$.txt"
> "$TEMP_FILE"

python3 -c "
import csv
import sys
from datetime import datetime

tsv_file = '$TSV_FILE'
sample_id = '$SAMPLE_ID'
gene_prefix = '$GENE_PREFIX'
output_file = '$TEMP_FILE'
qc_report_file = '$QC_REPORT'

genes = []
qc_failed_genes = []
qc_passed_genes = []
qc_no_data_genes = []
found_genes = []

try:
    with open(tsv_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        
        for row in reader:
            if len(row) < 15:
                continue
            
            # Kiểm tra sample ID và có nodes
            if row[2] == sample_id and row[14].strip() in ['YES', 'HAS_NODES'] and row[7] and row[7].strip() != 'N/A':
                gene_id = row[0].strip()
                gene_name = row[1].strip()
                node_sequence = row[7].strip()
                node_count = row[8].strip()
                orig_id = row[9].strip()
                val_id = row[10].strip()
                qc_status = row[11].strip()
                
                # Kiểm tra gene_id có đúng detected prefix không
                if gene_id.startswith(gene_prefix + '_'):
                    found_genes.append(gene_id)
                    full_name = f'{gene_id}_{gene_name}'
                    genes.append(f'{full_name}|{node_sequence}|{node_count}|{orig_id}|{val_id}|{qc_status}')
                    
                    # Phân loại QC
                    if qc_status.startswith('PASSED'):
                        qc_passed_genes.append(f'{gene_id} ({gene_name}): {orig_id}% vs {val_id}%')
                    elif qc_status.startswith('FAILED'):
                        if '(' in qc_status and ')' in qc_status:
                            diff = qc_status.split('(')[1].split(')')[0]
                        else:
                            diff = 'N/A'
                        qc_failed_genes.append(f'{gene_id} ({gene_name}): {orig_id}% vs {val_id}% - Diff: {diff}')
                    elif qc_status == 'NO_DATA':
                        qc_no_data_genes.append(f'{gene_id} ({gene_name}): No original identity')
    
    # Ghi genes vào temp file
    with open(output_file, 'w', encoding='utf-8') as f:
        for gene in genes:
            f.write(gene + '\\n')
    
    # Ghi QC report
    current_date = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(qc_report_file, 'w', encoding='utf-8') as f:
        f.write('=' * 50 + '\\n')
        f.write('GENE EXTRACTION QC REPORT\\n')
        f.write('=' * 50 + '\\n')
        f.write(f'Sample ID: {sample_id}\\n')
        f.write(f'Date: {current_date}\\n')
        f.write(f'Gene Prefix: {gene_prefix}_\\n')
        f.write(f'Total genes found: {len(genes)}\\n\\n')
        
        f.write('SUMMARY\\n')
        f.write('-' * 30 + '\\n')
        f.write(f'QC PASSED: {len(qc_passed_genes)}\\n')
        f.write(f'QC FAILED: {len(qc_failed_genes)}\\n')
        f.write(f'QC NO_DATA: {len(qc_no_data_genes)}\\n\\n')
        
        if qc_failed_genes:
            f.write('GENES WITH QC FAILED\\n')
            f.write('-' * 30 + '\\n')
            for gene in qc_failed_genes[:5]:
                f.write(f'  • {gene}\\n')
            if len(qc_failed_genes) > 5:
                f.write(f'  ... and {len(qc_failed_genes)-5} more\\n')
            f.write('\\n')
        
    print(f'Found {len(genes)} genes for sample {sample_id} with prefix {gene_prefix}_')
    print(f'QC PASSED: {len(qc_passed_genes)}')
    print(f'QC FAILED: {len(qc_failed_genes)}')
    print(f'QC NO_DATA: {len(qc_no_data_genes)}')
    
    if len(genes) == 0:
        print(f'❌ No genes found for this sample with prefix {gene_prefix}_')
        sys.exit(1)
        
except Exception as e:
    print(f'Error: {e}')
    sys.exit(1)
"

PYTHON_EXIT=$?
if [ $PYTHON_EXIT -ne 0 ]; then
    echo "❌ Error processing database"
    rm -f "$TEMP_FILE"
    exit 1
fi

GENES_COUNT=$(wc -l < "$TEMP_FILE" | tr -d ' ')
echo "✅ Found $GENES_COUNT genes to extract"

# Hiển thị genes
echo ""
echo "=== GENES TO EXTRACT ==="
counter=0
while IFS='|' read -r name nodes count orig val qc && [ $counter -lt 5 ]; do
    node_num=$(echo "$nodes" | tr ',' '\n' | wc -l | tr -d ' ')
    qc_symbol="✅"
    [[ $qc == FAILED* ]] && qc_symbol="❌"
    [[ $qc == "NO_DATA" ]] && qc_symbol="⚠️"
    gene_id=$(echo "$name" | cut -d'_' -f1)
    echo "  $qc_symbol $gene_id: $node_num nodes"
    ((counter++))
done < "$TEMP_FILE"

if [ $GENES_COUNT -gt 5 ]; then
    echo "  ... and $((GENES_COUNT - 5)) more"
fi

# Extract sequences
echo ""
echo "🔬 Extracting sequences from graph..."
echo "============================================="

EXTRACTED=0
FAILED=0
TOTAL_BP=0
CURRENT_DATE=$(date +"%Y-%m-%d %H:%M:%S")

while IFS='|' read -r name nodes count orig val qc; do
    gene_id=$(echo "$name" | cut -d'_' -f1)
    gene_name=$(echo "$name" | cut -d'_' -f2-)
    
    echo ""
    echo "🔬 Gene: $gene_id"
    echo "   Name: $gene_name"
    echo "   Nodes: $count"
    
    # Clean nodes
    CLEAN_NODES=$(echo "$nodes" | tr -cd '0-9,')
    
    # Tạo awk command
    AWK_CMD='
    BEGIN {
        split(nodes, arr, ",")
        for (i in arr) target[arr[i]] = 1
        seq = ""
    }
    /^S/ && $2 in target {
        seq = seq $3
    }
    END {
        if (seq != "") print seq
    }'
    
    # Extract sequence
    SEQUENCE=""
    if command -v vg &> /dev/null; then
        SEQUENCE=$(vg view "$GRAPH_GBZ" 2>/dev/null | awk -v nodes="$CLEAN_NODES" "$AWK_CMD")
    fi
    
    if [ -z "$SEQUENCE" ]; then
        echo "   ⚠️  Failed to extract"
        echo -e "$gene_id\t$gene_name\t$SAMPLE_ID\t$CURRENT_DATE\tFAILED\t0\t$count\t$orig\t$val\t$qc" >> "$EXTRACT_LOG"
        ((FAILED++))
        continue
    fi
    
    SEQ_LENGTH=$(echo "$SEQUENCE" | tr -d '\n' | wc -c | tr -d ' ')
    TOTAL_BP=$((TOTAL_BP + SEQ_LENGTH))
    
    # Ghi vào FASTA
    echo ">$name" >> "$OUTPUT_FILE"
    echo "$SEQUENCE" | fold -w 80 >> "$OUTPUT_FILE"
    
    # Ghi log
    echo -e "$gene_id\t$gene_name\t$SAMPLE_ID\t$CURRENT_DATE\tSUCCESS\t$SEQ_LENGTH\t$count\t$orig\t$val\t$qc" >> "$EXTRACT_LOG"
    
    echo "   ✅ Extracted: ${SEQ_LENGTH} bp"
    ((EXTRACTED++))
    
done < "$TEMP_FILE"

# Cleanup
rm -f "$TEMP_FILE"

# Kết quả
echo ""
echo "============================================="
echo "✅ EXTRACTION COMPLETED!"
echo "============================================="
echo "📁 Output file: $OUTPUT_FILE"
echo "📊 QC report: $QC_REPORT"
echo "📋 Extraction log: $EXTRACT_LOG"
echo ""
echo "=== EXTRACTION STATS ==="
echo "Sample: $SAMPLE_ID"
echo "Gene Prefix: ${GENE_PREFIX}_"
echo "Total genes: $GENES_COUNT"
echo "Successfully extracted: $EXTRACTED"
echo "Failed to extract: $FAILED"
echo "Total length: $TOTAL_BP bp"
echo "Success rate: $((EXTRACTED * 100 / GENES_COUNT))%"