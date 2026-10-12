#!/bin/bash

# =============================================
# CẤU HÌNH
# =============================================

GENES_DIR="/Users/shinra/SAMEA8218225/amr_gene_reads_v2/amr_genes_only"
OUTPUT_DIR="/Users/shinra/SAMEA8218225/amr_gene_assemblies"

# Tham số
CPU=2
MEM=4
MINLEN=100

# =============================================
# KHỞI TẠO
# =============================================

echo "AMR Gene Assembly Pipeline - NO all_contigs"
echo "==========================================="
echo "Input: $GENES_DIR"
echo "Output: $OUTPUT_DIR"
echo ""

# Tạo thư mục output
mkdir -p "$OUTPUT_DIR"

# Kiểm tra Docker
if ! command -v docker &> /dev/null; then
    echo "❌ ERROR: Docker not found!"
    exit 1
fi

# =============================================
# HÀM CHÍNH - CHỈ GIỮ 1 FILE CONTIGS
# =============================================

assemble_gene() {
    local gene_dir="$1"
    local gene_name=$(basename "$gene_dir")
    local outdir="$OUTPUT_DIR/$gene_name"
    
    echo "🔬 Processing: $gene_name"
    
    # Tìm file reads
    local r1_file=""
    local r2_file=""
    
    for file in "$gene_dir"/*; do
        local fname=$(basename "$file")
        if [[ "$fname" == *"_R1.fastq"* ]] || [[ "$fname" == *"_1.fastq"* ]]; then
            r1_file="$file"
        elif [[ "$fname" == *"_R2.fastq"* ]] || [[ "$fname" == *"_2.fastq"* ]]; then
            r2_file="$file"
        fi
    done
    
    # Kiểm tra file
    if [ -z "$r1_file" ] || [ ! -f "$r1_file" ]; then
        echo "   ❌ ERROR: R1 file not found"
        return 1
    fi
    
    echo "   Found:"
    echo "     R1: $(basename "$r1_file")"
    if [ -n "$r2_file" ] && [ -f "$r2_file" ]; then
        echo "     R2: $(basename "$r2_file")"
        local paired_mode=true
    else
        echo "     R2: Not found (single-end mode)"
        local paired_mode=false
    fi
    
    # Tạo thư mục output
    if [ -d "$outdir" ]; then
        echo "   Cleaning previous output..."
        rm -rf "$outdir"
    fi
    mkdir -p "$outdir"
    
    # KIỂM TRA FILE
    if [ ! -r "$r1_file" ]; then
        echo "   ❌ ERROR: Cannot read R1 file: $r1_file"
        return 1
    fi
    
    if [ "$paired_mode" = true ] && [ ! -r "$r2_file" ]; then
        echo "   ❌ ERROR: Cannot read R2 file: $r2_file"
        return 1
    fi
    
    # TẠO DOCKER COMMAND
    echo "   Running Shovill..."
    
    local docker_cmd="docker run --rm"
    docker_cmd="$docker_cmd --platform linux/amd64"
    docker_cmd="$docker_cmd --memory ${MEM}g"
    docker_cmd="$docker_cmd --cpus $CPU"
    docker_cmd="$docker_cmd -v '$gene_dir':/input"
    docker_cmd="$docker_cmd -v '$outdir':/output"
    docker_cmd="$docker_cmd staphb/shovill"
    docker_cmd="$docker_cmd shovill"
    docker_cmd="$docker_cmd --assembler spades"
    docker_cmd="$docker_cmd --trim"
    docker_cmd="$docker_cmd --cpus $CPU"
    docker_cmd="$docker_cmd --ram $MEM"
    docker_cmd="$docker_cmd --force"
    docker_cmd="$docker_cmd --nostitch"
    
    # Thêm R1/R2
    local r1_basename=$(basename "$r1_file")
    docker_cmd="$docker_cmd --R1 /input/$r1_basename"
    
    if [ "$paired_mode" = true ]; then
        local r2_basename=$(basename "$r2_file")
        docker_cmd="$docker_cmd --R2 /input/$r2_basename"
    fi
    
    docker_cmd="$docker_cmd --outdir /output"
    
    echo "   Command: shovill --nostitch --force"
    
    # Chạy Docker
    if eval "$docker_cmd" 2>&1 | tee "$outdir/shovill.log"; then
        # KIỂM TRA THÀNH CÔNG
        if tail -5 "$outdir/shovill.log" 2>/dev/null | grep -q "Done\." && [ -f "$outdir/contigs.fa" ]; then
            echo "   ✅ Assembly successful"
            
            # Tên file final
            local final_contig_file="${gene_name}_contigs.fa"
            
            # Cập nhật header trong FASTA và lưu thành file mới
            echo "   ✏️  Creating $final_contig_file with updated headers..."
            
            awk -v gene="$gene_name" '
                /^>/ {
                    if (NR > 1) print ""
                    split($0, parts, " ")
                    printf(">%s_%s", gene, substr(parts[1], 2))
                    for (i=2; i<=length(parts); i++) {
                        printf(" %s", parts[i])
                    }
                    print ""
                    next
                }
                { print }
            ' "$outdir/contigs.fa" > "$outdir/$final_contig_file"
            
            # XÓA FILE contigs.fa GỐC để chỉ còn 1 file
            echo "   🗑️  Removing original contigs.fa"
            rm -f "$outdir/contigs.fa"
            
            # Đếm contigs
            local contig_count=$(grep -c '^>' "$outdir/$final_contig_file" 2>/dev/null || echo "0")
            local max_len=$(awk '/^>/ {if (seqlen) print seqlen; seqlen=0; next} {seqlen+=length($0)} END {print seqlen}' "$outdir/$final_contig_file" 2>/dev/null | sort -nr | head -1)
            local total_len=$(grep -v '^>' "$outdir/$final_contig_file" | tr -d '\n' | wc -c 2>/dev/null || echo "0")
            
            echo "   📊 Stats: $contig_count contigs, max length: ${max_len:-0} bp, total: ${total_len:-0} bp"
            
            # Tạo file stats
            {
                echo "Gene: $gene_name"
                echo "Contigs: $contig_count"
                echo "Max contig length: $max_len bp"
                echo "Total assembly length: $total_len bp"
                echo "Date: $(date)"
                echo ""
                echo "Contig lengths (bp):"
                awk '/^>/ {if (seqlen) print seqlen; seqlen=0; next} {seqlen+=length($0)} END {print seqlen}' "$outdir/$final_contig_file" | sort -nr
                echo ""
                echo "FASTA headers:"
                grep '^>' "$outdir/$final_contig_file"
            } > "$outdir/${gene_name}_assembly_stats.txt"
            
            # Thông tin từ log
            if grep -q "It contains" "$outdir/shovill.log"; then
                echo "   📋 " $(grep "It contains" "$outdir/shovill.log" | sed 's/\[shovill\] //')
            fi
            
            return 0
        else
            echo "   ❌ ERROR: Assembly failed"
            echo "   Last 5 lines of log:"
            tail -5 "$outdir/shovill.log"
            return 1
        fi
    else
        echo "   ❌ ERROR: Docker command failed"
        return 1
    fi
}

# =============================================
# CHẠY ASSEMBLY
# =============================================

echo "Scanning for gene directories..."

# Lấy tất cả gene
gene_list=()
while IFS= read -r -d $'\0' dir; do
    gene_list+=("$dir")
done < <(find "$GENES_DIR" -maxdepth 1 -type d -name "*" ! -name "$(basename "$GENES_DIR")" -print0 | sort -z)

total_genes=${#gene_list[@]}

if [ $total_genes -eq 0 ]; then
    echo "❌ ERROR: No gene directories found"
    exit 1
fi

echo "Found $total_genes gene directories"
echo ""

# Bắt đầu
START_TIME=$SECONDS

gene_count=0
success_count=0
failed_genes=()

echo "Starting assembly of $total_genes genes..."
echo ""

for gene_dir in "${gene_list[@]}"; do
    if [ ! -d "$gene_dir" ]; then
        continue
    fi
    
    ((gene_count++))
    gene_name=$(basename "$gene_dir")
    
    echo ""
    echo "========================================"
    echo "GENE $gene_count/$total_genes: $gene_name"
    echo "========================================"
    
    if assemble_gene "$gene_dir"; then
        ((success_count++))
        echo "✅ Success"
    else
        failed_genes+=("$gene_name")
        echo "❌ Failed"
    fi
    
    # Progress
    elapsed_time=$((SECONDS - START_TIME))
    avg_time=$((elapsed_time / gene_count))
    remaining_genes=$((total_genes - gene_count))
    remaining_time=$((avg_time * remaining_genes))
    
    echo ""
    echo "📈 Progress: $gene_count/$total_genes"
    echo "⏱️  Elapsed: $(date -u -d @$elapsed_time +'%H:%M:%S')"
    echo "⏳ Remaining: ~$(date -u -d @$remaining_time +'%H:%M:%S')"
done

# Tổng kết
TOTAL_TIME=$((SECONDS - START_TIME))

echo ""
echo "========================================"
echo "FINAL SUMMARY"
echo "========================================"
echo "Total genes processed: $gene_count"
echo "Successful assemblies: $success_count"
echo "Failed assemblies: ${#failed_genes[@]}"
echo "Success rate: $((success_count * 100 / gene_count))%"
echo "Total time: $(date -u -d @$TOTAL_TIME +'%H:%M:%S')"

# =============================================
# TẠO BÁO CÁO (KHÔNG tạo all_contigs)
# =============================================

REPORT_FILE="$OUTPUT_DIR/assembly_report_$(date +%Y%m%d_%H%M%S).txt"

{
    echo "AMR Gene Assembly Report"
    echo "========================"
    echo "Date: $(date)"
    echo "Input: $GENES_DIR"
    echo "Output: $OUTPUT_DIR"
    echo "Total time: $(date -u -d @$TOTAL_TIME +'%H:%M:%S')"
    echo ""
    echo "SUMMARY"
    echo "-------"
    echo "Total genes: $total_genes"
    echo "Successful: $success_count"
    echo "Failed: ${#failed_genes[@]}"
    echo "Success rate: $((success_count * 100 / total_genes))%"
    echo ""
    
    if [ ${#failed_genes[@]} -gt 0 ]; then
        echo "FAILED GENES"
        echo "------------"
        for gene in "${failed_genes[@]}"; do
            echo "- $gene"
        done
        echo ""
    fi
    
    echo "SUCCESSFUL ASSEMBLIES"
    echo "---------------------"
    for gene_dir in "$OUTPUT_DIR"/*/; do
        gene_name=$(basename "$gene_dir")
        contig_file="${gene_name}_contigs.fa"
        
        if [ -f "$gene_dir/$contig_file" ]; then
            contig_count=$(grep -c '^>' "$gene_dir/$contig_file" 2>/dev/null || echo "0")
            total_len=$(grep -v '^>' "$gene_dir/$contig_file" | tr -d '\n' | wc -c 2>/dev/null || echo "0")
            echo "$gene_name: $contig_count contig(s), $total_len bp"
        fi
    done
    echo ""
    echo "OUTPUT STRUCTURE"
    echo "----------------"
    echo "Each gene directory contains:"
    echo "  - {gene_name}_contigs.fa       # Main contig file"
    echo "  - {gene_name}_assembly_stats.txt"
    echo "  - shovill.log"
    echo ""
    echo "For mapping contigs to graph:"
    echo "  Use individual gene directories"
    echo ""
    echo "Example:"
    echo "  $OUTPUT_DIR/HFJCAC_03307/HFJCAC_03307_contigs.fa"
} > "$REPORT_FILE"

echo ""
echo "📋 Report saved to: $(basename "$REPORT_FILE")"
echo ""
echo "✅ AMR GENE ASSEMBLY COMPLETED!"
echo ""
echo "📁 Output structure:"
echo "  $OUTPUT_DIR/"
echo "  ├── HFJCAC_00054/"
echo "  │   ├── HFJCAC_00054_contigs.fa"
echo "  │   └── HFJCAC_00054_assembly_stats.txt"
echo "  ├── HFJCAC_00093/"
echo "  │   ├── HFJCAC_00093_contigs.fa"
echo "  │   └── HFJCAC_00093_assembly_stats.txt"
echo "  └── ..."
echo ""
echo "🎯 Next step: Use gene directories for mapping to graph"