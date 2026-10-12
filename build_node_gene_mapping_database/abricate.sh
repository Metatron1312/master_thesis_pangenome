#!/bin/bash
# Script: run_abricate_complete.sh
# COMPLETE ABRICATE ANALYSIS - WITH AUTO-DETECT GENE PREFIX

# ========== CẤU HÌNH ĐƯỜNG DẪN ==========
ROOT_DIR="/Users/shinra/SAMEA8218225/amr_gene_assemblies"
DOCKER_IMAGE="staphb/abricate:latest"
OUTPUT_BASE_DIR="/Users/shinra/SAMEA8218225/amr_gene_contig_abricate_results"
SAMPLE_NAME="SAMEA8218225"

# ========== TỰ ĐỘNG DETECT GENE PREFIX ==========
echo "🔍 Auto-detecting gene prefix in: $ROOT_DIR"

# Kiểm tra thư mục tồn tại
if [ ! -d "$ROOT_DIR" ]; then
    echo "❌ ERROR: Directory not found: $ROOT_DIR"
    exit 1
fi

# Lấy tất cả directories, tìm pattern có dạng ABC_123
PREFIX=$(ls -1 "$ROOT_DIR" 2>/dev/null | grep -E '^[A-Z]{4,}_[0-9]+$' | head -1 | cut -d'_' -f1)

if [ -z "$PREFIX" ]; then
    # Thử tìm pattern khác (4+ chữ cái, gạch dưới, rồi số)
    PREFIX=$(ls -1 "$ROOT_DIR" 2>/dev/null | grep -E '^[A-Z]{2,}_[0-9]+' | head -1 | cut -d'_' -f1)
fi

if [ -z "$PREFIX" ]; then
    echo "⚠️  Could not auto-detect gene prefix, using default: EGPNMC"
    PREFIX="EGPNMC"
fi

echo "✅ Using gene prefix: ${PREFIX}_"

# ========== CHẠY ABRICATE HOÀN CHỈNH ==========
echo "========================================="
echo "COMPLETE ABRICATE ANALYSIS - ${SAMPLE_NAME}"
echo "========================================="

# Tạo thư mục output đơn giản
OUTPUT_DIR="${OUTPUT_BASE_DIR}/${SAMPLE_NAME}"
mkdir -p "${OUTPUT_DIR}"

SUMMARY_FILE="${OUTPUT_DIR}/abricate_summary.txt"
COMBINED_RESULTS="${OUTPUT_DIR}/${SAMPLE_NAME}_abricate_results.tab"

# Hàm tạo header
print_header() {
    echo "======================================================="
    echo "$1"
    echo "======================================================="
}

# Bắt đầu log
print_header "ABRICATE ANALYSIS STARTED - $(date)"
echo "Sample: ${SAMPLE_NAME}"
echo "Root directory: ${ROOT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Gene prefix detected: ${PREFIX}_"
echo ""

# Hàm tìm file FASTA
find_fasta_for_gene() {
    local gene_dir="$1"
    local gene_id="$2"
    
    # Ưu tiên 1: file chuẩn
    if [ -f "${gene_dir}/${gene_id}_contigs.fa" ]; then
        echo "${gene_dir}/${gene_id}_contigs.fa"
        return 0
    fi
    
    # Ưu tiên 2: file với .fa extension
    if [ -f "${gene_dir}/${gene_id}_contigs.fasta" ]; then
        echo "${gene_dir}/${gene_id}_contigs.fasta"
        return 0
    fi
    
    # Ưu tiên 3: tìm bất kỳ file .fa hoặc .fasta
    local found_file=$(find "${gene_dir}" -maxdepth 1 \( -name "*.fa" -o -name "*.fasta" \) -type f 2>/dev/null | head -1)
    
    if [ -n "${found_file}" ]; then
        echo "${found_file}"
        return 0
    fi
    
    # Ưu tiên 4: tìm bất kỳ file nào chứa "contig"
    found_file=$(find "${gene_dir}" -maxdepth 1 -name "*contig*" -type f 2>/dev/null | head -1)
    
    if [ -n "${found_file}" ]; then
        echo "${found_file}"
        return 0
    fi
    
    return 1
}

# Tìm các gene directories
echo "Scanning for gene directories with pattern: ${PREFIX}_*"
genes=()
while IFS= read -r -d $'\0' dir; do
    genes+=("$(basename "${dir}")")
done < <(find "${ROOT_DIR}" -maxdepth 1 -type d -name "${PREFIX}_*" -print0 | sort -z)

total_genes=${#genes[@]}

if [ "${total_genes}" -eq 0 ]; then
    echo "❌ ERROR: No gene directories found matching ${PREFIX}_* in ${ROOT_DIR}"
    echo "Available items:"
    ls -la "${ROOT_DIR}/" | head -20
    
    # Thử liệt kê tất cả directories
    echo ""
    echo "All directories found:"
    find "${ROOT_DIR}" -maxdepth 1 -type d -name "*" | sort
    
    # Gợi ý
    echo ""
    echo "💡 SUGGESTIONS:"
    echo "1. Check if the directory path is correct"
    echo "2. Gene prefix might be different. Available prefixes:"
    ls -1 "${ROOT_DIR}" | grep -E '^[A-Z]+_[0-9]' | cut -d'_' -f1 | sort -u | head -5
    exit 1
fi

echo "✅ Found ${total_genes} gene directories"

# Tạo header cho file kết quả tổng hợp
echo "# Abricate Complete Results - Sample: ${SAMPLE_NAME}" > "${COMBINED_RESULTS}"
echo "# Analysis Date: $(date)" >> "${COMBINED_RESULTS}"
echo "# Database: CARD (built-in)" >> "${COMBINED_RESULTS}"
echo "# Settings: --minid 80 --mincov 80" >> "${COMBINED_RESULTS}"
echo -e "GENE_ID\tFILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE" >> "${COMBINED_RESULTS}"

# Biến thống kê
processed=0
with_hits=0
total_hits=0
all_amr_genes=()

print_header "PROCESSING GENES"

# Xử lý từng gene
gene_index=0
for gene_id in "${genes[@]}"; do
    ((gene_index++))
    
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  GENE ${gene_index}/${total_genes}: ${gene_id}"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    gene_dir="${ROOT_DIR}/${gene_id}"
    
    if [ ! -d "${gene_dir}" ]; then
        echo "❌ ERROR: Gene directory not found: ${gene_dir}"
        continue
    fi
    
    # Tìm file FASTA
    echo "🔍 Looking for FASTA file..."
    input_file=$(find_fasta_for_gene "${gene_dir}" "${gene_id}")
    
    if [ -z "${input_file}" ] || [ ! -f "${input_file}" ]; then
        echo "❌ ERROR: No FASTA file found for ${gene_id}"
        echo "   Available files in ${gene_dir}:"
        ls -la "${gene_dir}/" 2>/dev/null | head -10 || echo "   (empty directory)"
        continue
    fi
    
    echo "✅ Found FASTA file: $(basename "${input_file}")"
    echo "   Size: $(ls -lh "${input_file}" | awk '{print $5}')"
    
    # Tạo file output cho gene này
    output_file="${OUTPUT_DIR}/${gene_id}_abricate_results.tab"
    
    # Chạy Abricate
    echo "🔬 Running Abricate on $(basename "${input_file}")..."
    
    docker_run_start=$(date +%s)
    
    docker run --rm \
        --platform linux/amd64 \
        -v "${gene_dir}:/data" \
        -w "/data" \
        "${DOCKER_IMAGE}" \
        abricate --db card \
                 --minid 80 \
                 --mincov 80 \
                 "$(basename "${input_file}")" > "${output_file}" 2>&1
    
    docker_exit_code=$?
    docker_run_end=$(date +%s)
    docker_run_time=$((docker_run_end - docker_run_start))
    
    echo "⏱️  Abricate completed in ${docker_run_time}s, exit code: ${docker_exit_code}"
    
    # Kiểm tra kết quả
    if [ ${docker_exit_code} -eq 0 ] && [ -s "${output_file}" ]; then
        # Đếm số hits (dòng có tab separator)
        hits=$(grep -c $'\t' "${output_file}")
        
        if [ "${hits}" -gt 0 ]; then
            echo "✅ SUCCESS: Found ${hits} AMR hits"
            
            # Xử lý từng hit
            hit_count=0
            while IFS=$'\t' read -r -a fields; do
                if [ ${#fields[@]} -ge 14 ]; then
                    gene_name="${fields[5]}"
                    coverage="${fields[9]}"
                    identity="${fields[10]}"
                    product="${fields[13]}"
                    resistance="${fields[14]}"
                    
                    # Thêm vào combined results
                    echo -e "${gene_id}\t$(echo "${fields[*]}" | tr ' ' '\t')" >> "${COMBINED_RESULTS}"
                    
                    # Ghi log chi tiết
                    if [ "${gene_name}" != "GENE" ]; then
                        echo "   📊 Hit $((++hit_count)): ${gene_name}"
                        echo "      Identity: ${identity}%, Coverage: ${coverage}%"
                        if [ -n "${product}" ] && [ "${product}" != "." ]; then
                            echo "      Product: ${product}"
                        fi
                        all_amr_genes+=("${gene_name}")
                    else
                        echo "   ⚠️  Hit $((++hit_count)): Gene name labeled as 'GENE'"
                        echo "      Identity: ${identity}%, Coverage: ${coverage}%"
                        echo "      Product field: ${product}"
                        
                        # Cố gắng extract gene name từ product
                        if [[ "${product}" == *"armA"* ]] || [[ "${product}" == *"ArmA"* ]]; then
                            echo "      ↳ Extracted: armA from product field"
                            all_amr_genes+=("armA")
                        elif [[ "${product}" == *"ade"* ]] || [[ "${product}" == *"Ade"* ]]; then
                            gene_extracted=$(echo "${product}" | grep -o "ade[A-Z]" | head -1 | tr '[:upper:]' '[:lower:]' 2>/dev/null || echo "ade")
                            echo "      ↳ Extracted: ${gene_extracted} from product field"
                            all_amr_genes+=("${gene_extracted}")
                        elif [[ "${product}" == *"tet"* ]] || [[ "${product}" == *"Tet"* ]]; then
                            gene_extracted=$(echo "${product}" | grep -oi "tet[A-Z0-9()\-]*" | head -1)
                            echo "      ↳ Extracted: ${gene_extracted} from product field"
                            all_amr_genes+=("${gene_extracted}")
                        fi
                    fi
                fi
            done < <(grep $'\t' "${output_file}")
            
            ((with_hits++))
            total_hits=$((total_hits + hits))
            
        else
            echo "ℹ️  No AMR hits found"
            
            # Kiểm tra xem có message "Found 0 genes" không
            if grep -q "Found 0 genes" "${output_file}"; then
                echo "   Note: Abricate explicitly reported 'Found 0 genes'"
            fi
        fi
    else
        echo "❌ ERROR: Abricate failed or produced empty output"
        
        # Hiển thị error nếu có
        if [ -f "${output_file}" ] && [ -s "${output_file}" ]; then
            echo "Error output (first 5 lines):"
            head -5 "${output_file}"
        fi
    fi
    
    ((processed++))
    
    # Nghỉ giữa các gene để tránh overload
    if [ ${gene_index} -lt ${total_genes} ]; then
        echo ""
        echo "⏳ Waiting 1 second before next gene..."
        sleep 1
    fi
done

# ========== PHÂN TÍCH KẾT QUẢ ==========

print_header "ANALYZING RESULTS"

# Đếm hits thực tế từ file tổng hợp (bỏ header lines)
actual_hits=0
if [ -f "${COMBINED_RESULTS}" ]; then
    actual_hits=$(grep -c $'\t' "${COMBINED_RESULTS}")
    # Trừ đi header line nếu có
    if [ "${actual_hits}" -gt 0 ]; then
        actual_hits=$((actual_hits - 1))
    fi
fi

# Tạo summary file
cat > "${SUMMARY_FILE}" << EOF
COMPLETE ABRICATE ANALYSIS - ${SAMPLE_NAME}
==========================================
Analysis Date: $(date)
Database: CARD (built-in)
Settings: --minid 80 --mincov 80
Gene prefix: ${PREFIX}_

STATISTICS:
-----------
Total gene directories: ${total_genes}
Genes processed: ${processed}
Genes with AMR hits: ${with_hits}
Total AMR hits found: ${actual_hits}

OUTPUT FILES:
-------------
1. ${SAMPLE_NAME}_abricate_results.tab - All AMR gene hits combined
2. abricate_summary.txt - This summary file
3. *gene_id*_abricate_results.tab - Individual gene results

EOF

# Phân tích chi tiết nếu có hits
if [ ${actual_hits} -gt 0 ]; then
    echo "DETAILED ANALYSIS:" >> "${SUMMARY_FILE}"
    echo "-----------------" >> "${SUMMARY_FILE}"
    
    # Đếm tần suất AMR genes
    echo "" >> "${SUMMARY_FILE}"
    echo "AMR GENE FREQUENCY:" >> "${SUMMARY_FILE}"
    echo "------------------" >> "${SUMMARY_FILE}"
    
    grep $'\t' "${COMBINED_RESULTS}" | tail -n +2 | cut -f7 | sort | uniq -c | sort -rn | while read count gene; do
        if [ "${gene}" != "GENE" ]; then
            echo "${gene}: ${count} hits" >> "${SUMMARY_FILE}"
        fi
    done
    
    # Genes labeled as "GENE"
    gene_count=$(grep $'\t' "${COMBINED_RESULTS}" | tail -n +2 | cut -f7 | grep -c "GENE")
    if [ ${gene_count} -gt 0 ]; then
        echo "" >> "${SUMMARY_FILE}"
        echo "Hits labeled as 'GENE': ${gene_count}" >> "${SUMMARY_FILE}"
        echo "These may need manual inspection of the PRODUCT field." >> "${SUMMARY_FILE}"
        
        # Liệt kê các product field cho genes labeled as "GENE"
        echo "" >> "${SUMMARY_FILE}"
        echo "PRODUCT fields for 'GENE' entries:" >> "${SUMMARY_FILE}"
        grep $'\t' "${COMBINED_RESULTS}" | awk -F'\t' '$7 == "GENE" {print $15}' | sort | uniq -c | sort -rn | head -10 | while read count product; do
            if [ -n "${product}" ] && [ "${product}" != "." ]; then
                echo "  ${product}: ${count}" >> "${SUMMARY_FILE}"
            fi
        done
    fi
    
    # High confidence hits (>95% identity and coverage)
    echo "" >> "${SUMMARY_FILE}"
    echo "HIGH CONFIDENCE HITS (>95% identity AND >95% coverage):" >> "${SUMMARY_FILE}"
    echo "------------------------------------------------------" >> "${SUMMARY_FILE}"
    
    high_conf_count=0
    grep $'\t' "${COMBINED_RESULTS}" | tail -n +2 | while IFS=$'\t' read -r gene_id file sequence start end strand gene coverage coverage_map gaps pct_coverage identity database accession product resistance; do
        if (( $(echo "${identity} > 95" | bc -l 2>/dev/null || echo "0") == 1 )) && (( $(echo "${pct_coverage} > 95" | bc -l 2>/dev/null || echo "0") == 1 )); then
            echo "${gene_id}: ${gene} (${identity}% identity, ${pct_coverage}% coverage)" >> "${SUMMARY_FILE}"
            ((high_conf_count++))
        fi
    done
    
    if [ ${high_conf_count} -eq 0 ]; then
        echo "No hits meeting high confidence criteria" >> "${SUMMARY_FILE}"
    else
        echo "Total high confidence hits: ${high_conf_count}" >> "${SUMMARY_FILE}"
    fi
    
    # Check for specific known genes
    echo "" >> "${SUMMARY_FILE}"
    echo "SPECIFIC GENE CHECK:" >> "${SUMMARY_FILE}"
    echo "-------------------" >> "${SUMMARY_FILE}"
    
    # Danh sách các gene quan trọng cần kiểm tra
    important_genes=("armA" "tet(B)" "OXA-66" "mphE" "msrE" "abeM" "abeS" "LpsB" "adeA" "adeB" "adeC" "adeF" "adeG" "adeH")
    
    for important_gene in "${important_genes[@]}"; do
        if grep $'\t' "${COMBINED_RESULTS}" | tail -n +2 | cut -f7 | grep -qi "^${important_gene}$"; then
            echo "✅ ${important_gene}: DETECTED in GENE column" >> "${SUMMARY_FILE}"
        elif grep $'\t' "${COMBINED_RESULTS}" | tail -n +2 | cut -f15 | grep -qi "${important_gene}"; then
            echo "⚠️  ${important_gene}: Found in PRODUCT field" >> "${SUMMARY_FILE}"
        else
            echo "❌ ${important_gene}: NOT DETECTED" >> "${SUMMARY_FILE}"
        fi
    done
    
    # Resistance classes
    echo "" >> "${SUMMARY_FILE}"
    echo "RESISTANCE CLASSES:" >> "${SUMMARY_FILE}"
    echo "------------------" >> "${SUMMARY_FILE}"
    
    grep $'\t' "${COMBINED_RESULTS}" | tail -n +2 | cut -f16 | tr ';' '\n' | sort | uniq -c | sort -rn | head -10 | while read count class; do
        if [ -n "${class}" ] && [ "${class}" != "." ] && [ "${class}" != "-" ]; then
            echo "${class}: ${count}" >> "${SUMMARY_FILE}"
        fi
    done
else
    echo "NO AMR GENES DETECTED" >> "${SUMMARY_FILE}"
    echo "--------------------" >> "${SUMMARY_FILE}"
    echo "Abricate did not find any antimicrobial resistance genes." >> "${SUMMARY_FILE}"
    echo "This could be due to:" >> "${SUMMARY_FILE}"
    echo "1. Database limitations (CARD may be outdated)" >> "${SUMMARY_FILE}"
    echo "2. Sequence quality issues" >> "${SUMMARY_FILE}"
    echo "3. Gene sequences not in the database" >> "${SUMMARY_FILE}"
    echo "4. Thresholds too high (--minid 80 --mincov 80)" >> "${SUMMARY_FILE}"
fi

# ========== HIỂN THỊ KẾT QUẢ ==========
print_header "ANALYSIS COMPLETED"

echo ""
echo "📊 FINAL SUMMARY:"
echo "================="
echo "Sample: ${SAMPLE_NAME}"
echo "Total genes analyzed: ${total_genes}"
echo "Genes processed: ${processed}"
echo "Genes with AMR hits: ${with_hits}"
echo "Total AMR hits found: ${actual_hits}"
echo ""
echo "📁 OUTPUT DIRECTORY: ${OUTPUT_DIR}"
echo ""
echo "📄 MAIN OUTPUT FILES:"
echo "1. ${SAMPLE_NAME}_abricate_results.tab - Combined results"
echo "2. abricate_summary.txt - Analysis summary"
echo "3. Individual gene result files (*_abricate_results.tab)"
echo ""

# Hiển thị vài kết quả
if [ ${actual_hits} -gt 0 ]; then
    echo "🔍 SAMPLE OF RESULTS:"
    echo "-------------------"
    grep $'\t' "${COMBINED_RESULTS}" | tail -n +2 | head -5 | while IFS=$'\t' read -r gene_id file sequence start end strand gene coverage coverage_map gaps pct_coverage identity database accession product resistance; do
        echo "  ${gene_id}: ${gene} (${identity}% identity, ${pct_coverage}% coverage)"
    done
    
    # Kiểm tra cụ thể cho armA
    echo ""
    echo "🧬 CHECKING FOR armA:"
    echo "-------------------"
    if grep $'\t' "${COMBINED_RESULTS}" | tail -n +2 | cut -f7 | grep -qi "^armA$"; then
        armA_hits=$(grep $'\t' "${COMBINED_RESULTS}" | tail -n +2 | cut -f7 | grep -ci "^armA$")
        echo "✅ armA detected ${armA_hits} time(s) by Abricate"
        grep $'\t' "${COMBINED_RESULTS}" | awk -F'\t' 'tolower($7) == "arma" {print $1 ": " $11 "% identity, " $10 "% coverage"}' | head -3
    elif grep $'\t' "${COMBINED_RESULTS}" | tail -n +2 | cut -f15 | grep -qi "arma"; then
        echo "⚠️  armA mentioned in PRODUCT field"
        grep $'\t' "${COMBINED_RESULTS}" | awk -F'\t' 'tolower($15) ~ /arma/ {print $1 ": " $15}' | head -3
    else
        echo "❌ armA NOT detected"
    fi
else
    echo "⚠️  No AMR genes detected in any sequences"
    echo "   Consider:"
    echo "   - Lowering thresholds: --minid 70 --mincov 70"
    echo "   - Checking sequence quality"
    echo "   - Trying different database: --db ncbi, resfinder, etc."
fi

echo ""
echo "========================================="
echo "✅ ABRICATE ANALYSIS COMPLETE"
echo "========================================="
echo ""
echo "💡 RECOMMENDATIONS:"
echo "1. Check individual result files for detailed information"
echo "2. For genes labeled as 'GENE', inspect the PRODUCT field"
echo "3. Verify high-confidence hits manually"
echo "4. Consider running with different databases if needed"
echo ""
echo "📋 Summary file: ${SUMMARY_FILE}"
echo "📊 Combined results: ${COMBINED_RESULTS}"