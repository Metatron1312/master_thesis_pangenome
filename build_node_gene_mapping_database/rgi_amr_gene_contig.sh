#!/bin/bash
# Script: run_rgi_clean.sh
# RGI ANALYSIS WITH CLEAN OUTPUT AND AUTO-DETECT

# ========== CẤU HÌNH ==========
ROOT_DIR="/Users/shinra/SAMEA8218225/amr_gene_assemblies"
DOCKER_IMAGE="quay.io/biocontainers/rgi:6.0.5--pyh05cac1d_0"
OUTPUT_BASE_DIR="/Users/shinra/SAMEA8218225/amr_gene_contig_rgi_results"
SAMPLE_NAME="SAMEA8218225"

# ========== TỰ ĐỘNG DETECT GENE PREFIX ==========
echo "🔍 Auto-detecting gene prefix in: $ROOT_DIR"

if [ ! -d "$ROOT_DIR" ]; then
    echo "❌ ERROR: Directory not found: $ROOT_DIR"
    exit 1
fi

# Lấy gene prefix từ directories
PREFIX=$(ls -1 "$ROOT_DIR" 2>/dev/null | grep -E '^[A-Z]{4,}_[0-9]+$' | head -1 | cut -d'_' -f1)

if [ -z "$PREFIX" ]; then
    PREFIX=$(ls -1 "$ROOT_DIR" 2>/dev/null | grep -E '^[A-Z]{2,}_[0-9]+' | head -1 | cut -d'_' -f1)
fi

if [ -z "$PREFIX" ]; then
    echo "⚠️  Could not auto-detect gene prefix, checking available directories..."
    ls -la "$ROOT_DIR" | head -10
    echo "💡 Using default: EGPNMC"
    PREFIX="EGPNMC"
fi

echo "✅ Using gene prefix: ${PREFIX}_"

# ========== CHẠY RGI ==========

echo "========================================="
echo "RGI CLEAN ANALYSIS - ${SAMPLE_NAME}"
echo "========================================="

# Tạo thư mục output
OUTPUT_DIR="${OUTPUT_BASE_DIR}/${SAMPLE_NAME}"
mkdir -p "${OUTPUT_DIR}"

SUMMARY_FILE="${OUTPUT_DIR}/rgi_summary.txt"
COMBINED_RESULTS="${OUTPUT_DIR}/${SAMPLE_NAME}_rgi_results.tab"

# Bắt đầu log
echo "RGI ANALYSIS - ${SAMPLE_NAME}" > "${SUMMARY_FILE}"
echo "==================================" >> "${SUMMARY_FILE}"
echo "Analysis Date: $(date)" >> "${SUMMARY_FILE}"
echo "Gene prefix: ${PREFIX}_" >> "${SUMMARY_FILE}"
echo "Settings: --include_loose --low_quality" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

# Hàm tìm file FASTA
find_fasta_for_gene() {
    local gene_dir="$1"
    local gene_id="$2"
    
    # Ưu tiên 1: file chuẩn
    if [ -f "${gene_dir}/${gene_id}_contigs.fa" ]; then
        echo "${gene_dir}/${gene_id}_contigs.fa"
        return 0
    fi
    
    # Ưu tiên 2: file với .fasta extension
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
    
    # Ưu tiên 4: tìm file có chứa "contig"
    found_file=$(find "${gene_dir}" -maxdepth 1 -name "*contig*" -type f 2>/dev/null | head -1)
    
    if [ -n "${found_file}" ]; then
        echo "${found_file}"
        return 0
    fi
    
    return 1
}

# Tìm các gene directories
echo "Finding gene directories with pattern: ${PREFIX}_*..."
genes=()
while IFS= read -r -d $'\0' dir; do
    genes+=("$(basename "${dir}")")
done < <(find "${ROOT_DIR}" -maxdepth 1 -type d -name "${PREFIX}_*" -print0 | sort -z)

total_genes=${#genes[@]}

if [ "${total_genes}" -eq 0 ]; then
    echo "❌ ERROR: No gene directories found matching ${PREFIX}_*"
    echo "Available directories:"
    find "${ROOT_DIR}" -maxdepth 1 -type d -name "*" | sort
    exit 1
fi

echo "✅ Found ${total_genes} gene directories"

# Tạo header cho file kết quả tổng hợp (đúng cột từ file RGI output)
echo -e "GENE_ID\tBest_Hit_ARO\tBest_Identities(%)\tDrug_Class\tResistance_Mechanism\tCut_Off\tStart\tStop\tOrientation\tContig\tARO\tModel_Type" > "${COMBINED_RESULTS}"

# Biến thống kê
processed=0
with_hits=0
total_hits=0
all_amr_genes=()

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
        echo "❌ ERROR: Directory not found: ${gene_dir}"
        continue
    fi
    
    # Tìm file FASTA
    input_file=$(find_fasta_for_gene "${gene_dir}" "${gene_id}")
    
    if [ -z "${input_file}" ] || [ ! -f "${input_file}" ]; then
        echo "❌ ERROR: No FASTA file found for ${gene_id}"
        echo "   Available files in ${gene_dir}:"
        ls -la "${gene_dir}/" 2>/dev/null | head -10 || echo "   (empty directory)"
        continue
    fi
    
    echo "✅ Found FASTA file: $(basename "${input_file}")"
    echo "   Size: $(ls -lh "${input_file}" | awk '{print $5}')"
    
    # Tạo thư mục output tạm thời cho RGI
    TMP_RGI_DIR="/tmp/rgi_${gene_id}_$(date +%s)"
    mkdir -p "${TMP_RGI_DIR}"
    
    # Chạy RGI trong thư mục tạm
    echo "🔬 Running RGI..."
    
    docker_run_start=$(date +%s)
    
    docker run --rm \
        --platform linux/amd64 \
        -v "${gene_dir}:/input" \
        -v "${TMP_RGI_DIR}:/output" \
        -w "/input" \
        "${DOCKER_IMAGE}" \
        rgi main \
        -i "$(basename "${input_file}")" \
        -o "/output/${gene_id}_rgi" \
        -t contig \
        --include_loose \
        --low_quality \
        --num_threads 1 > "${TMP_RGI_DIR}/rgi.log" 2>&1
    
    docker_exit_code=$?
    docker_run_end=$(date +%s)
    docker_run_time=$((docker_run_end - docker_run_start))
    
    echo "⏱️  RGI completed in ${docker_run_time}s, exit code: ${docker_exit_code}"
    
    if [ ${docker_exit_code} -eq 0 ]; then
        result_file="${TMP_RGI_DIR}/${gene_id}_rgi.txt"
        
        if [ -f "${result_file}" ] && [ -s "${result_file}" ]; then
            # Đếm hits (bỏ header)
            hits=$(tail -n +2 "${result_file}" | wc -l | tr -d ' ')
            echo "✅ Found ${hits} AMR hits"
            
            # Parse kết quả với đúng cột
            hit_count=0
            tail -n +2 "${result_file}" | while IFS=$'\t' read -r line; do
                # Tách các cột
                IFS=$'\t' read -r -a fields <<< "${line}"
                
                if [ ${#fields[@]} -ge 16 ]; then
                    best_hit_aro="${fields[8]}"
                    best_identities="${fields[9]}"
                    drug_class="${fields[14]}"
                    resistance_mech="${fields[15]}"
                    cut_off="${fields[5]}"
                    start="${fields[2]}"
                    stop="${fields[3]}"
                    orientation="${fields[4]}"
                    contig="${fields[1]}"
                    aro="${fields[10]}"
                    model_type="${fields[11]}"
                    
                    # Thêm vào combined results
                    echo -e "${gene_id}\t${best_hit_aro}\t${best_identities}\t${drug_class}\t${resistance_mech}\t${cut_off}\t${start}\t${stop}\t${orientation}\t${contig}\t${aro}\t${model_type}" >> "${COMBINED_RESULTS}"
                    
                    # Thêm vào danh sách genes
                    all_amr_genes+=("${best_hit_aro}")
                    
                    # Hiển thị chi tiết hits
                    ((hit_count++))
                    if [ ${hit_count} -le 3 ]; then
                        echo "   📊 Hit ${hit_count}: ${best_hit_aro}"
                        echo "      Identity: ${best_identities}%, Cutoff: ${cut_off}"
                        if [ -n "${drug_class}" ] && [ "${drug_class}" != "." ]; then
                            echo "      Drug Class: ${drug_class}"
                        fi
                    fi
                fi
            done
            
            if [ ${hit_count} -gt 3 ]; then
                echo "   ... and $((hit_count - 3)) more hits"
            fi
            
            # Copy chỉ file kết quả quan trọng đến output chính
            FINAL_OUTPUT_DIR="${OUTPUT_DIR}/${gene_id}"
            mkdir -p "${FINAL_OUTPUT_DIR}"
            
            cp "${result_file}" "${FINAL_OUTPUT_DIR}/"
            cp "${TMP_RGI_DIR}/rgi.log" "${FINAL_OUTPUT_DIR}/" 2>/dev/null
            
            ((with_hits++))
            total_hits=$((total_hits + hits))
            
            # Ghi vào summary
            echo "${gene_id}: ${hits} hits" >> "${SUMMARY_FILE}"
        else
            echo "ℹ️  No AMR hits found"
            
            # Kiểm tra log
            if [ -f "${TMP_RGI_DIR}/rgi.log" ]; then
                if grep -q "Found 0 genes" "${TMP_RGI_DIR}/rgi.log"; then
                    echo "   Note: RGI explicitly reported 'Found 0 genes'"
                fi
            fi
            
            echo "${gene_id}: 0 hits" >> "${SUMMARY_FILE}"
        fi
        
        ((processed++))
        
        # Dọn dẹp thư mục tạm
        rm -rf "${TMP_RGI_DIR}"
    else
        echo "❌ RGI failed"
        
        # Hiển thị error log
        if [ -f "${TMP_RGI_DIR}/rgi.log" ] && [ -s "${TMP_RGI_DIR}/rgi.log" ]; then
            echo "Error log (first 5 lines):"
            head -5 "${TMP_RGI_DIR}/rgi.log"
        fi
        
        echo "${gene_id}: ERROR" >> "${SUMMARY_FILE}"
        rm -rf "${TMP_RGI_DIR}"
    fi
    
    # Nghỉ giữa các gene
    if [ ${gene_index} -lt ${total_genes} ]; then
        echo ""
        echo "⏳ Waiting 1 second before next gene..."
        sleep 1
    fi
done

# ========== PHÂN TÍCH KẾT QUẢ ==========

# Đếm hits thực tế từ file tổng hợp
actual_hits=0
if [ -f "${COMBINED_RESULTS}" ]; then
    actual_hits=$(tail -n +2 "${COMBINED_RESULTS}" | wc -l)
fi

# Hoàn thành summary
echo "" >> "${SUMMARY_FILE}"
echo "SUMMARY STATISTICS" >> "${SUMMARY_FILE}"
echo "==================" >> "${SUMMARY_FILE}"
echo "Total gene directories: ${total_genes}" >> "${SUMMARY_FILE}"
echo "Genes processed: ${processed}" >> "${SUMMARY_FILE}"
echo "Genes with AMR hits: ${with_hits}" >> "${SUMMARY_FILE}"
echo "Total AMR hits: ${actual_hits}" >> "${SUMMARY_FILE}"

# Phân tích theo cutoff
echo "" >> "${SUMMARY_FILE}"
echo "RESULTS BY CUTOFF STRINGENCY:" >> "${SUMMARY_FILE}"
echo "=============================" >> "${SUMMARY_FILE}"

strict_hits=$(grep -c $'\tStrict\t' "${COMBINED_RESULTS}" 2>/dev/null || echo "0")
loose_hits=$(grep -c $'\tLoose\t' "${COMBINED_RESULTS}" 2>/dev/null || echo "0")
perfect_hits=$(grep -c $'\tPerfect\t' "${COMBINED_RESULTS}" 2>/dev/null || echo "0")

echo "Perfect cutoff hits: ${perfect_hits}" >> "${SUMMARY_FILE}"
echo "Strict cutoff hits: ${strict_hits}" >> "${SUMMARY_FILE}"
echo "Loose cutoff hits: ${loose_hits}" >> "${SUMMARY_FILE}"

# Phân tích AMR genes với identity cao (>95%)
echo "" >> "${SUMMARY_FILE}"
echo "HIGH CONFIDENCE HITS (>95% identity):" >> "${SUMMARY_FILE}"
echo "===================================" >> "${SUMMARY_FILE}"

high_conf_hits=0
if [ -f "${COMBINED_RESULTS}" ] && [ ${actual_hits} -gt 0 ]; then
    tail -n +2 "${COMBINED_RESULTS}" | while IFS=$'\t' read -r gene_id best_hit_aro best_identities drug_class resistance_mech cut_off start stop orientation contig aro model_type; do
        # Loại bỏ % sign nếu có
        identity_num=$(echo "${best_identities}" | sed 's/%//g')
        
        if (( $(echo "${identity_num} > 95" | bc -l 2>/dev/null || echo "0") == 1 )); then
            echo "${gene_id}: ${best_hit_aro} (${best_identities}%, ${cut_off})" >> "${SUMMARY_FILE}"
            ((high_conf_hits++))
        fi
    done
fi

if [ ${high_conf_hits} -eq 0 ]; then
    echo "No hits with >95% identity" >> "${SUMMARY_FILE}"
fi

# Top AMR genes
echo "" >> "${SUMMARY_FILE}"
echo "TOP AMR GENES (all hits):" >> "${SUMMARY_FILE}"
echo "========================" >> "${SUMMARY_FILE}"

if [ -f "${COMBINED_RESULTS}" ] && [ ${actual_hits} -gt 0 ]; then
    tail -n +2 "${COMBINED_RESULTS}" | cut -f2 | sort | uniq -c | sort -rn | head -10 | while read count gene; do
        echo "${gene}: ${count} hits" >> "${SUMMARY_FILE}"
    done
else
    echo "No AMR genes detected" >> "${SUMMARY_FILE}"
fi

# Top drug classes
echo "" >> "${SUMMARY_FILE}"
echo "TOP DRUG CLASSES:" >> "${SUMMARY_FILE}"
echo "================" >> "${SUMMARY_FILE}"

if [ -f "${COMBINED_RESULTS}" ] && [ ${actual_hits} -gt 0 ]; then
    tail -n +2 "${COMBINED_RESULTS}" | cut -f4 | sort | uniq -c | sort -rn | head -5 | while read count class; do
        if [ -n "${class}" ] && [ "${class}" != "." ] && [ "${class}" != "-" ]; then
            echo "${class}: ${count}" >> "${SUMMARY_FILE}"
        fi
    done
fi

# ========== HIỂN THỊ KẾT QUẢ ==========
echo ""
echo "========================================="
echo "📊 ANALYSIS COMPLETED"
echo "========================================="
echo ""
echo "FINAL STATISTICS:"
echo "================="
echo "Sample: ${SAMPLE_NAME}"
echo "Gene prefix: ${PREFIX}_"
echo "Total gene directories: ${total_genes}"
echo "Genes processed: ${processed}"
echo "Genes with AMR hits: ${with_hits}"
echo "Total AMR hits: ${actual_hits}"
echo ""
echo "Cutoff stringency:"
echo "  Perfect: ${perfect_hits}"
echo "  Strict: ${strict_hits}"
echo "  Loose: ${loose_hits}"
echo ""
echo "High confidence hits (>95% identity): ${high_conf_hits}"
echo ""
echo "📁 Output directory: ${OUTPUT_DIR}"
echo ""
echo "📄 MAIN OUTPUT FILES:"
echo "1. ${SAMPLE_NAME}_rgi_results.tab - All results combined (clean format)"
echo "2. rgi_summary.txt - Analysis summary"
echo "3. ${with_hits} gene directories with individual results"
echo ""

# Hiển thị vài kết quả chất lượng cao
if [ ${actual_hits} -gt 0 ]; then
    echo "🔍 HIGH QUALITY RESULTS (Strict/Perfect cutoff, >90% identity):"
    echo "=========================================================="
    
    if [ -f "${COMBINED_RESULTS}" ]; then
        tail -n +2 "${COMBINED_RESULTS}" | awk -F'\t' '
            ($6 == "Strict" || $6 == "Perfect") && $3+0 > 90 {
                printf "  %s: %s (%s%% identity, %s) - %s\n", $1, $2, $3, $6, $4
            }
        ' | head -10
        
        if [ ${perfect_hits} -gt 0 ]; then
            echo ""
            echo "🏆 PERFECT CUTOFF HITS:"
            tail -n +2 "${COMBINED_RESULTS}" | awk -F'\t' '
                $6 == "Perfect" {
                    printf "  %s: %s (%s%% identity)\n", $1, $2, $3
                }
            ' | head -5
        fi
    fi
else
    echo "⚠️  No AMR hits detected in any sequences"
    echo ""
    echo "💡 TROUBLESHOOTING SUGGESTIONS:"
    echo "1. Check FASTA file format and content"
    echo "2. Try lowering quality thresholds"
    echo "3. Verify RGI database is properly installed"
    echo "4. Check individual gene directories for issues"
fi

echo ""
echo "========================================="
echo "✅ RGI CLEAN ANALYSIS COMPLETED!"
echo "========================================="
echo ""
echo "📋 Summary file: ${SUMMARY_FILE}"
echo "📊 Combined results: ${COMBINED_RESULTS}"
echo ""
echo "💡 Next steps:"
echo "1. Review high-confidence hits"
echo "2. Check cutoff stringency distribution"
echo "3. Compare with Abricate results for validation"