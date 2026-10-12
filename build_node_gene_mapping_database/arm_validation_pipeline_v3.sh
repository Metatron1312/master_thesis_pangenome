#!/bin/bash
# AMR VALIDATION PIPELINE v3 - WITH CUTOFF COLUMNS & GENE_MATCH_QC

# =============================================================================
# CONFIGURATION
# =============================================================================
SAMPLE_ID="SAMEA8218226"
BASE_DIR="/Users/shinra"

# Database
DB_V9="${BASE_DIR}/Bioinformatics/WGS_AMR/data_public/Graph_annotation/amr_gene_nodes_database_v10.tsv"
ORIGINAL_RGI_DIR="${BASE_DIR}/Bioinformatics/WGS_AMR/data_public/Rgi_results"

# Sequence
EXTRACTED_FASTA="${BASE_DIR}/sequence_extraction/${SAMPLE_ID}/${SAMPLE_ID}_AMR_gene_sequences.fasta"

# Current work directory
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
CURRENT_WORK_DIR="${BASE_DIR}/node_gene_mapping_database_validation/${SAMPLE_ID}_${TIMESTAMP}"

# Output directories
ABRICATE_DIR="${CURRENT_WORK_DIR}/abricate_results"
RGI_DIR="${CURRENT_WORK_DIR}/rgi_results" 
VALIDATION_DIR="${CURRENT_WORK_DIR}/validation_results"

# Master database
MASTER_VALIDATION_DB="${BASE_DIR}/node_gene_mapping_database_validation/MASTER_VALIDATION_DATABASE_v9.tsv"

# =============================================================================
# FUNCTIONS
# =============================================================================

detect_gene_prefix() {
    local prefix=""
    
    if [ -f "$DB_V9" ]; then
        local first_gene=$(awk -F'\t' -v s="$SAMPLE_ID" '$3 == s {print $1; exit}' "$DB_V9" 2>/dev/null)
        if [ -n "$first_gene" ]; then
            prefix=$(echo "$first_gene" | grep -o '^[A-Z]\+_[0-9]\+' | cut -d'_' -f1)
            if [ -n "$prefix" ]; then
                echo "$prefix"
                return 0
            fi
        fi
    fi
    
    if [[ "$SAMPLE_ID" == SAMEA8* ]]; then
        echo "EGPNMC"
    else
        echo "HFJCAC"
    fi
}

print_header() {
    echo ""
    echo "========================================="
    echo "$1"
    echo "========================================="
}

print_step() { echo "🔹 $1"; }
print_success() { echo "✅ $1"; }
print_error() { echo "❌ $1"; }
print_warning() { echo "⚠️  $1"; }

normalize_gene() {
    echo "$1" | tr '[:upper:]' '[:lower:]' | sed 's/[^a-z0-9]//g'
}

gene_names_match() {
    local gene1="$1"; local gene2="$2"
    local norm1=$(normalize_gene "$gene1"); local norm2=$(normalize_gene "$gene2")
    if [[ "$norm1" == "$norm2" ]]; then return 0; fi
    if [[ "$norm1" == *"$norm2"* ]] || [[ "$norm2" == *"$norm1"* ]]; then return 0; fi
    return 1
}

get_original_rgi_cutoff() {
    local sample_id="$1"
    local orig_txt="${ORIGINAL_RGI_DIR}/${sample_id}/rgi_results.txt"
    if [ -f "$orig_txt" ]; then
        grep -m1 -oE "Perfect|Strict|Loose" "$orig_txt" 2>/dev/null | head -1
    else
        echo "N/A"
    fi
}

get_gene_match_qc() {
    local sample_id="$1"
    local gene_id="$2"
    local gene_name="$3"
    
    local rgi_file="${BASE_DIR}/${sample_id}/amr_gene_contig_rgi_results/${sample_id}/${gene_id}/${gene_id}_rgi.txt"
    
    if [ ! -f "$rgi_file" ]; then
        echo "NO_RGI_FILE"
        return
    fi
    
    # Dùng Python để tính match value
    python3 << EOF
import re
import sys

rgi_file = "$rgi_file"
gene_name = "$gene_name"

best_gene = ""
best_pct = 0

try:
    with open(rgi_file, 'r') as f:
        lines = f.readlines()
    
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 21:
            try:
                pct = float(parts[20])
                if pct > best_pct:
                    best_pct = pct
                    best_gene = parts[8]
            except:
                pass
    
    if best_gene:
        norm_db = re.sub(r'[^a-z0-9]', '', gene_name.lower())
        norm_rgi = re.sub(r'[^a-z0-9]', '', best_gene.lower())
        
        if norm_rgi == norm_db or norm_db in norm_rgi or norm_rgi in norm_db:
            print("MATCH")
        else:
            print(f"DIFFERENT({best_gene})")
    else:
        print("NO_RGI_HIT")
except Exception as e:
    print("NO_DATA")
EOF
}

# =============================================================================
# STEP 1: SETUP ENVIRONMENT
# =============================================================================
setup_environment() {
    print_header "SETUP ENVIRONMENT"
    
    echo "🔍 Auto-detecting gene prefix..."
    GENE_PREFIX=$(detect_gene_prefix | tail -1)
    echo "✅ Detected gene prefix: ${GENE_PREFIX}_"
    
    rm -rf "$CURRENT_WORK_DIR"
    mkdir -p "$ABRICATE_DIR" "$RGI_DIR" "$VALIDATION_DIR"
    
    if [ ! -f "$DB_V9" ]; then
        print_error "Database not found: $DB_V9"
        exit 1
    fi
    
    if [ ! -f "$EXTRACTED_FASTA" ]; then
        print_error "FASTA not found: $EXTRACTED_FASTA"
        exit 1
    fi
    
    GENES_IN_DB=$(awk -F'\t' -v s="$SAMPLE_ID" '$3 == s {count++} END {print count}' "$DB_V9")
    echo "📊 Sample $SAMPLE_ID has $GENES_IN_DB genes in database"
    
    FASTA_COUNT=$(grep -c "^>" "$EXTRACTED_FASTA" 2>/dev/null || echo "0")
    echo "📊 FASTA has $FASTA_COUNT sequences"
    
    if [ -f "$EXTRACTED_FASTA" ] && [ $FASTA_COUNT -gt 0 ]; then
        echo "📋 FASTA sample headers:"
        head -3 "$EXTRACTED_FASTA" | grep "^>" | sed 's/^/  /'
    fi
    
    print_success "Environment setup complete"
}

# =============================================================================
# STEP 2: RUN RGI
# =============================================================================
run_rgi_and_create_lookup() {
    print_header "RUNNING RGI"
    
    if [ ! -f "$EXTRACTED_FASTA" ] || [ ! -s "$EXTRACTED_FASTA" ]; then
        print_error "FASTA file is empty or missing: $EXTRACTED_FASTA"
        return 1
    fi
    
    RGI_TEMP="/tmp/rgi_${SAMPLE_ID}_${TIMESTAMP}"
    mkdir -p "$RGI_TEMP"
    
    print_step "Running RGI with Docker..."
    
    docker run --rm \
        --platform linux/amd64 \
        -v "$(dirname "$EXTRACTED_FASTA"):/input" \
        -v "$RGI_TEMP:/output" \
        -w "/input" \
        quay.io/biocontainers/rgi:6.0.5--pyh05cac1d_0 \
        rgi main \
        -i "$(basename "$EXTRACTED_FASTA")" \
        -o "/output/rgi_results" \
        -t contig \
        --include_loose \
        --low_quality \
        --num_threads 1 > "${RGI_TEMP}/rgi.log" 2>&1
    
    RGI_RAW="${RGI_DIR}/${SAMPLE_ID}_rgi_raw.txt"
    
    if [ -f "${RGI_TEMP}/rgi_results.txt" ]; then
        cp "${RGI_TEMP}/rgi_results.txt" "$RGI_RAW"
        print_step "RGI completed successfully"
    else
        print_error "RGI failed - no output file"
        tail -10 "${RGI_TEMP}/rgi.log"
        return 1
    fi
    
    print_step "Loading database gene names for matching..."
    
    DB_GENE_FILE="${RGI_DIR}/db_genes.txt"
    awk -F'\t' -v s="$SAMPLE_ID" '$3 == s {print $1 "\t" $2}' "$DB_V9" > "$DB_GENE_FILE"
    
    print_step "Parsing RGI results..."
    
    TEMP_ALL="${RGI_DIR}/all_rgi_entries.tmp"
    > "$TEMP_ALL"
    
    GENE_PREFIX=$(detect_gene_prefix | tail -1)
    
    tail -n +2 "$RGI_RAW" | while IFS= read -r line; do
        IFS=$'\t' read -ra FIELDS <<< "$line"
        
        if [ ${#FIELDS[@]} -ge 16 ]; then
            contig="${FIELDS[1]}"
            cut_off="${FIELDS[5]}"
            best_hit_aro="${FIELDS[8]}"
            best_identities="${FIELDS[9]}"
            
            gene_id=$(echo "$contig" | grep -o "${GENE_PREFIX}_[0-9]*" | head -1)
            if [ -z "$gene_id" ]; then
                gene_id=$(echo "$contig" | grep -o '[A-Z]\+_[0-9]\+' | head -1)
            fi
            
            if [ -n "$gene_id" ] && [ -n "$best_hit_aro" ]; then
                db_gene=$(grep "^$gene_id" "$DB_GENE_FILE" | cut -f2)
                
                gene_match="NO"
                if [ -n "$db_gene" ]; then
                    if gene_names_match "$db_gene" "$best_hit_aro"; then
                        gene_match="YES"
                    fi
                fi
                
                clean_id=$(echo "$best_identities" | sed 's/%//g')
                safe_gene=$(echo "$best_hit_aro" | tr '|' '_')
                safe_db_gene=$(echo "$db_gene" | tr '|' '_')
                
                echo "${gene_id}|${safe_gene}|${clean_id}|${cut_off}|${gene_match}|${safe_db_gene}" >> "$TEMP_ALL"
            fi
        fi
    done
    
    print_step "Applying gene selection logic..."
    
    RGI_LOOKUP="${RGI_DIR}/rgi_lookup.txt"
    echo "# Gene_ID|RGI_Gene|RGI_Identity|RGI_Cutoff|Gene_Match" > "$RGI_LOOKUP"
    
    cut -d'|' -f1 "$TEMP_ALL" | sort -u | while read gene_id; do
        grep "^${gene_id}|" "$TEMP_ALL" > "${TEMP_ALL}.${gene_id}"
        
        best_score=-999999
        best_entry=""
        
        while IFS='|' read -r gid rgi_gene rgi_id rgi_cut gene_match db_gene; do
            int_rgi_id=$(echo "$rgi_id * 100" | bc 2>/dev/null | cut -d'.' -f1)
            [ -z "$int_rgi_id" ] && int_rgi_id=0
            
            match_score=0
            [ "$gene_match" = "YES" ] && match_score=10000
            
            cut_score=0
            [ "$rgi_cut" = "Perfect" ] && cut_score=1000
            [ "$rgi_cut" = "Strict" ] && cut_score=500
            [ "$rgi_cut" = "Loose" ] && cut_score=100
            
            total_score=$((match_score + cut_score + int_rgi_id))
            
            if [ "$total_score" -gt "$best_score" ]; then
                best_score="$total_score"
                best_entry="$gid|$rgi_gene|$rgi_id|$rgi_cut|$gene_match"
            fi
        done < "${TEMP_ALL}.${gene_id}"
        
        [ -n "$best_entry" ] && echo "$best_entry" >> "$RGI_LOOKUP"
        rm -f "${TEMP_ALL}.${gene_id}"
    done
    
    RGI_COUNT=$(tail -n +2 "$RGI_LOOKUP" | wc -l)
    print_success "RGI processing complete: $RGI_COUNT hits"
    
    echo ""
    echo "📋 RGI Results Sample:"
    head -5 "$RGI_LOOKUP" | while IFS='|' read -r gene_id rgi_gene rgi_id rgi_cut gene_match; do
        if [[ "$gene_id" == "#"* ]]; then
            echo "$gene_id|$rgi_gene|$rgi_id|$rgi_cut|$gene_match"
        else
            printf "%-15s %-60s %-10s %-10s %-10s\n" "$gene_id" "$rgi_gene" "$rgi_id" "$rgi_cut" "$gene_match"
        fi
    done
    
    rm -f "$TEMP_ALL" "$DB_GENE_FILE"
    rm -rf "$RGI_TEMP"
}

# =============================================================================
# STEP 3: RUN ABRICATE
# =============================================================================
run_abricate_and_create_lookup() {
    print_header "RUNNING ABRICATE"
    
    if [ ! -f "$EXTRACTED_FASTA" ] || [ ! -s "$EXTRACTED_FASTA" ]; then
        print_error "FASTA file is empty or missing: $EXTRACTED_FASTA"
        return 1
    fi
    
    ABRICATE_RAW="${ABRICATE_DIR}/${SAMPLE_ID}_abricate_raw.txt"
    
    print_step "Running Abricate..."
    
    docker run --rm \
        --platform linux/amd64 \
        -v "$(dirname "$EXTRACTED_FASTA"):/data" \
        -w "/data" \
        staphb/abricate:latest \
        abricate --db card \
                 --minid 80 \
                 --mincov 50 \
                 "$(basename "$EXTRACTED_FASTA")" > "$ABRICATE_RAW" 2>&1
    
    ABRICATE_LOOKUP="${ABRICATE_DIR}/abricate_lookup.txt"
    echo "# Gene_ID|Abricate_Gene|Abricate_Identity|Abricate_Coverage" > "$ABRICATE_LOOKUP"
    
    print_step "Parsing Abricate results..."
    
    GENE_PREFIX=$(detect_gene_prefix | tail -1)
    
    grep -v "^#" "$ABRICATE_RAW" | while IFS=$'\t' read -r file sequence start end strand gene coverage coverage_map gaps pct_coverage identity database accession product resistance; do
        gene_id=$(echo "$sequence" | grep -o "${GENE_PREFIX}_[0-9]*" | head -1)
        if [ -z "$gene_id" ]; then
            gene_id=$(echo "$sequence" | grep -o '[A-Z]\+_[0-9]\+' | head -1)
        fi
        
        if [ -n "$gene_id" ] && [ -n "$gene" ]; then
            safe_gene=$(echo "$gene" | tr '|' '_')
            
            if grep -q "^${gene_id}|" "$ABRICATE_LOOKUP"; then
                existing_line=$(grep "^${gene_id}|" "$ABRICATE_LOOKUP")
                existing_id=$(echo "$existing_line" | cut -d'|' -f3)
                if [ "$(echo "$identity > $existing_id" | bc -l 2>/dev/null || echo "0")" -eq 1 ]; then
                    grep -v "^${gene_id}|" "$ABRICATE_LOOKUP" > "${ABRICATE_LOOKUP}.tmp"
                    mv "${ABRICATE_LOOKUP}.tmp" "$ABRICATE_LOOKUP"
                    echo "${gene_id}|${safe_gene}|${identity}|${pct_coverage}" >> "$ABRICATE_LOOKUP"
                fi
            else
                echo "${gene_id}|${safe_gene}|${identity}|${pct_coverage}" >> "$ABRICATE_LOOKUP"
            fi
        fi
    done
    
    AB_COUNT=$(tail -n +2 "$ABRICATE_LOOKUP" | wc -l)
    print_success "Abricate processing complete: $AB_COUNT hits"
    
    echo ""
    echo "📋 Abricate Results Sample:"
    head -5 "$ABRICATE_LOOKUP" | while IFS='|' read -r gene_id gene identity coverage; do
        if [[ "$gene_id" == "#"* ]]; then
            echo "$gene_id|$gene|$identity|$coverage"
        else
            printf "%-15s %-40s %-10s %-10s\n" "$gene_id" "$gene" "$identity" "$coverage"
        fi
    done
}

# =============================================================================
# STEP 4: CREATE VALIDATION REPORT (THÊM GENE_MATCH_QC)
# =============================================================================
create_validation_report() {
    print_header "CREATING VALIDATION REPORT"
    
    RGI_LOOKUP="${RGI_DIR}/rgi_lookup.txt"
    ABRICATE_LOOKUP="${ABRICATE_DIR}/abricate_lookup.txt"
    VALIDATION_REPORT="${VALIDATION_DIR}/${SAMPLE_ID}_validation_report.tsv"
    
    # HEADER: THÊM CỘT Gene_Match_QC
    echo -e "Gene_ID\tGene_Name\tSample_ID\tValidation_Date\tOriginal_RGI_Identity\tOriginal_RGI_Cutoff\tValidation_Identity\tValidation_RGI_Cutoff\tGene_Match_QC\tQC_Status\tAbricate_Gene\tAbricate_Identity\tAbricate_Coverage\tRGI_Gene\tRGI_Identity\tRGI_Cutoff\tValidation_Status\tGene_Match_Status\tHas_Nodes" > "$VALIDATION_REPORT"
    
    TOTAL_GENES=0
    BOTH_VALIDATED=0
    ONE_VALIDATED=0
    NONE_VALIDATED=0
    
    print_step "Processing database genes..."
    
    awk -F'\t' -v sample="$SAMPLE_ID" '$3 == sample' "$DB_V9" | while IFS=$'\t' read -r gene_id gene_name sample_id original_contig gene_start gene_end gene_length node_sequence node_count original_rgi_identity validation_identity qc_status source cut_off has_nodes; do
        
        ((TOTAL_GENES++))
        
        # Lấy Original Cutoff
        ORIG_CUT=$(get_original_rgi_cutoff "$sample_id")
        [ -z "$ORIG_CUT" ] && ORIG_CUT="N/A"
        
        # Lấy Gene_Match_QC từ RGI của contig ráp
        GENE_MATCH_QC=$(get_gene_match_qc "$sample_id" "$gene_id" "$gene_name")
        
        # Xác định QC_Status mới
        if [[ "$GENE_MATCH_QC" == DIFFERENT* ]] || [[ "$GENE_MATCH_QC" == "NO_RGI_HIT" ]]; then
            QC_STATUS="FAILED"
        else
            QC_STATUS="$qc_status"
        fi
        
        # Lấy Abricate results
        AB_LINE=$(grep "^${gene_id}|" "$ABRICATE_LOOKUP" 2>/dev/null | head -1)
        if [ -n "$AB_LINE" ]; then
            IFS='|' read -r gid ab_gene ab_id ab_cov <<< "$AB_LINE"
            AB_GENE="$ab_gene"; AB_ID="$ab_id"; AB_COV="$ab_cov"
        else
            AB_GENE="NOT_FOUND"; AB_ID="0"; AB_COV="0"
        fi
        
        # Lấy RGI results
        RGI_LINE=$(grep "^${gene_id}|" "$RGI_LOOKUP" 2>/dev/null | head -1)
        if [ -n "$RGI_LINE" ]; then
            IFS='|' read -r gid rgi_gene rgi_id rgi_cut rgi_match <<< "$RGI_LINE"
            RGI_GENE="$rgi_gene"; RGI_ID="$rgi_id"; RGI_CUT="$rgi_cut"
        else
            RGI_GENE="NOT_FOUND"; RGI_ID="0"; RGI_CUT="N/A"; rgi_match="NO"
        fi
        
        # Validation Status
        VALIDATION_STATUS="NOT_VALIDATED"
        if [ "$AB_GENE" != "NOT_FOUND" ] && [ "$RGI_GENE" != "NOT_FOUND" ]; then
            VALIDATION_STATUS="BOTH_VALIDATED"
            ((BOTH_VALIDATED++))
        elif [ "$AB_GENE" != "NOT_FOUND" ] || [ "$RGI_GENE" != "NOT_FOUND" ]; then
            VALIDATION_STATUS="ONE_VALIDATED"
            ((ONE_VALIDATED++))
        else
            ((NONE_VALIDATED++))
        fi
        
        # Gene Match Status
        GENE_MATCH_STATUS="NO_MATCH"
        ab_match="NO"
        if [ "$AB_GENE" != "NOT_FOUND" ]; then
            ab_gene_clean=$(echo "$AB_GENE" | tr '_' ' ')
            gene_names_match "$gene_name" "$ab_gene_clean" && ab_match="YES"
        fi
        
        if [ -z "$rgi_match" ]; then
            rgi_match="NO"
            if [ "$RGI_GENE" != "NOT_FOUND" ]; then
                gene_names_match "$gene_name" "$RGI_GENE" && rgi_match="YES"
            fi
        fi
        
        if [ "$ab_match" = "YES" ] && [ "$rgi_match" = "YES" ]; then
            GENE_MATCH_STATUS="BOTH_MATCH"
        elif [ "$ab_match" = "YES" ]; then
            GENE_MATCH_STATUS="ABRICATE_MATCH"
        elif [ "$rgi_match" = "YES" ]; then
            GENE_MATCH_STATUS="RGI_MATCH"
        elif [ "$VALIDATION_STATUS" != "NOT_VALIDATED" ]; then
            GENE_MATCH_STATUS="DIFFERENT_GENE"
        fi
        
        VALIDATION_DATE=$(date +"%Y-%m-%d")
        
        # Ghi dòng với cột Gene_Match_QC
        echo -e "$gene_id\t$gene_name\t$sample_id\t$VALIDATION_DATE\t$original_rgi_identity\t$ORIG_CUT\t$validation_identity\t$cut_off\t$GENE_MATCH_QC\t$QC_STATUS\t$AB_GENE\t$AB_ID\t$AB_COV\t$RGI_GENE\t$RGI_ID\t$RGI_CUT\t$VALIDATION_STATUS\t$GENE_MATCH_STATUS\t$has_nodes" >> "$VALIDATION_REPORT"
        
    done
    
    print_success "Validation report created: $TOTAL_GENES genes"
    echo "  - BOTH_VALIDATED: $BOTH_VALIDATED"
    echo "  - ONE_VALIDATED: $ONE_VALIDATED"
    echo "  - NONE_VALIDATED: $NONE_VALIDATED"
    
    echo ""
    echo "📊 Genes found by RGI: $(tail -n +2 "$RGI_LOOKUP" | wc -l 2>/dev/null || echo "0")"
    echo "📊 Genes found by Abricate: $(tail -n +2 "$ABRICATE_LOOKUP" | wc -l 2>/dev/null || echo "0")"
    
    echo ""
    echo "📋 Validation Report Sample (First 5 genes):"
    if [ -f "$VALIDATION_REPORT" ] && [ $(wc -l < "$VALIDATION_REPORT" 2>/dev/null) -gt 1 ]; then
        head -6 "$VALIDATION_REPORT" | column -t -s $'\t'
    else
        echo "  No data in validation report"
    fi
}

# =============================================================================
# STEP 5: UPDATE MASTER DATABASE
# =============================================================================
update_master_database() {
    print_header "UPDATING MASTER DATABASE"
    
    VALIDATION_REPORT="${VALIDATION_DIR}/${SAMPLE_ID}_validation_report.tsv"
    
    if [ ! -f "$VALIDATION_REPORT" ]; then
        print_error "Validation report not found"
        return 1
    fi
    
    VALIDATED_COUNT=$(tail -n +2 "$VALIDATION_REPORT" | awk -F'\t' '$14 != "NOT_VALIDATED" {count++} END {print count}')
    
    if [ "$VALIDATED_COUNT" -eq 0 ]; then
        print_warning "No genes validated! Skipping master database update."
        return 0
    fi
    
    if [ ! -f "$MASTER_VALIDATION_DB" ]; then
        print_step "Creating new master database..."
        mkdir -p "$(dirname "$MASTER_VALIDATION_DB")"
        head -1 "$VALIDATION_REPORT" > "$MASTER_VALIDATION_DB"
    fi
    
    print_step "Removing old entries for $SAMPLE_ID..."
    if [ -f "$MASTER_VALIDATION_DB" ]; then
        grep -v "$SAMPLE_ID" "$MASTER_VALIDATION_DB" > "${MASTER_VALIDATION_DB}.tmp" 2>/dev/null || true
        if [ -s "${MASTER_VALIDATION_DB}.tmp" ]; then
            mv "${MASTER_VALIDATION_DB}.tmp" "$MASTER_VALIDATION_DB"
        else
            head -1 "$VALIDATION_REPORT" > "$MASTER_VALIDATION_DB"
        fi
    fi
    
    print_step "Adding new validation results..."
    tail -n +2 "$VALIDATION_REPORT" >> "$MASTER_VALIDATION_DB"
    
    # Xóa dòng trống
    sed -i '' '/^[[:space:]]*$/d' "$MASTER_VALIDATION_DB"
    
    MASTER_COUNT=$(tail -n +2 "$MASTER_VALIDATION_DB" | wc -l 2>/dev/null || echo "0")
    MASTER_SAMPLES=$(tail -n +2 "$MASTER_VALIDATION_DB" | cut -f3 | sort -u | wc -l 2>/dev/null || echo "0")
    
    print_success "Master database updated"
    echo "  Total entries: $MASTER_COUNT"
    echo "  Unique samples: $MASTER_SAMPLES"
    echo "  File: $MASTER_VALIDATION_DB"
}

# =============================================================================
# STEP 6: CREATE STATISTICS
# =============================================================================
create_statistics() {
    print_header "CREATING STATISTICS"
    
    VALIDATION_REPORT="${VALIDATION_DIR}/${SAMPLE_ID}_validation_report.tsv"
    STATS_FILE="${VALIDATION_DIR}/${SAMPLE_ID}_validation_statistics.txt"
    
    if [ ! -f "$VALIDATION_REPORT" ] || [ $(wc -l < "$VALIDATION_REPORT" 2>/dev/null) -le 1 ]; then
        print_error "Validation report is empty or missing"
        return 1
    fi
    
    TOTAL=$(tail -n +2 "$VALIDATION_REPORT" | wc -l)
    BOTH_V=$(grep -c "BOTH_VALIDATED" "$VALIDATION_REPORT")
    BOTH_M=$(grep -c "BOTH_MATCH" "$VALIDATION_REPORT")
    DIFF_GENES=$(grep -c "DIFFERENT_GENE" "$VALIDATION_REPORT")
    QC_PASSED=$(grep -c "PASSED" "$VALIDATION_REPORT")
    
    cat > "$STATS_FILE" << EOF
VALIDATION STATISTICS - $SAMPLE_ID
====================================
Date: $(date)
Sample: $SAMPLE_ID

OVERALL RESULTS:
---------------
Total Genes: $TOTAL
Genes with Nodes: $TOTAL

VALIDATION SUCCESS:
------------------
Both Tools Validated: $BOTH_V
Gene Name Matches: $BOTH_M
Different Genes Detected: $DIFF_GENES
QC Passed: $QC_PASSED

VALIDATION STATUS:
------------------
$(tail -n +2 "$VALIDATION_REPORT" | cut -f16 | sort | uniq -c | awk '{printf "  %-20s: %d genes\\n", \$2, \$1}')

GENE MATCH STATUS:
------------------
$(tail -n +2 "$VALIDATION_REPORT" | cut -f17 | sort | uniq -c | awk '{printf "  %-20s: %d genes\\n", \$2, \$1}')

QC STATUS:
----------
$(tail -n +2 "$VALIDATION_REPORT" | cut -f10 | sort | uniq -c | awk '{printf "  %-20s: %d genes\\n", \$2, \$1}')

FILES:
------
Validation Report: $(basename "$VALIDATION_REPORT")
Master Database: $(basename "$MASTER_VALIDATION_DB")
This Report: $(basename "$STATS_FILE")
EOF
    
    if [ "$DIFF_GENES" -gt 0 ]; then
        echo "" >> "$STATS_FILE"
        echo "DIFFERENT GENES DETECTED:" >> "$STATS_FILE"
        echo "-------------------------" >> "$STATS_FILE"
        grep "DIFFERENT_GENE" "$VALIDATION_REPORT" | head -10 | awk -F'\t' '{printf "  %s: DB=%s → RGI=%s (Identity: %s, Cutoff: %s)\\n", \$1, \$2, \$13, \$14, \$15}' >> "$STATS_FILE"
    fi
    
    print_success "Statistics created: $STATS_FILE"
    
    echo ""
    echo "📊 SUMMARY:"
    echo "  Total Genes: $TOTAL"
    echo "  Both Validated: $BOTH_V"
    echo "  Both Match: $BOTH_M"
    echo "  Different Genes: $DIFF_GENES"
    echo "  QC Passed: $QC_PASSED"
}

# =============================================================================
# MAIN
# =============================================================================
main() {
    print_header "AMR VALIDATION PIPELINE v3 - WITH CUTOFF COLUMNS & GENE_MATCH_QC"
    echo "Sample: $SAMPLE_ID"
    echo "Timestamp: $TIMESTAMP"
    echo "Work Directory: $CURRENT_WORK_DIR"
    echo ""
    
    setup_environment
    run_rgi_and_create_lookup
    run_abricate_and_create_lookup
    create_validation_report
    update_master_database
    create_statistics
    
    print_header "✅ VALIDATION PIPELINE COMPLETED"
    echo ""
    echo "📁 OUTPUT FILES:"
    echo "  • Validation Report: ${VALIDATION_DIR}/${SAMPLE_ID}_validation_report.tsv"
    echo "  • RGI Results: ${RGI_DIR}/rgi_lookup.txt"
    echo "  • Abricate Results: ${ABRICATE_DIR}/abricate_lookup.txt"
    echo "  • Statistics: ${VALIDATION_DIR}/${SAMPLE_ID}_validation_statistics.txt"
    echo "  • Master Database: $MASTER_VALIDATION_DB"
}

main