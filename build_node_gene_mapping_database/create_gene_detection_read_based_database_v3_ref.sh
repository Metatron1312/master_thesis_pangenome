#!/bin/bash
# =============================================================================
# PIPELINE REFERENCE BJAB07104 - HOÀN CHỈNH (GIỮ NGUYÊN SCRIPT 5 GỐC)
# =============================================================================

set -e
set -u

SAMPLE_ID="BJAB07104"
BASE_DIR="/Users/shinra"
SCRIPT_DIR="${BASE_DIR}/gene_detection_read_based_script"

# Tạo thư mục
SAMPLE_DIR="${BASE_DIR}/${SAMPLE_ID}"
mkdir -p "$SAMPLE_DIR"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="${SAMPLE_DIR}/pipeline_${TIMESTAMP}.log"
ERROR_LOG="${SAMPLE_DIR}/pipeline_errors_${TIMESTAMP}.log"

# Scripts gốc (GIỮ NGUYÊN)
SCRIPT1="${SCRIPT_DIR}/map_read_to_contig_AMR_only_v4.py"
SCRIPT2="${SCRIPT_DIR}/shovill_arm_gene_reads.sh"
SCRIPT3="${SCRIPT_DIR}/vg_giraffe_amr_gene_contigs.py"
SCRIPT4_1="${SCRIPT_DIR}/rgi_amr_gene_contig.sh"
SCRIPT4_2="${SCRIPT_DIR}/abricate.sh"
SCRIPT5="${SCRIPT_DIR}/extract_node_path_of_arm_gene_v12.py"
SCRIPT6="${SCRIPT_DIR}/extract_sequence_from_node_gene_database_v3.sh"
SCRIPT7="${SCRIPT_DIR}/arm_validation_pipeline_v3.sh"

# Database files
BAKTA_DB="${BASE_DIR}/Bioinformatics/WGS_AMR/data_public/Graph_annotation/bakta_rgi_integration_107_samples_v5.tsv"
REFERENCE_CONTIG="${BASE_DIR}/Bioinformatics/WGS_AMR/data_public/combined_contigs_and_refseq_v2/BJAB07104.filtered.complete.fasta"
NON_AMR_DB="${BASE_DIR}/Bioinformatics/WGS_AMR/data_public/Graph_annotation/non_amr_reference_database_v1.tsv"

# Graph files cho NON-AMR
GRAPH_DIR="${BASE_DIR}/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/vg_graph/graph/106_strains_with_refseq_complete_vg_autoindex"
GRAPH_GBZ="${GRAPH_DIR}/complete_106_strains_with_refseq_pangenome.giraffe.gbz"
MIN_INDEX="${GRAPH_DIR}/complete_106_strains_with_refseq_pangenome.shortread.withzip.min"
DIST_INDEX="${GRAPH_DIR}/complete_106_strains_with_refseq_pangenome.dist"

# Đường dẫn dữ liệu
ASSEMBLY_DIR="${SAMPLE_DIR}/amr_gene_assemblies"
MAPPED_DIR="${SAMPLE_DIR}/arm_contigs_mapped_to_graph"
RGI_RESULTS_DIR="${SAMPLE_DIR}/amr_gene_contig_rgi_results"
ABRICATE_RESULTS_DIR="${SAMPLE_DIR}/amr_gene_contig_abricate_results"
SEQUENCE_EXTRACTION_DIR="${BASE_DIR}/sequence_extraction/${SAMPLE_ID}"
NON_AMR_DIR="${SAMPLE_DIR}/non_amr_data"

mkdir -p "$ASSEMBLY_DIR" "$MAPPED_DIR" "$RGI_RESULTS_DIR" "$ABRICATE_RESULTS_DIR"
mkdir -p "$SEQUENCE_EXTRACTION_DIR" "$NON_AMR_DIR"

# ==================== FUNCTIONS ====================
log() { echo -e "$1" 2>&1 | tee -a "$LOG_FILE"; }
log_error() { echo -e "❌ ERROR: $1" 2>&1 | tee -a "$LOG_FILE" "$ERROR_LOG"; }
log_success() { echo -e "✅ $1" 2>&1 | tee -a "$LOG_FILE"; }
log_step() { 
    echo "" 2>&1 | tee -a "$LOG_FILE"
    echo "=========================================" 2>&1 | tee -a "$LOG_FILE"
    echo "🚀 $1" 2>&1 | tee -a "$LOG_FILE"
    echo "=========================================" 2>&1 | tee -a "$LOG_FILE"
}

# ==================== THAY THẾ STEP 1-2 ====================
setup_reference_assembly() {
    log_step "THAY THẾ STEP 1-2: TẠO FASTA TỪ BAKTA DATABASE"
    
    python3 << EOF
import csv, os, re

print("📖 Đọc reference contig...")
contig_seq = {}
with open("$REFERENCE_CONTIG", 'r') as f:
    header = None
    seq = []
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if header:
                contig_seq[header] = ''.join(seq)
            header = line[1:].split()[0]
            seq = []
        else:
            seq.append(line)
    if header:
        contig_seq[header] = ''.join(seq)

print("🔍 Lấy AMR genes từ bakta database...")
amr_count = 0

with open("$BAKTA_DB", 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if row.get('Sample_ID') != "$SAMPLE_ID":
            continue
        gene_id = row.get('GeneID', '')
        amr_val = row.get('AMR_Gene', '').strip()
        
        # AMR nếu có ARO number
        if not (amr_val and re.match(r'^\d+$', amr_val)):
            continue
        
        start = int(row.get('Start', 0))
        end = int(row.get('End', 0))
        strand = row.get('Strand', '+')
        contig_name = row.get('Original_Contig', 'contig_1')
        product = row.get('Product', '')
        gene_name = row.get('Gene_Name', '') or product[:50]
        
        if contig_name in contig_seq and start > 0 and end > 0:
            seq = contig_seq[contig_name][start-1:end]
            if strand == '-':
                comp = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
                seq = ''.join(comp.get(b, b) for b in seq[::-1])
            
            # Tạo thư mục cho gene này (giống step 2 output)
            gene_dir = os.path.join("$ASSEMBLY_DIR", gene_id)
            os.makedirs(gene_dir, exist_ok=True)
            with open(os.path.join(gene_dir, f"{gene_id}_contigs.fa"), 'w') as f:
                f.write(f">{gene_id} {gene_name}\n")
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + "\n")
            amr_count += 1

print(f"✅ Đã tạo FASTA cho {amr_count} AMR genes")
EOF

    local cnt=$(find "$ASSEMBLY_DIR" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l | tr -d ' ')
    log_success "Đã tạo FASTA cho $cnt AMR genes"
    echo "$cnt" > "$SAMPLE_DIR/amr_count.txt"
}

# ==================== HÀM CHẠY SCRIPT (GIỮ NGUYÊN) ====================
run_script() {
    local script=$1
    local step_name=$2
    local need_sample=$3
    
    log_step "$step_name"
    log "Script: $(basename $script)"
    
    local tmp_script="/tmp/$(basename $script).$$.tmp"
    cp "$script" "$tmp_script"
    
    # Thay thế sample ID
    local old_sample=$(grep -Eo 'SAMEA[0-9]{7}|SAMN[0-9]{8}' "$script" 2>/dev/null | head -1)
    if [ -n "$old_sample" ] && [ "$old_sample" != "$SAMPLE_ID" ]; then
        log "🔄 Thay thế: $old_sample → $SAMPLE_ID"
        sed -i '' "s/${old_sample}/${SAMPLE_ID}/g" "$tmp_script"
        sed -i '' "s|/${old_sample}/|/${SAMPLE_ID}/|g" "$tmp_script"
        sed -i '' "s/SAMPLE_ID=\"${old_sample}\"/SAMPLE_ID=\"${SAMPLE_ID}\"/g" "$tmp_script"
    fi
    
    chmod +x "$tmp_script"
    cd "$SAMPLE_DIR"
    
    if [[ "$script" == *.py ]]; then
        python3 "$tmp_script" 2>&1 | tee -a "$LOG_FILE"
    else
        if [ "$need_sample" == "yes" ]; then
            bash "$tmp_script" "$SAMPLE_ID" 2>&1 | tee -a "$LOG_FILE"
        else
            bash "$tmp_script" 2>&1 | tee -a "$LOG_FILE"
        fi
    fi
    
    rm -f "$tmp_script"
    log_success "$step_name hoàn thành"
}

# ==================== STEP 8: TẠO NON-AMR DATABASE ====================
step8_non_amr() {
    log_step "STEP 8: TẠO NON-AMR REFERENCE DATABASE"
    
    log "📋 Đang tạo NON-AMR genes từ bakta database..."
    
    # Tạo FASTA cho NON-AMR genes
    python3 << EOF
import csv, os, re

print("📖 Đọc reference contig...")
contig_seq = {}
with open("$REFERENCE_CONTIG", 'r') as f:
    header = None
    seq = []
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if header:
                contig_seq[header] = ''.join(seq)
            header = line[1:].split()[0]
            seq = []
        else:
            seq.append(line)
    if header:
        contig_seq[header] = ''.join(seq)

print("🔍 Đọc bakta database và lấy NON-AMR genes...")
non_amr_genes = []
amr_set = set()

# Đọc danh sách AMR genes
with open("$BAKTA_DB", 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if row.get('Sample_ID') == "$SAMPLE_ID":
            amr_val = row.get('AMR_Gene', '').strip()
            if amr_val and re.match(r'^\d+$', amr_val):
                amr_set.add(row.get('GeneID', ''))

# Đọc NON-AMR genes
with open("$BAKTA_DB", 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if row.get('Sample_ID') != "$SAMPLE_ID":
            continue
        gene_id = row.get('GeneID', '')
        if not gene_id or gene_id in amr_set:
            continue
        
        start = int(row.get('Start', 0))
        end = int(row.get('End', 0))
        strand = row.get('Strand', '+')
        contig = row.get('Original_Contig', 'contig_1')
        name = row.get('Gene_Name', '') or row.get('Product', '')[:50]
        
        if contig in contig_seq and start > 0 and end > 0:
            seq = contig_seq[contig][start-1:end]
            if strand == '-':
                comp = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
                seq = ''.join(comp.get(b, b) for b in seq[::-1])
            non_amr_genes.append((gene_id, name, contig, start, end, len(seq), seq))

print(f"📊 Tìm thấy {len(non_amr_genes)} NON-AMR genes")

# Tạo FASTA cho NON-AMR
with open("$NON_AMR_DIR/non_amr.fasta", 'w') as f:
    for gid, name, contig, start, end, length, seq in non_amr_genes:
        f.write(f">{gid}|{name}|{contig}|{start}|{end}|{length}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")

# Lưu metadata
with open("$NON_AMR_DIR/non_amr_metadata.tsv", 'w') as f:
    f.write("gene_id\tgene_name\tcontig\tstart\tend\tlength\n")
    for gid, name, contig, start, end, length, _ in non_amr_genes:
        f.write(f"{gid}\t{name}\t{contig}\t{start}\t{end}\t{length}\n")

print(f"✅ Đã tạo FASTA cho {len(non_amr_genes)} NON-AMR genes")
EOF

    # Map NON-AMR genes
    log "🗺️ Đang map NON-AMR genes vào graph..."
    
    local fasta="$NON_AMR_DIR/non_amr.fasta"
    local fastq="$NON_AMR_DIR/temp.fastq"
    local gam="$NON_AMR_DIR/non_amr.gam"
    
    python3 << EOF
with open("$fasta", 'r') as f, open("$fastq", 'w') as out:
    lines = f.readlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith('>'):
            header = lines[i][1:].split('|')[0]
            i += 1
            seq = ""
            while i < len(lines) and not lines[i].startswith('>'):
                seq += lines[i].strip()
                i += 1
            out.write(f"@{header}\n{seq}\n+\n{'F'*len(seq)}\n")
        else:
            i += 1
print("✅ Chuyển sang FASTQ")
EOF
    
    vg giraffe -Z "$GRAPH_GBZ" -m "$MIN_INDEX" -d "$DIST_INDEX" -f "$fastq" -t 8 > "$gam" 2>> "$LOG_FILE"
    rm -f "$fastq"
    
    # Extract nodes và tạo database
    log "📁 Đang tạo NON-AMR database..."
    
    python3 << EOF
import csv, os, json, subprocess, re
from datetime import datetime

metadata = {}
with open("$NON_AMR_DIR/non_amr_metadata.tsv", 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        metadata[row['gene_id']] = row

gene_nodes = {}
gam_file = "$gam"
if os.path.exists(gam_file) and os.path.getsize(gam_file) > 0:
    result = subprocess.run(f"vg view -a '{gam_file}' -j 2>/dev/null", shell=True, capture_output=True, text=True, timeout=120)
    if result.stdout:
        for line in result.stdout.split('\n'):
            if not line.strip():
                continue
            try:
                data = json.loads(line)
                if 'name' in data and 'path' in data:
                    gid = data['name']
                    nodes = []
                    for m in data['path'].get('mapping', []):
                        if 'position' in m and 'node_id' in m['position']:
                            nodes.append(str(m['position']['node_id']))
                    if nodes:
                        gene_nodes[gid] = list(set(nodes))
            except:
                pass

os.makedirs(os.path.dirname("$NON_AMR_DB"), exist_ok=True)
with open("$NON_AMR_DB", 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['gene_id', 'gene_name', 'sample_id', 'original_contig', 
                                           'gene_start', 'gene_end', 'gene_length', 'node_sequence', 
                                           'node_count', 'extraction_date'], delimiter='\t')
    writer.writeheader()
    written = 0
    for gid, meta in metadata.items():
        nodes = gene_nodes.get(gid, [])
        writer.writerow({
            'gene_id': gid,
            'gene_name': meta['gene_name'],
            'sample_id': "$SAMPLE_ID",
            'original_contig': meta['contig'],
            'gene_start': meta['start'],
            'gene_end': meta['end'],
            'gene_length': meta['length'],
            'node_sequence': ','.join(nodes) if nodes else 'N/A',
            'node_count': len(nodes),
            'extraction_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        })
        written += 1

print(f"✅ Đã ghi {written} NON-AMR genes vào database")
EOF
    
    log_success "NON-AMR database đã tạo xong"
}

# ==================== MAIN PIPELINE ====================
main() {
    echo "========================================="
    echo "🚀 AMR GENE PIPELINE - REFERENCE VERSION"
    echo "========================================="
    echo "Sample ID: $SAMPLE_ID"
    echo "Log file: $LOG_FILE"
    echo "========================================="
    echo ""
    echo "⚠️  LƯU Ý: Step 1-2 được thay thế bằng việc tạo FASTA từ bakta database"
    echo "   Các step 3-7 chạy HOÀN TOÀN GIỐNG pipeline AMR gốc"
    echo "   Khi được hỏi 'Continue? (y/N):' -> gõ y và Enter"
    echo "   Step 8: Tạo NON-AMR database tham khảo"
    echo ""
    
    # Thay thế step 1-2
    setup_reference_assembly
    
    # ==================== STEP 3 ====================
    run_script "$SCRIPT3" "STEP 3: MAP AMR CONTIGS TO GRAPH" "no"
    
    # ==================== STEP 4.1 ====================
    run_script "$SCRIPT4_1" "STEP 4.1: RGI ANALYSIS" "no"
    
    # ==================== STEP 4.2 ====================
    run_script "$SCRIPT4_2" "STEP 4.2: ABRICATE ANALYSIS" "no"
    
    # ==================== STEP 5 ====================
    run_script "$SCRIPT5" "STEP 5: UPDATE AMR DATABASE" "no"
    
    # ==================== STEP 6 ====================
    run_script "$SCRIPT6" "STEP 6: EXTRACT SEQUENCES" "yes"
    
    # ==================== STEP 7 ====================
    run_script "$SCRIPT7" "STEP 7: VALIDATION" "no"
    
    # ==================== STEP 8: NON-AMR DATABASE ====================
    step8_non_amr
    
    # ==================== COMPLETION ====================
    echo ""
    echo "========================================="
    echo "🎉 PIPELINE COMPLETED!"
    echo "========================================="
    echo "Sample: $SAMPLE_ID"
    echo ""
    echo "📊 KẾT QUẢ:"
    echo "   🔴 AMR Database: $NODE_DB"
    echo "      (Đã qua step 3-7 như pipeline gốc)"
    echo ""
    echo "   ⚪ NON-AMR Database: $NON_AMR_DB"
    echo "      (Tham khảo - đã map vào graph)"
    echo ""
    echo "📄 Log file: $LOG_FILE"
    echo "========================================="
}

# Run
main