#!/bin/bash

# Script: process_vcf_files.sh
# Description: Nén, index, sửa header và tạo manifest cho VCF files
# Usage: ./process_vcf_files.sh

# ===== CẤU HÌNH =====
VCF_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/vg_graph/msa2vcf_v19"
OUTPUT_DIR="$VCF_DIR/processed"
MANIFEST_FILE="$OUTPUT_DIR/vcf_manifest.txt"
BGZIP_THREADS=8

# ===== KIỂM TRA CÔNG CỤ =====
check_tools() {
    echo "🔧 Kiểm tra công cụ..."
    
    if ! command -v bgzip &> /dev/null; then
        echo "❌ bgzip không được cài đặt. Cài đặt: conda install -c bioconda htslib"
        exit 1
    fi
    
    if ! command -v tabix &> /dev/null; then
        echo "❌ tabix không được cài đặt. Cài đặt: conda install -c bioconda htslib"
        exit 1
    fi
    
    if ! command -v bcftools &> /dev/null; then
        echo "❌ bcftools không được cài đặt. Cài đặt: conda install -c bioconda bcftools"
        exit 1
    fi
    
    echo "✅ Tất cả công cụ đã sẵn sàng"
}

# ===== TẠO THƯ MỤC =====
create_directories() {
    echo "📁 Tạo thư mục output..."
    mkdir -p "$OUTPUT_DIR"
    echo "✅ Output directory: $OUTPUT_DIR"
}

# ===== KIỂM TRA VÀ LỌC VCF FILES =====
filter_vcf_files() {
    echo "🔍 Kiểm tra và lọc VCF files..."
    
    local total_files=$(find "$VCF_DIR" -maxdepth 1 -name "block_*.vcf" | wc -l)
    local valid_files=0
    local empty_files=0
    local no_sample_files=0
    local temp_dir="$OUTPUT_DIR/temp_valid"
    
    mkdir -p "$temp_dir"
    
    for vcf_file in "$VCF_DIR"/block_*.vcf; do
        if [[ -f "$vcf_file" ]]; then
            local base_name=$(basename "$vcf_file" .vcf)
            local is_valid=0
            
            # Kiểm tra 1: File có chứa variants không (có dòng không phải comment)
            local variant_count=$(grep -vc "^#" "$vcf_file" 2>/dev/null || echo 0)
            
            # Kiểm tra 2: File có samples không
            local sample_count=0
            if grep -q "^#CHROM" "$vcf_file"; then
                local chrom_line=$(grep "^#CHROM" "$vcf_file")
                IFS=$'\t' read -ra headers <<< "$chrom_line"
                sample_count=$(( ${#headers[@]} - 9 ))
            fi
            
            # QUYẾT ĐỊNH: File được coi là valid nếu:
            # 1. Có ít nhất 1 variant VÀ ít nhất 1 sample
            if [[ $variant_count -gt 0 && $sample_count -gt 0 ]]; then
                # File có variants và samples -> VALID
                cp "$vcf_file" "$temp_dir/${base_name}.vcf"
                ((valid_files++))
                is_valid=1
                echo "✅ VALID: $base_name (variants: $variant_count, samples: $sample_count)"
            elif [[ $variant_count -eq 0 ]]; then
                # File không có variants -> EMPTY (bỏ qua và XOÁ file gốc)
                echo "🗑️  EMPTY: $base_name - Xoá file gốc"
                rm "$vcf_file"
                ((empty_files++))
            else
                # File có variants nhưng không có samples -> INVALID (bỏ qua)
                echo "⚠️  NO_SAMPLES: $base_name (variants: $variant_count) - Bỏ qua"
                ((no_sample_files++))
            fi
            
            # Hiển thị progress
            local processed=$((valid_files + empty_files + no_sample_files))
            if [[ $((processed % 50)) -eq 0 ]]; then
                echo "📊 Đã xử lý: $processed/$total_files files"
            fi
        fi
    done
    
    echo "=========================================="
    echo "📊 KẾT QUẢ LỌC:"
    echo "   ✅ Valid files: $valid_files"
    echo "   🗑️  Empty files (đã xoá): $empty_files"
    echo "   ⚠️  No-sample files: $no_sample_files"
    echo "   📁 Total processed: $total_files"
    echo "=========================================="
    
    # Nếu không có file valid nào, dừng script
    if [[ $valid_files -eq 0 ]]; then
        echo "❌ Không có VCF files nào valid. Dừng script."
        rm -rf "$temp_dir"
        exit 1
    fi
    
    return $valid_files
}

# ===== SỬA HEADER VCF =====
fix_vcf_headers() {
    echo "🔧 Sửa header VCF files để tương thích với vg..."
    
    local total_files=$(find "$OUTPUT_DIR/temp_valid" -maxdepth 1 -name "block_*.vcf" | wc -l)
    local processed=0
    local temp_dir="$OUTPUT_DIR/temp_fixed"
    
    mkdir -p "$temp_dir"
    
    for vcf_file in "$OUTPUT_DIR/temp_valid"/block_*.vcf; do
        if [[ -f "$vcf_file" ]]; then
            local base_name=$(basename "$vcf_file" .vcf)
            local fixed_file="$temp_dir/${base_name}.vcf"
            
            # Sửa header VCF
            fix_single_vcf_header "$vcf_file" "$fixed_file"
            
            ((processed++))
            
            # Hiển thị progress
            if [[ $((processed % 50)) -eq 0 ]]; then
                echo "📊 Đã sửa header: $processed/$total_files files"
            fi
        fi
    done
    
    # Xóa thư mục temp_valid
    rm -rf "$OUTPUT_DIR/temp_valid"
    
    echo "✅ Đã sửa header cho $processed VCF files"
}

# ===== SỬA HEADER CHO TỪNG VCF =====
fix_single_vcf_header() {
    local input_file="$1"
    local output_file="$2"
    local temp_file="${output_file}.tmp"
    
    # Extract block ID từ filename
    local block_id=$(basename "$input_file" .vcf | sed 's/block_//')
    
    # Tạo file tạm
    > "$temp_file"
    
    # Đọc và sửa header
    while IFS= read -r line; do
        if [[ "$line" == \#\#contig* ]]; then
            # Sửa thành format đơn giản: ##contig=<ID=block_id,length=...>
            if [[ "$line" =~ length=([0-9]+) ]]; then
                local length="${BASH_REMATCH[1]}"
                echo "##contig=<ID=$block_id,length=$length>" >> "$temp_file"
            else
                echo "##contig=<ID=$block_id>" >> "$temp_file"
            fi
        elif [[ "$line" == \#\#fileformat* ]]; then
            echo "##fileformat=VCFv4.2" >> "$temp_file"
        elif [[ "$line" == \#CHROM* ]]; then
            # Đảm bảo header line không có khoảng trắng thừa
            echo "$line" | sed 's/[[:space:]]*$//' >> "$temp_file"
        else
            # Giữ nguyên các dòng header khác
            echo "$line" >> "$temp_file"
        fi
    done < "$input_file"
    
    # Di chuyển file tạm thành output
    mv "$temp_file" "$output_file"
}

# ===== CHUẨN HÓA SAMPLE NAMES =====
normalize_sample_names() {
    echo "👥 Chuẩn hóa sample names trong VCF files..."
    
    local processed=0
    local total_files=$(find "$OUTPUT_DIR/temp_fixed" -maxdepth 1 -name "block_*.vcf" | wc -l)
    
    for vcf_file in "$OUTPUT_DIR/temp_fixed"/block_*.vcf; do
        if [[ -f "$vcf_file" ]]; then
            local base_name=$(basename "$vcf_file" .vcf)
            local normalized_file="$OUTPUT_DIR/${base_name}.vcf"
            
            # Chuẩn hóa sample names
            normalize_vcf_samples "$vcf_file" "$normalized_file"
            
            ((processed++))
            
            if [[ $((processed % 50)) -eq 0 ]]; then
                echo "📊 Đã chuẩn hóa: $processed/$total_files files"
            fi
        fi
    done
    
    # Xóa thư mục tạm
    rm -rf "$OUTPUT_DIR/temp_fixed"
    
    echo "✅ Đã chuẩn hóa sample names cho $processed files"
}

# ===== CHUẨN HÓA SAMPLE NAMES TRONG VCF =====
normalize_vcf_samples() {
    local input_file="$1"
    local output_file="$2"
    local temp_file="${output_file}.tmp"
    
    # Tạo file tạm
    cp "$input_file" "$temp_file"
    
    # Tìm dòng #CHROM (header cuối)
    local header_line=$(grep "^#CHROM" "$temp_file")
    
    if [[ -n "$header_line" ]]; then
        # Tách các cột
        IFS=$'\t' read -r -a columns <<< "$header_line"
        
        # Chuẩn hóa sample names (cột từ 9 trở đi)
        for ((i=9; i<${#columns[@]}; i++)); do
            # Thay thế khoảng trắng và ký tự đặc biệt
            columns[i]=$(echo "${columns[i]}" | sed -E 's/[^a-zA-Z0-9_.-]/_/g')
            
            # Đảm bảo không rỗng
            if [[ -z "${columns[i]}" ]]; then
                columns[i]="sample_$((i-8))"
            fi
        done
        
        # Tạo header line mới
        local new_header=$(IFS=$'\t'; echo "${columns[*]}")
        
        # Thay thế header line trong file
        sed -i.tmp "s/^#CHROM.*/$new_header/" "$temp_file"
        rm -f "${temp_file}.tmp"
    fi
    
    # Di chuyển file tạm thành output
    mv "$temp_file" "$output_file"
}

# ===== KIỂM TRA VCF FILES TRƯỚC KHI NÉN =====
validate_vcf_files() {
    echo "🔍 Kiểm tra VCF files trước khi nén..."
    
    local total_files=$(find "$OUTPUT_DIR" -maxdepth 1 -name "block_*.vcf" | wc -l)
    local valid_files=0
    local invalid_files=0
    
    for vcf_file in "$OUTPUT_DIR"/block_*.vcf; do
        if [[ -f "$vcf_file" ]]; then
            # Kiểm tra với bcftools
            if bcftools view -h "$vcf_file" &> /dev/null; then
                # Kiểm tra có variants và samples không
                local variant_count=$(grep -vc "^#" "$vcf_file")
                local sample_count=$(bcftools view -h "$vcf_file" | grep "^#CHROM" | head -1 | awk '{print NF-9}')
                
                if [[ $variant_count -gt 0 && $sample_count -gt 0 ]]; then
                    ((valid_files++))
                else
                    echo "❌ File không đủ điều kiện: $(basename "$vcf_file") (variants: $variant_count, samples: $sample_count)"
                    rm "$vcf_file"
                    ((invalid_files++))
                fi
            else
                echo "❌ File lỗi: $(basename "$vcf_file")"
                rm "$vcf_file"
                ((invalid_files++))
            fi
        fi
    done
    
    echo "✅ Kiểm tra hoàn tất: $valid_files files hợp lệ, $invalid_files files bị loại bỏ"
}

# ===== NÉN VCF FILES =====
compress_vcf_files() {
    echo "🗜️  Nén VCF files..."
    
    local total_files=$(find "$OUTPUT_DIR" -maxdepth 1 -name "block_*.vcf" | wc -l)
    local processed=0
    
    for vcf_file in "$OUTPUT_DIR"/block_*.vcf; do
        if [[ -f "$vcf_file" ]]; then
            local base_name=$(basename "$vcf_file" .vcf)
            local compressed_file="$OUTPUT_DIR/${base_name}.vcf.gz"
            
            # Nén file VCF
            bgzip -c -@ $BGZIP_THREADS "$vcf_file" > "$compressed_file"
            
            # Index file đã nén
            tabix -p vcf "$compressed_file"
            
            # Xóa file VCF gốc đã sửa
            rm "$vcf_file"
            
            ((processed++))
            
            # Hiển thị progress
            if [[ $((processed % 100)) -eq 0 ]]; then
                echo "📊 Đã nén: $processed/$total_files files"
            fi
        fi
    done
    
    echo "✅ Đã nén và index $processed VCF files"
}

# ===== KIỂM TRA HEADER VCF ĐÃ NÉN =====
verify_vcf_headers() {
    echo "🔍 Kiểm tra header VCF files đã nén..."
    
    local checked=0
    local errors=0
    
    for vcf_file in "$OUTPUT_DIR"/block_*.vcf.gz; do
        if [[ -f "$vcf_file" ]]; then
            # Kiểm tra header với bcftools
            if bcftools view -h "$vcf_file" &> /dev/null; then
                # Kiểm tra số samples và variants
                local sample_count=$(bcftools view -h "$vcf_file" | grep "^#CHROM" | head -1 | awk '{print NF-9}')
                local variant_count=$(bcftools view "$vcf_file" | grep -vc "^#")
                
                if [[ $sample_count -gt 0 && $variant_count -gt 0 ]]; then
                    ((checked++))
                else
                    echo "❌ File không đủ điều kiện: $(basename "$vcf_file") (variants: $variant_count, samples: $sample_count)"
                    ((errors++))
                fi
            else
                echo "❌ File lỗi: $(basename "$vcf_file")"
                ((errors++))
            fi
            
            if [[ $((checked % 100)) -eq 0 ]]; then
                echo "📊 Đã kiểm tra: $checked files"
            fi
        fi
    done
    
    echo "✅ Đã kiểm tra $checked VCF files, $errors files lỗi"
}

# ===== TẠO MANIFEST FILE =====
create_manifest() {
    echo "📝 Tạo manifest file..."
    
    # Xóa manifest file cũ nếu tồn tại
    rm -f "$MANIFEST_FILE"
    
    # Tạo manifest với đường dẫn đầy đủ đến các file .vcf.gz
    find "$OUTPUT_DIR" -name "block_*.vcf.gz" | sort > "$MANIFEST_FILE"
    
    local file_count=$(wc -l < "$MANIFEST_FILE")
    echo "✅ Manifest created: $MANIFEST_FILE"
    echo "   Contains $file_count VCF files"
}

# ===== KIỂM TRA KẾT QUẢ =====
verify_results() {
    echo "🔍 Kiểm tra kết quả cuối cùng..."
    
    local total_original=$(find "$VCF_DIR" -maxdepth 1 -name "block_*.vcf" | wc -l)
    local total_gz=$(find "$OUTPUT_DIR" -name "block_*.vcf.gz" | wc -l)
    local total_tbi=$(find "$OUTPUT_DIR" -name "block_*.vcf.gz.tbi" | wc -l)
    
    echo "📊 Thống kê:"
    echo "   Original VCF files: $total_original"
    echo "   Compressed VCF files: $total_gz"
    echo "   Index files: $total_tbi"
    
    if [[ $total_gz -eq $total_tbi ]]; then
        echo "✅ Tất cả files đã được xử lý thành công!"
    else
        echo "⚠️  Có sự khác biệt trong số lượng files"
    fi
    
    # Hiển thị 5 dòng đầu của manifest
    echo ""
    echo "📋 Manifest preview (first 5 lines):"
    head -5 "$MANIFEST_FILE"
    
    # Hiển thị header của 1 file để kiểm tra
    echo ""
    echo "🔍 Kiểm tra header của 1 file VCF:"
    local sample_file=$(find "$OUTPUT_DIR" -name "block_*.vcf.gz" | head -1)
    if [[ -n "$sample_file" ]]; then
        echo "File: $(basename "$sample_file")"
        bcftools view -h "$sample_file" | head -8
        echo "..."
        echo "Sample count: $(bcftools view -h "$sample_file" | grep "^#CHROM" | head -1 | awk '{print NF-9}')"
        echo "Variant count: $(bcftools view "$sample_file" | grep -vc "^#")"
    fi
}

# ===== TẠO SCRIPT MERGE =====
create_merge_script() {
    echo "🔄 Tạo script merge VCF..."
    
    local merge_script="$OUTPUT_DIR/merge_vcfs.sh"
    
    cat > "$merge_script" << EOF
#!/bin/bash

# Script: merge_vcfs.sh
# Description: Merge tất cả VCF files đã nén
# Usage: ./merge_vcfs.sh

VCF_DIR="$OUTPUT_DIR"
MANIFEST="\$VCF_DIR/vcf_manifest.txt"
OUTPUT_VCF="\$VCF_DIR/merged_variants.vcf.gz"

echo "🔄 Merging VCF files..."

# Kiểm tra bcftools
if ! command -v bcftools &> /dev/null; then
    echo "❌ bcftools không được cài đặt. Cài đặt: conda install -c bioconda bcftools"
    exit 1
fi

# Kiểm tra manifest
if [[ ! -f "\$MANIFEST" ]]; then
    echo "❌ Manifest file không tồn tại: \$MANIFEST"
    exit 1
fi

# Kiểm tra số files
local file_count=\$(wc -l < "\$MANIFEST")
echo "📁 Số VCF files để merge: \$file_count"

# Sử dụng bcftools concat để merge
bcftools concat \\
    --file-list "\$MANIFEST" \\
    --output-type z \\
    --output "\$OUTPUT_VCF"

if [[ \$? -eq 0 ]]; then
    echo "✅ Đã merge VCF files: \$OUTPUT_VCF"
    
    # Index file merged
    tabix -p vcf "\$OUTPUT_VCF"
    echo "✅ Đã index merged VCF"
    
    # Thống kê
    echo "📊 Thống kê merged VCF:"
    echo "   Variant count: \$(bcftools view "\$OUTPUT_VCF" | grep -v "^#" | wc -l)"
    echo "   Number of contigs: \$(bcftools view -h "\$OUTPUT_VCF" | grep "^##contig" | wc -l)"
    echo "   Sample count: \$(bcftools view -h "\$OUTPUT_VCF" | grep "^#CHROM" | head -1 | awk '{print NF-9}')"
else
    echo "❌ Lỗi khi merge VCF files"
    exit 1
fi
EOF

    chmod +x "$merge_script"
    echo "✅ Merge script created: $merge_script"
}

# ===== DỌN DẸP FILE GỐC TRỐNG =====
cleanup_empty_files() {
    echo "🧹 Dọn dẹp files trống trong thư mục gốc..."
    
    local empty_count=0
    local total_remaining=0
    
    # Đếm số file trống còn lại và xóa chúng
    for vcf_file in "$VCF_DIR"/block_*.vcf; do
        if [[ -f "$vcf_file" ]]; then
            local variant_count=$(grep -vc "^#" "$vcf_file" 2>/dev/null || echo 0)
            if [[ $variant_count -eq 0 ]]; then
                echo "🗑️  Xoá file trống: $(basename "$vcf_file")"
                rm "$vcf_file"
                ((empty_count++))
            else
                ((total_remaining++))
            fi
        fi
    done
    
    echo "✅ Đã xoá $empty_count files trống"
    echo "📁 Còn lại $total_remaining files trong thư mục gốc"
}

# ===== MAIN =====
main() {
    echo "🚀 BẮT ĐẦU XỬ LÝ VCF FILES - STRICT VERSION"
    echo "=========================================="
    echo "🔧 CHẾ ĐỘ:"
    echo "   ✅ CHỈ giữ files có variants VÀ samples"
    echo "   🗑️  XOÁ files trống (không có variants)"
    echo "   ⚠️  Bỏ qua files có variants nhưng không có samples"
    echo "=========================================="
    
    START_TIME=$(date +%s)
    
    # Thực hiện các bước
    check_tools
    create_directories
    filter_vcf_files
    fix_vcf_headers
    normalize_sample_names
    validate_vcf_files
    compress_vcf_files
    verify_vcf_headers
    create_manifest
    verify_results
    create_merge_script
    cleanup_empty_files
    
    END_TIME=$(date +%s)
    TOTAL_TIME=$((END_TIME - START_TIME))
    
    echo ""
    echo "🎉 HOÀN THÀNH!"
    echo "=========================================="
    echo "📁 Output directory: $OUTPUT_DIR"
    echo "📋 Manifest file: $MANIFEST_FILE"
    echo "⏱️  Thời gian: $((TOTAL_TIME / 60)) phút $((TOTAL_TIME % 60)) giây"
    echo ""
    echo "🔜 NEXT STEPS:"
    echo "1. Chạy merge script:"
    echo "   cd $OUTPUT_DIR && ./merge_vcfs.sh"
    echo "2. Tạo graph với vg construct:"
    echo "   vg construct -r reference.fasta -v merged_variants.vcf.gz -m 32 > graph.vg"
    echo ""
    echo "⚠️  LƯU Ý: Tất cả files trống đã bị xoá khỏi thư mục gốc"
}

# Chạy main function
main
