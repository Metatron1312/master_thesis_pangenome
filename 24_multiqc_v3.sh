#!/bin/bash

# Cấu hình
DATA_ROOT="/Users/shinra/gene_detection_v2/validate_database/SAMN49741230/qc_checked/short_reads"
WORK_DIR="multiqc_temp"
REPORT_DIR="/Users/shinra/gene_detection_v2/validate_database/SAMN49741230/qc_checked/short_reads/multiqc_report"
mkdir -p $WORK_DIR/short_reads $WORK_DIR/long_reads $REPORT_DIR

# Hàm trích xuất metadata chính xác
function get_metadata {
    local path=$1
    rel_path=${path#$DATA_ROOT/}
    IFS='/' read -ra parts <<< "$rel_path"
    
    country="${parts[0]}"
    project="No_Project"
    [[ ${#parts[@]} -ge 2 ]] && project="${parts[1]}"
    
    echo -e "$country\t$project"
}

# Xử lý từng loại read
for read_type in short long; do
    # Tạo metadata
    echo -e "Sample\tCountry\tProject" > $WORK_DIR/${read_type}_reads/metadata.tsv
    
    # Quét và xử lý file
    find "$DATA_ROOT" -name "*_fastqc.zip" -path "*${read_type}_reads*" | while read file; do
        # Lấy metadata chính xác
        read -r country project <<< $(get_metadata "$file")
        sample=$(basename "$file" | sed 's/_fastqc\.zip//')
        
        # Ghi metadata
        echo -e "${sample}\t${country}\t${project}" >> $WORK_DIR/${read_type}_reads/metadata.tsv
        
        # Tạo symlink với tên chuẩn
        ln -sf "$(realpath "$file")" "$WORK_DIR/${read_type}_reads/${sample}_fastqc.zip"
    done

    # Tạo config với cấu hình đặc biệt
    read_type_upper=$(echo "$read_type" | tr '[:lower:]' '[:upper:]')
    cat <<EOF > $WORK_DIR/${read_type}_reads/multiqc_config.yaml
title: "BÁO CÁO ${read_type_upper}-READS"
use_filename_as_sample_name: false
strict_names: true
fastqc_config:
  max_samples: 0
  filter_samples: false
  override_sample_names: true  # Quan trọng: Đồng bộ tên mẫu
table_columns_visible:
  Country: true
  Project: true
report_header_info:
  - "Sample": "%(sample_name)s"
  - "Country": "%(Country)s"
  - "Project": "%(Project)s"
custom_logo: /path/to/logo.png
report_section_order:
  fastqc:
    order: -1000
custom_content:
  - name: "Hướng dẫn"
    content: |
      <h3>Dữ liệu ${read_type}-reads</h3>
      <p>Lọc theo Country/Project trong bảng General Statistics</p>
EOF
done

# Chạy MultiQC với cấu hình đồng bộ
for read_type in short long; do
    multiqc $WORK_DIR/${read_type}_reads/ \
        -o $REPORT_DIR/${read_type}_reads \
        --filename "${read_type}_reads_report" \
        --sample-names $WORK_DIR/${read_type}_reads/metadata.tsv \
        --config $WORK_DIR/${read_type}_reads/multiqc_config.yaml \
        --interactive --force \
        --cl-config "override_sample_names: true"  # Đảm bảo ghi đè tên mẫu
done

# Dọn dẹp
rm -rf $WORK_DIR

# Thông báo kết quả
echo "=== BÁO CÁO ĐÃ ĐƯỢC TẠO ==="
echo "- Short-reads: file://$(realpath $REPORT_DIR)/short_reads/short_reads_report.html"
echo "- Long-reads: file://$(realpath $REPORT_DIR)/long_reads/long_reads_report.html"
