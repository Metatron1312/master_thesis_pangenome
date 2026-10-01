#!/bin/bash

# ------------------------------------------
# CẤU HÌNH CHUNG
# ------------------------------------------
INPUT_DIR="/Users/shinra/gene_detection_v2/validate_database/SAMN49741230/raw_reads"                       # Thư mục đầu vào (truyền qua tham số khi chạy script)
FASTQC_RESULT_DIR="/Users/shinra/gene_detection_v2/validate_database/SAMN49741230/qc_checked"   # Thư mục output
THREADS=8                           # Số luồng xử lý
# Ví dụ thực hiện: #./fastqc_pipeline.sh /đường/dẫn/đến/thư/mục/quốc_gia kết quả sẽ được tạo tại nơi chạy script
# ------------------------------------------
# HÀM PHÁT HIỆN VÀ ĐỔI TÊN LONG READS
# ------------------------------------------
mkdir -p $FASTQC_RESULT_DIR
function detect_long_reads() {
    echo "=== PHÁT HIỆN LONG READS ==="
    
    if [ ! -d "$INPUT_DIR" ]; then
        echo "Lỗi: Thư mục không tồn tại - $INPUT_DIR"
        exit 1
    fi

    find "$INPUT_DIR" -type f -name "*_1.fastq.gz" | while read r1_file; do
        r2_file=$(echo "$r1_file" | sed 's/_1\.fastq\.gz$/_2.fastq.gz/')
        
        if [ ! -f "$r2_file" ]; then
            new_name=$(echo "$r1_file" | sed 's/_1\.fastq\.gz$/_longread.fastq.gz/')
            echo "Đổi tên: $(basename $r1_file) -> $(basename $new_name)"
            mv "$r1_file" "$new_name"
        fi
    done
    
    echo "Hoàn thành phát hiện long reads!"
    echo "--------------------------------"
}

# ------------------------------------------
# HÀM CHẠY FASTQC
# ------------------------------------------
function run_fastqc() {
    echo "=== CHẠY FASTQC ==="
    
    mkdir -p "${FASTQC_RESULT_DIR}"/{short_reads,long_reads}
    
    # Xử lý short reads
    echo "Đang xử lý short reads (R1 và R2)..."
    find "$INPUT_DIR" -name "*_1.fastq.gz" | parallel -j $THREADS "fastqc {} -o ${FASTQC_RESULT_DIR}/short_reads"
    find "$INPUT_DIR" -name "*_2.fastq.gz" | parallel -j $THREADS "fastqc {} -o ${FASTQC_RESULT_DIR}/short_reads"
    
    # Xử lý long reads
    echo "Đang xử lý long reads..."
    find "$INPUT_DIR" -name "*_longread.fastq.gz" | parallel -j $((THREADS/2)) "fastqc {} -o ${FASTQC_RESULT_DIR}/long_reads"
    
    echo "Hoàn thành chạy FASTQC!"
    echo "-----------------------"
}

# ------------------------------------------
# HÀM KIỂM TRA VÀ CHẠY LẠI FASTQC
# ------------------------------------------
function check_fastqc() {
    echo "=== KIỂM TRA KẾT QUẢ FASTQC ==="
    
    # Kiểm tra short reads
    find "$INPUT_DIR" -name "*_[12].fastq.gz" | while read file; do
        sample_name=$(basename "$file" .fastq.gz)
        output_file="${FASTQC_RESULT_DIR}/short_reads/${sample_name}_fastqc.html"
        
        if [ ! -f "$output_file" ]; then
            echo "[THIẾU] Chạy lại: $sample_name"
            fastqc "$file" -o "${FASTQC_RESULT_DIR}/short_reads" -t $((THREADS/2))
        fi
    done
    
    # Kiểm tra long reads
    find "$INPUT_DIR" -name "*_longread.fastq.gz" | while read file; do
        sample_name=$(basename "$file" .fastq.gz)
        output_file="${FASTQC_RESULT_DIR}/long_reads/${sample_name}_fastqc.html"
        
        if [ ! -f "$output_file" ]; then
            echo "[THIẾU] Chạy lại: $sample_name"
            fastqc "$file" -o "${FASTQC_RESULT_DIR}/long_reads" -t $((THREADS/2))
        fi
    done
    
    echo "Hoàn thành kiểm tra FASTQC!"
    echo "--------------------------"
}

# ------------------------------------------
# MENU CHÍNH
# ------------------------------------------
function main() {
    if [ -z "$INPUT_DIR" ]; then
        echo "Cách dùng: $0 /đường/dẫn/đến/thư/mục/dữ/liệu"
        exit 1
    fi

    echo "=== CHƯƠNG TRÌNH XỬ LÝ FASTQC ==="
    echo "Thư mục đầu vào: $INPUT_DIR"
    echo "1. Phát hiện và đổi tên long reads"
    echo "2. Chạy FASTQC ban đầu"
    echo "3. Kiểm tra và chạy lại FASTQC"
    echo "4. Chạy toàn bộ quy trình"
    echo "0. Thoát"
    
    read -p "Chọn chức năng (0-4): " choice
    
    case $choice in
        1) detect_long_reads ;;
        2) run_fastqc ;;
        3) check_fastqc ;;
        4) 
            detect_long_reads
            run_fastqc
            check_fastqc
            ;;
        0) exit 0 ;;
        *) echo "Lựa chọn không hợp lệ!" ;;
    esac
}

# Chạy chương trình
main
