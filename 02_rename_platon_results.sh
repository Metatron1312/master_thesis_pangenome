#!/bin/bash

# Script đổi tên PERMANENT - XÓA FILE GỐC chứa '.temp'
# (Chỉ giữ lại file mới không có '.temp' trong tên)

# Cấu hình
TARGET_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/Platon_results"  # Thay bằng đường dẫn thực tế

# ---- KHÔNG CHỈNH SỬA PHẦN DƯỚI ĐÂY ----
echo "🔍 Đang quét thư mục: $TARGET_DIR"
total_renamed=0
total_skipped=0

# Xử lý đệ quy
find "$TARGET_DIR" -type f -name "*temp*" | while read -r old_path; do
    # Tạo tên mới (loại bỏ ALL '.temp')
    new_path="${old_path//.temp/}"
    
    # Kiểm tra trùng lặp
    if [ "$old_path" != "$new_path" ]; then
        if [ ! -e "$new_path" ]; then
            # Đổi tên PERMANENT
            if mv -v "$old_path" "$new_path"; then
                ((total_renamed++))
                echo "✅ Đã đổi: $old_path → $new_path"
            else
                echo "❌ Lỗi khi đổi $old_path" >&2
            fi
        else
            echo "⚠️  Bỏ qua: $new_path đã tồn tại" >&2
            ((total_skipped++))
        fi
    fi
done

echo "📊 Kết quả:"
echo "- Đã đổi tên: $total_renamed file"
echo "- Bỏ qua: $total_skipped file (đích đã tồn tại)"
echo "🎉 Hoàn thành!"
