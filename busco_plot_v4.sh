#!/bin/bash

# ================= CẤU HÌNH CHUẨN =================
OUTPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/busco_results"
PLOT_DIR="$OUTPUT_DIR/busco_plots_final"
SUMMARY_DIR="$OUTPUT_DIR/busco_summaries_$(date +%s)"
BUSCO_IMAGE="ezlabgva/busco:v6.0.0_cv1"
# ==================================================

# 1. Tạo thư mục làm việc
mkdir -p "$PLOT_DIR" "$SUMMARY_DIR"
echo -e "📂 Thư mục kết quả sẽ được lưu tại: $PLOT_DIR"

# 2. Tập hợp tất cả file json từ các thư mục con
echo -e "\n🔍 Đang thu thập các file kết quả BUSCO..."
find "$OUTPUT_DIR" -type f -name "short_summary.*.json" -exec cp {} "$SUMMARY_DIR" \;
TOTAL_FILES=$(ls -1 "$SUMMARY_DIR"/*.json 2>/dev/null | wc -l)

if [ "$TOTAL_FILES" -eq 0 ]; then
    echo -e "\n❌ Không tìm thấy file JSON nào trong các thư mục con"
    echo "ℹ️ Hãy kiểm tra cấu trúc thư mục:"
    find "$OUTPUT_DIR" -maxdepth 2 -type d -name "busco_*" | head -5
    exit 1
fi

# 3. Tạo plot với các tùy chọn tối ưu
echo -e "\n📊 Đang tạo plot cho $TOTAL_FILES mẫu..."
docker run --rm --platform linux/amd64 \
  -v "$SUMMARY_DIR":/summaries \
  -v "$PLOT_DIR":/plots \
  -v "$OUTPUT_DIR":/original_results \
  "$BUSCO_IMAGE" \
  bash -c "
  # Chạy lệnh plot chính thức
  busco --plot /summaries --out_path /plots || true
  
  # Nếu không tìm thấy plot ở vị trí mong đợi, tìm kiếm trong toàn bộ container
  if [ ! -f /plots/busco_figure*.png ]; then
    echo 'ℹ️ Đang tìm kiếm plot trong container...'
    find / -name 'busco_figure*.png' -exec cp {} /plots/ \; 2>/dev/null || true
  fi
  "

# 4. Tìm kiếm kết quả ở cả 3 vị trí có thể
RESULTS_FOUND=0
POSSIBLE_PATHS=(
  "$PLOT_DIR"
  "$SUMMARY_DIR" 
  "$OUTPUT_DIR"
)

for path in "${POSSIBLE_PATHS[@]}"; do
  if ls "$path"/busco_figure*.png 1> /dev/null 2>&1; then
    ((RESULTS_FOUND++))
    echo -e "\n✅ Tìm thấy plot tại: $path"
    cp "$path"/busco_figure*.png "$PLOT_DIR/" 2>/dev/null
    ls -lh "$PLOT_DIR"/busco_figure*.png | head -5
  fi
done

# 5. Xử lý khi không tìm thấy
if [ "$RESULTS_FOUND" -eq 0 ]; then
  echo -e "\n❌ Không tìm thấy file plot nào. Nguyên nhân có thể:"
  echo "1. Vấn đề với thư viện matplotlib trong container"
  echo "2. File JSON không đúng định dạng chuẩn"
  echo "3. Lỗi phiên bản BUSCO"
  
  # Tạo plot thủ công từ file JSON
  echo -e "\n🛠 Đang thử tạo plot thủ công từ file JSON..."
  docker run --rm --platform linux/amd64 \
    -v "$SUMMARY_DIR":/data \
    -v "$PLOT_DIR":/plots \
    python:3.9-slim \
    bash -c "
    pip install matplotlib pandas --quiet > /dev/null 2>&1 &&
    python -c \"
import glob, json, os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Thu thập dữ liệu
data = []
for f in glob.glob('/data/short_summary.*.json'):
    try:
        with open(f) as jf:
            d = json.load(jf)
        sample = os.path.basename(f).split('.')[-2]
        data.append({
            'Sample': sample[:30],  # Giới hạn độ dài tên
            'Complete': d['results']['Complete'],
            'Fragmented': d['results']['Fragmented'],
            'Missing': d['results']['Missing']
        })
    except Exception as e:
        print(f'Error processing {f}: {str(e)}')

# Tạo DataFrame và sắp xếp
df = pd.DataFrame(data).set_index('Sample').sort_index()

# Tạo PDF chuyên nghiệp
pdf_path = '/plots/busco_results.pdf'
with PdfPages(pdf_path) as pdf:
    # Tạo trang tổng quan
    plt.figure(figsize=(12, 8))
    df.mean().plot(kind='pie', autopct='%1.1f%%', 
                  colors=['#1f77b4', '#ff7f0e', '#d62728'],
                  title='Average BUSCO Results (All Samples)')
    pdf.savefig(bbox_inches='tight')
    plt.close()
    
    # Tạo từng trang chi tiết
    samples_per_page = 10
    total_pages = (len(df) + samples_per_page - 1) // samples_per_page
    
    for page in range(total_pages):
        start_idx = page * samples_per_page
        end_idx = (page + 1) * samples_per_page
        df_page = df.iloc[start_idx:end_idx]
        
        plt.figure(figsize=(12, 8))
        df_page[['Complete', 'Fragmented', 'Missing']].plot(
            kind='bar', stacked=True, 
            color=['#1f77b4', '#ff7f0e', '#d62728']
        )
        
        plt.title(f'BUSCO Results - Samples {start_idx+1} to {min(end_idx, len(df))}')
        plt.ylabel('Percentage (%)')
        plt.ylim(0, 100)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        pdf.savefig()
        plt.close()
    
    # Thêm metadata
    metadata = pdf.infodict()
    metadata['Title'] = 'BUSCO Results Report'
    metadata['Author'] = 'Automated BUSCO Analysis'
    metadata['Subject'] = 'Genome Assessment Results'
    metadata['Keywords'] = 'BUSCO WGS Genomics'

print(f'✅ Đã tạo file PDF chuyên nghiệp: {pdf_path}')
\""
fi

# 6. Tạo file PDF tổng hợp từ các plot PNG (nếu có)
if [ "$RESULTS_FOUND" -gt 0 ]; then
  echo -e "\n📄 Đang tạo file PDF tổng hợp từ các plot..."
  docker run --rm --platform linux/amd64 \
    -v "$PLOT_DIR":/plots \
    python:3.9-slim \
    bash -c "
    pip install matplotlib pillow reportlab --quiet > /dev/null 2>&1 &&
    python -c \"
import os
import glob
from PIL import Image
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter, landscape
from reportlab.lib.utils import ImageReader

# Thu thập tất cả các file plot
plot_files = sorted(glob.glob('/plots/busco_figure*.png'))
if not plot_files:
    exit(0)

# Tạo file PDF
pdf_path = '/plots/busco_plots_combined.pdf'
c = canvas.Canvas(pdf_path, pagesize=landscape(letter))

# Thiết lập thông số
margin = 40
page_width, page_height = landscape(letter)
usable_width = page_width - 2 * margin
usable_height = page_height - 2 * margin

# Thêm trang tiêu đề
c.setFont('Helvetica-Bold', 24)
c.drawCentredString(page_width/2, page_height/2, 'BUSCO Results Visualization')
c.setFont('Helvetica', 16)
c.drawCentredString(page_width/2, page_height/2 - 40, 
                   f'Total Samples: {len(plot_files)}')
c.showPage()

# Thêm từng ảnh vào PDF
for i, img_path in enumerate(plot_files):
    try:
        img = Image.open(img_path)
        img_width, img_height = img.size
        
        # Tính toán tỷ lệ để vừa trang
        scale = min(usable_width/img_width, usable_height/img_height) * 0.9
        new_width = img_width * scale
        new_height = img_height * scale
        
        # Đặt ảnh vào giữa trang
        x = (page_width - new_width) / 2
        y = (page_height - new_height) / 2
        
        # Thêm ảnh
        c.drawImage(ImageReader(img), x, y, width=new_width, height=new_height)
        c.setFont('Helvetica', 10)
        c.drawString(margin, margin, f'Page {i+1} - {os.path.basename(img_path)}')
        c.showPage()
    except Exception as e:
        print(f'Error processing {img_path}: {str(e)}')

# Thêm trang tổng kết
c.setFont('Helvetica-Bold', 18)
c.drawCentredString(page_width/2, page_height/2, '--- End of Report ---')
c.setFont('Helvetica', 12)
c.drawCentredString(page_width/2, page_height/2 - 30, 
                   'Generated by Automated BUSCO Analysis Pipeline')
c.save()

print(f'✅ Đã tạo file PDF tổng hợp: {pdf_path}')
\""
fi

# 7. Xóa thư mục summary
echo -e "\n🧹 Đang dọn dẹp thư mục summary..."
if [ -d "$SUMMARY_DIR" ]; then
    rm -rf "$SUMMARY_DIR"
    echo "✅ Đã xóa thư mục summary: $SUMMARY_DIR"
else
    echo "ℹ️ Không tìm thấy thư mục summary để xóa"
fi

# 8. Hiển thị kết quả cuối cùng
echo -e "\n📝 KẾT QUẢ CUỐI CÙNG:"
find "$PLOT_DIR" -type f \( -name "*.png" -o -name "*.pdf" \) | while read -r file; do
    echo -e "   - ${file}"
done

echo -e "\n🎉 HOÀN TẤT! Tất cả kết quả đã được lưu tại: $PLOT_DIR"