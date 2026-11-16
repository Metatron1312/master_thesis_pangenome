docker run --platform linux/amd64 \
  -v /Users/shinra/Bioinformatics/WGS_AMR/data_public/combined_contigs:/data \
  -v /Users/shinra/Mash_output:/output \
  quay.io/biocontainers/mash:2.3--ha61e061_0 \
  bash -c "
  # 1. Tạo sketch từ danh sách file
  find /data -name '*.filtered.complete.fasta' > /output/file_list.txt
  mash sketch -p 4 -s 10000 -o /output/all_genomes -l /output/file_list.txt
  
  # 2. Tính khoảng cách với xử lý p-value
  (echo -e 'Sample1\tSample2\tDistance\tP-value\tShared_hashes'
  mash dist -p 4 /output/all_genomes.msh /output/all_genomes.msh | \
  awk 'BEGIN{OFS=\"\t\"} {
    # Xử lý tên file để lấy ID mẫu
    split(\$1,a,\"/\"); id1=a[length(a)]; sub(/\\.filtered\\.complete\\.fasta$/,\"\",id1);
    split(\$2,b,\"/\"); id2=b[length(b)]; sub(/\\.filtered\\.complete\\.fasta$/,\"\",id2);
    
    # Xử lý p-value: thay giá trị 0 bằng giá trị nhỏ nhất có thể
    pval = \$4;
    if (pval == 0) pval = 2.22507e-308; # Giá trị double nhỏ nhất
    
    if (id1 != id2) print id1, id2, \$3, pval, \$5
  }') > /output/divergence_matrix.tsv
  
  # 3. Thống kê chi tiết
  awk 'BEGIN{FS=\"\t\"; total=0; sum=0; max=0; min=1; sig=0; lt5=0; gte5_lt10=0; gte10_lte20=0; gt20=0} 
  NR>1 {
    dist=\$3; 
    total++; 
    sum+=dist; 
    if (dist>max) max=dist; 
    if (dist<min) min=dist;
    if (\$4 < 0.05) sig++;
    if (dist < 0.05) lt5++;
    else if (dist < 0.1) gte5_lt10++;
    else if (dist <= 0.2) gte10_lte20++;
    else gt20++;
  } END {
    if(total==0) {print \"Không có dữ liệu thống kê\" > \"/output/divergence_stats.txt\"; exit};
    mean=sum/total;
    print \"THỐNG KÊ PHÂN TÍCH DIVERGENCE\" > \"/output/divergence_stats.txt\";
    print \"================================\" >> \"/output/divergence_stats.txt\";
    print \"Tổng số mẫu:\", NR >> \"/output/divergence_stats.txt\";
    print \"Tổng cặp so sánh:\", total >> \"/output/divergence_stats.txt\";
    print \"Giá trị trung bình:\", sprintf(\"%.6f\", mean) >> \"/output/divergence_stats.txt\";
    print \"Khoảng cách nhỏ nhất:\", sprintf(\"%.6f\", min) >> \"/output/divergence_stats.txt\";
    print \"Khoảng cách lớn nhất:\", sprintf(\"%.6f\", max) >> \"/output/divergence_stats.txt\";
    
    print \"\nPHÂN BỐ KHOẢNG CÁCH\" >> \"/output/divergence_stats.txt\";
    printf \"- Dưới 5%%:\t\t%6d (%.2f%%)\\n\", lt5, (lt5/total)*100 >> \"/output/divergence_stats.txt\";
    printf \"- Từ 5%% đến dưới 10%%:\t%6d (%.2f%%)\\n\", gte5_lt10, (gte5_lt10/total)*100 >> \"/output/divergence_stats.txt\";
    printf \"- Từ 10%% đến 20%%:\t%6d (%.2f%%)\\n\", gte10_lte20, (gte10_lte20/total)*100 >> \"/output/divergence_stats.txt\";
    printf \"- Trên 20%%:\t\t%6d (%.2f%%)\\n\", gt20, (gt20/total)*100 >> \"/output/divergence_stats.txt\";
    
    print \"\nÝ NGHĨA THỐNG KÊ\" >> \"/output/divergence_stats.txt\";
    print \"Cặp có p-value < 0.05:\", sig, sprintf(\"(%.2f%%)\", (sig/total)*100) >> \"/output/divergence_stats.txt\";
    print \"\nGhi chú: Các giá trị p-value bằng 0 đã được thay thế bằng giá trị nhỏ nhất có thể (2.22507e-308)\" >> \"/output/divergence_stats.txt\";
  }' /output/divergence_matrix.tsv"