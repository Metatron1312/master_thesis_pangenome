import os
import glob
import pandas as pd

def clean_scientific_name(name):
    """Xử lý khoảng trắng trong scientific_name"""
    return ' '.join(str(name).strip().split())

# Thay đổi đường dẫn này thành thư mục chứa file report của bạn
report_dir = "/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur/kraken2_results"
all_reports = glob.glob(f"{report_dir}/*.report")

combined_df = pd.DataFrame()

for file in all_reports:
    try:
        filename = os.path.basename(file)
        
        # Trích xuất biosample_id (AB_Precipitated_Purple_20)
        biosample_id = '_'.join(filename.split('_')[:4])
        
        # Trích xuất file_type (plasmid hoặc chromosome)
        file_type = filename.split('_')[-1].split('.')[0]
        
        df = pd.read_csv(file, sep='\t', header=None,
                        names=['percentage', 'clade_fragments', 'direct_fragments',
                               'rank_code', 'taxonomy_id', 'scientific_name'])
        
        # Xử lý khoảng trắng
        df['scientific_name'] = df['scientific_name'].apply(clean_scientific_name)
        
        df['biosample_id'] = biosample_id
        df['file_type'] = file_type
        combined_df = pd.concat([combined_df, df], ignore_index=True)
        
        print(f"✅ Processed: {filename}")
        
    except Exception as e:
        print(f"⚠️ Error processing {filename}: {str(e)}")
        continue

if not combined_df.empty:
    # Sắp xếp lại các cột
    combined_df = combined_df[['biosample_id', 'file_type', 'percentage', 'clade_fragments',
                             'direct_fragments', 'rank_code', 'taxonomy_id', 'scientific_name']]
    
    # Lưu file CSV
    output_path = os.path.join(report_dir, "combined_reports_v2.csv")
    combined_df.to_csv(output_path, index=False)
    print(f"\n🎉 Success! Combined report saved to:\n{output_path}")
    
    # Hiển thị preview
    print("\nPreview of the first 5 rows:")
    print(combined_df.head())
else:
    print("\n❌ No valid data was processed. Please check your input files.")