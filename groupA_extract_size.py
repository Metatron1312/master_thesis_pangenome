import os
import re

# ============================================================
# CẤU HÌNH
# ============================================================

# File .log của hai cây (chứa thông tin alignment đầy đủ)
CORE_GENE_LOG = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/iqtree/core_gene/core_gene_tree.log"
CORE_GENOME_LOG = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/iqtree/core_genome/core_genome_tree.log"

# File output
OUTPUT_DIR = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/comparison"
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_TSV = os.path.join(OUTPUT_DIR, "groupA_size.tsv")

# ============================================================
# HÀM TRÍCH XUẤT
# ============================================================

def extract_size_info(log_file, label):
    """
    Trích xuất số mẫu, số cột, distinct patterns từ file .log
    Format: "Alignment has 107 sequences with 1311366 columns, 30489 distinct patterns"
    """
    if not os.path.exists(log_file):
        print(f"LỖI: Không tìm thấy {log_file}")
        return None
    
    with open(log_file) as f:
        content = f.read()
    
    match = re.search(
        r'Alignment has (\d+) sequences with (\d+) columns, (\d+) distinct patterns',
        content
    )
    
    if not match:
        print(f"LỖI: Không tìm thấy thông tin alignment trong {log_file}")
        return None
    
    n_sequences = int(match.group(1))
    n_columns = int(match.group(2))
    n_patterns = int(match.group(3))
    
    return {
        'label': label,
        'n_sequences': n_sequences,
        'n_columns': n_columns,
        'n_patterns': n_patterns,
    }

# ============================================================
# TRÍCH XUẤT
# ============================================================

print("=" * 60)
print("KÍCH THƯỚC DỮ LIỆU")
print("=" * 60)

core_gene_info = extract_size_info(CORE_GENE_LOG, "core_gene")
core_genome_info = extract_size_info(CORE_GENOME_LOG, "core_genome")

if core_gene_info is None or core_genome_info is None:
    print("LỖI: Không thể trích xuất thông tin!")
    exit(1)

# ============================================================
# IN KẾT QUẢ
# ============================================================

print(f"\n{'Thông số':<25} {'Core Gene':>15} {'Core Genome':>15} {'Khác biệt':>15}")
print("-" * 75)

for key, label in [
    ('n_sequences', 'Số mẫu'),
    ('n_columns', 'Số cột'),
    ('n_patterns', 'Distinct patterns'),
]:
    cg = core_gene_info[key]
    cge = core_genome_info[key]
    diff = cge - cg
    diff_pct = (diff / cg * 100) if cg != 0 else 0
    print(f"{label:<25} {cg:>15,} {cge:>15,} {diff:>+10,} ({diff_pct:+.1f}%)")

# ============================================================
# LƯU FILE TSV
# ============================================================

with open(OUTPUT_TSV, "w") as f:
    f.write("metric\tcore_gene\tcore_genome\tdifference\tpercent_change\n")
    for key, label in [
        ('n_sequences', 'Số mẫu'),
        ('n_columns', 'Số cột'),
        ('n_patterns', 'Distinct patterns'),
    ]:
        cg = core_gene_info[key]
        cge = core_genome_info[key]
        diff = cge - cg
        diff_pct = (diff / cg * 100) if cg != 0 else 0
        f.write(f"{label}\t{cg}\t{cge}\t{diff}\t{diff_pct:.2f}\n")

print(f"\nĐã lưu kết quả vào: {OUTPUT_TSV}")