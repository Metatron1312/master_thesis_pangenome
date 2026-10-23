import os
import re

# ============================================================
# CẤU HÌNH
# ============================================================

CORE_GENE_LOG = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/iqtree/core_gene/core_gene_tree.log"
CORE_GENOME_LOG = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/iqtree/core_genome/core_genome_tree.log"

OUTPUT_DIR = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/comparison"
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_TSV = os.path.join(OUTPUT_DIR, "groupB_quality.tsv")

# ============================================================
# HÀM TRÍCH XUẤT
# ============================================================

def extract_quality_info(log_file, label):
    """
    Trích xuất thông tin chất lượng alignment từ file .log
    Format:
      Alignment has 107 sequences with 1311366 columns, 30489 distinct patterns
      60134 parsimony-informative, 14638 singleton sites, 1236594 constant sites
    """
    if not os.path.exists(log_file):
        print(f"LỖI: Không tìm thấy {log_file}")
        return None
    
    with open(log_file) as f:
        content = f.read()
    
    # Số cột và distinct patterns
    match1 = re.search(
        r'Alignment has (\d+) sequences with (\d+) columns, (\d+) distinct patterns',
        content
    )
    
    # Parsimony-informative, singleton, constant
    match2 = re.search(
        r'(\d+) parsimony-informative, (\d+) singleton sites, (\d+) constant sites',
        content
    )
    
    if not match1 or not match2:
        print(f"LỖI: Không tìm thấy thông tin trong {log_file}")
        return None
    
    n_columns = int(match1.group(2))
    n_patterns = int(match1.group(3))
    n_parsimony = int(match2.group(1))
    n_singleton = int(match2.group(2))
    n_constant = int(match2.group(3))
    
    return {
        'label': label,
        'n_columns': n_columns,
        'n_patterns': n_patterns,
        'n_parsimony': n_parsimony,
        'n_singleton': n_singleton,
        'n_constant': n_constant,
    }

# ============================================================
# TRÍCH XUẤT
# ============================================================

print("=" * 70)
print("NHÓM B: CHẤT LƯỢNG ALIGNMENT")
print("=" * 70)

core_gene_info = extract_quality_info(CORE_GENE_LOG, "core_gene")
core_genome_info = extract_quality_info(CORE_GENOME_LOG, "core_genome")

if core_gene_info is None or core_genome_info is None:
    print("LỖI: Không thể trích xuất thông tin!")
    exit(1)

# ============================================================
# IN KẾT QUẢ
# ============================================================

print(f"\n{'Thông số':<30} {'Core Gene':>15} {'Core Genome':>15} {'Khác biệt':>15}")
print("-" * 80)

for key, label in [
    ('n_parsimony', 'Parsimony-informative'),
    ('n_constant', 'Constant sites'),
    ('n_singleton', 'Singleton sites'),
    ('n_columns', 'Tổng số cột'),
]:
    cg = core_gene_info[key]
    cge = core_genome_info[key]
    diff = cge - cg
    diff_pct = (diff / cg * 100) if cg != 0 else 0
    print(f"{label:<30} {cg:>15,} {cge:>15,} {diff:>+10,} ({diff_pct:+.1f}%)")

# ============================================================
# TÍNH TỶ LỆ
# ============================================================

print(f"\n{'Tỷ lệ':<30} {'Core Gene':>15} {'Core Genome':>15}")
print("-" * 65)

for key, label in [
    ('n_parsimony', 'Parsimony-informative / cột'),
    ('n_constant', 'Constant / cột'),
    ('n_singleton', 'Singleton / cột'),
]:
    cg = core_gene_info[key]
    cge = core_genome_info[key]
    cg_pct = cg / core_gene_info['n_columns'] * 100
    cge_pct = cge / core_genome_info['n_columns'] * 100
    print(f"{label:<30} {cg_pct:>14.2f}% {cge_pct:>14.2f}%")

# ============================================================
# LƯU FILE TSV
# ============================================================

with open(OUTPUT_TSV, "w") as f:
    f.write("metric\tcore_gene\tcore_genome\tdifference\tpercent_change\n")
    for key, label in [
        ('n_parsimony', 'Parsimony-informative'),
        ('n_constant', 'Constant sites'),
        ('n_singleton', 'Singleton sites'),
        ('n_columns', 'Tổng số cột'),
    ]:
        cg = core_gene_info[key]
        cge = core_genome_info[key]
        diff = cge - cg
        diff_pct = (diff / cg * 100) if cg != 0 else 0
        f.write(f"{label}\t{cg}\t{cge}\t{diff}\t{diff_pct:.2f}\n")

print(f"\nĐã lưu kết quả vào: {OUTPUT_TSV}")
