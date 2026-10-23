import os
import re

# ============================================================
# CẤU HÌNH
# ============================================================

CORE_GENE_LOG = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/iqtree/core_gene/core_gene_tree.log"
CORE_GENOME_LOG = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/iqtree/core_genome/core_genome_tree.log"

CORE_GENE_TREEFILE = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/iqtree/core_gene/core_gene_tree.treefile"
CORE_GENOME_TREEFILE = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/iqtree/core_genome/core_genome_tree.treefile"

OUTPUT_DIR = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/comparison"
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_TSV = os.path.join(OUTPUT_DIR, "groupC_bootstrap.tsv")

# ============================================================
# HÀM TRÍCH XUẤT BOOTSTRAP TỪ TREEFILE
# ============================================================

def extract_bootstrap_from_treefile(treefile):
    """Trích xuất tất cả giá trị bootstrap từ file .treefile"""
    if not os.path.exists(treefile):
        return None
    
    with open(treefile) as f:
        content = f.read()
    
    bootstraps = re.findall(r'\)(\d+)[:\,]', content)
    bootstraps = [int(b) for b in bootstraps]
    
    return bootstraps

# ============================================================
# HÀM TRÍCH XUẤT CORRELATION VÀ TREE LENGTH TỪ LOG
# ============================================================

def extract_correlation(log_file):
    """Trích xuất bootstrap correlation coefficient"""
    if not os.path.exists(log_file):
        return None
    
    with open(log_file) as f:
        content = f.read()
    
    match = re.search(r'Bootstrap correlation coefficient of split occurrence frequencies:\s*([\d.]+)', content)
    return float(match.group(1)) if match else None

def extract_tree_length(log_file):
    """Trích xuất tổng chiều dài cây"""
    if not os.path.exists(log_file):
        return None
    
    with open(log_file) as f:
        content = f.read()
    
    match = re.search(r'Total tree length:\s*([\d.]+)', content)
    return float(match.group(1)) if match else None

# ============================================================
# TRÍCH XUẤT
# ============================================================

print("=" * 75)
print("NHÓM C: BOOTSTRAP + TREE LENGTH")
print("=" * 75)

cg_bootstraps = extract_bootstrap_from_treefile(CORE_GENE_TREEFILE)
cge_bootstraps = extract_bootstrap_from_treefile(CORE_GENOME_TREEFILE)

if cg_bootstraps is None or cge_bootstraps is None:
    print("LỖI: Không thể trích xuất bootstrap!")
    exit(1)

cg_corr = extract_correlation(CORE_GENE_LOG)
cge_corr = extract_correlation(CORE_GENOME_LOG)
cg_tree_len = extract_tree_length(CORE_GENE_LOG)
cge_tree_len = extract_tree_length(CORE_GENOME_LOG)

# ============================================================
# TÍNH TOÁN
# ============================================================

def analyze_bootstraps(bootstraps):
    n_total = len(bootstraps)
    n_ge95 = sum(1 for b in bootstraps if b >= 95)
    n_70_94 = sum(1 for b in bootstraps if 70 <= b < 95)
    n_lt70 = sum(1 for b in bootstraps if b < 70)
    mean_bs = sum(bootstraps) / n_total if n_total > 0 else 0
    return {
        'n_total': n_total,
        'n_ge95': n_ge95,
        'n_70_94': n_70_94,
        'n_lt70': n_lt70,
        'mean': mean_bs,
    }

cg_stats = analyze_bootstraps(cg_bootstraps)
cge_stats = analyze_bootstraps(cge_bootstraps)

# ============================================================
# IN KẾT QUẢ
# ============================================================

print(f"\n{'Thông số':<35} {'Core Gene':>15} {'Core Genome':>15} {'Khác biệt':>15}")
print("-" * 85)

# Correlation
print(f"{'Bootstrap correlation':<35} {cg_corr:>15.4f} {cge_corr:>15.4f} {cge_corr - cg_corr:>+10.4f}")

# Tree length
print(f"{'Tổng chiều dài cây':<35} {cg_tree_len:>15.4f} {cge_tree_len:>15.4f} {cge_tree_len - cg_tree_len:>+10.4f}")

# Tổng số nhánh
print(f"{'Tổng số nhánh':<35} {cg_stats['n_total']:>15,} {cge_stats['n_total']:>15,} {cge_stats['n_total'] - cg_stats['n_total']:>+10,}")

# Bootstrap >= 95
diff = cge_stats['n_ge95'] - cg_stats['n_ge95']
pct = diff / cg_stats['n_ge95'] * 100 if cg_stats['n_ge95'] > 0 else 0
print(f"{'Số nhánh bootstrap >= 95':<35} {cg_stats['n_ge95']:>15,} {cge_stats['n_ge95']:>15,} {diff:>+10,} ({pct:+.1f}%)")

# Bootstrap 70-94
diff = cge_stats['n_70_94'] - cg_stats['n_70_94']
pct = diff / cg_stats['n_70_94'] * 100 if cg_stats['n_70_94'] > 0 else 0
print(f"{'Số nhánh bootstrap 70-94':<35} {cg_stats['n_70_94']:>15,} {cge_stats['n_70_94']:>15,} {diff:>+10,} ({pct:+.1f}%)")

# Bootstrap < 70
diff = cge_stats['n_lt70'] - cg_stats['n_lt70']
pct = diff / cg_stats['n_lt70'] * 100 if cg_stats['n_lt70'] > 0 else 0
print(f"{'Số nhánh bootstrap < 70':<35} {cg_stats['n_lt70']:>15,} {cge_stats['n_lt70']:>15,} {diff:>+10,} ({pct:+.1f}%)")

# Bootstrap trung bình
diff = cge_stats['mean'] - cg_stats['mean']
pct = diff / cg_stats['mean'] * 100 if cg_stats['mean'] > 0 else 0
print(f"{'Bootstrap trung bình':<35} {cg_stats['mean']:>15.2f} {cge_stats['mean']:>15.2f} {diff:>+10.2f} ({pct:+.1f}%)")

# ============================================================
# LƯU FILE TSV
# ============================================================

with open(OUTPUT_TSV, "w") as f:
    f.write("metric\tcore_gene\tcore_genome\tdifference\tpercent_change\n")
    
    # Correlation
    f.write(f"Bootstrap correlation\t{cg_corr}\t{cge_corr}\t{cge_corr - cg_corr}\t0\n")
    
    # Tree length
    f.write(f"Tổng chiều dài cây\t{cg_tree_len}\t{cge_tree_len}\t{cge_tree_len - cg_tree_len}\t0\n")
    
    # Bootstrap counts
    f.write(f"Tổng số nhánh\t{cg_stats['n_total']}\t{cge_stats['n_total']}\t{cge_stats['n_total'] - cg_stats['n_total']}\t0\n")
    f.write(f"Số nhánh bootstrap >= 95\t{cg_stats['n_ge95']}\t{cge_stats['n_ge95']}\t{cge_stats['n_ge95'] - cg_stats['n_ge95']}\t0\n")
    f.write(f"Số nhánh bootstrap 70-94\t{cg_stats['n_70_94']}\t{cge_stats['n_70_94']}\t{cge_stats['n_70_94'] - cg_stats['n_70_94']}\t0\n")
    f.write(f"Số nhánh bootstrap < 70\t{cg_stats['n_lt70']}\t{cge_stats['n_lt70']}\t{cge_stats['n_lt70'] - cg_stats['n_lt70']}\t0\n")
    f.write(f"Bootstrap trung bình\t{cg_stats['mean']:.2f}\t{cge_stats['mean']:.2f}\t{cge_stats['mean'] - cg_stats['mean']:.2f}\t0\n")

print(f"\nĐã lưu kết quả vào: {OUTPUT_TSV}")