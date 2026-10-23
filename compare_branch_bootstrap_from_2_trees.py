import os
import re
from Bio import Phylo
from io import StringIO

# ============================================================
# CẤU HÌNH
# ============================================================

CORE_GENE_TREEFILE = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/iqtree/core_gene/core_gene_tree.treefile"
CORE_GENOME_TREEFILE = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/iqtree/core_genome/core_genome_tree.treefile"

OUTPUT_DIR = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/comparison"
os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_TSV = os.path.join(OUTPUT_DIR, "bootstrap_both.tsv")

REFERENCE_MAPPING = {
    'BJAB07104': 'NC_021726.1',
}

# ============================================================
# HÀM ĐỌC CÂY
# ============================================================

def read_tree(treefile):
    with open(treefile) as f:
        content = f.read()
    return Phylo.read(StringIO(content), "newick")

# ============================================================
# HÀM LẤY BOOTSTRAP CỦA CLADE NHỎ NHẤT VÀ CLADE CHA
# ============================================================

def get_bootstrap_info(tree, sample_name):
    """
    Trả về:
    - bs_closest: bootstrap của clade nhỏ nhất chứa mẫu
    - n_closest: số mẫu trong clade nhỏ nhất
    - bs_parent: bootstrap của clade cha
    - n_parent: số mẫu trong clade cha
    """
    # Bước 1: Tìm clade nhỏ nhất
    smallest = None
    min_n = float('inf')
    
    for clade in tree.find_clades():
        if clade.is_terminal():
            continue
        terminals = [c.name for c in clade.get_terminals() if c.name]
        if sample_name in terminals:
            n = len(terminals)
            if n < min_n:
                min_n = n
                smallest = clade
    
    if smallest is None:
        return None
    
    # Bước 2: Tìm clade cha
    parent = None
    for clade in tree.find_clades():
        if smallest in clade.clades:
            parent = clade
            break
    
    # Bước 3: Lấy thông tin
    bs_closest = smallest.confidence if smallest.confidence is not None else 0
    n_closest = min_n
    
    if parent is not None:
        bs_parent = parent.confidence if parent.confidence is not None else 0
        parent_terminals = [c.name for c in parent.get_terminals() if c.name]
        n_parent = len(parent_terminals)
    else:
        bs_parent = 0
        n_parent = n_closest
    
    return {
        'bs_closest': bs_closest,
        'n_closest': n_closest,
        'bs_parent': bs_parent,
        'n_parent': n_parent,
    }

# ============================================================
# ĐỌC CÂY
# ============================================================

print("=" * 100)
print("SO SÁNH BOOTSTRAP: CLADE NHỎ NHẤT VS CLADE CHA")
print("=" * 100)

tree_cg = read_tree(CORE_GENE_TREEFILE)
tree_cge = read_tree(CORE_GENOME_TREEFILE)

tree_cg.root_with_outgroup('BJAB07104')
tree_cge.root_with_outgroup('NC_021726.1')

samples_cg = [c.name for c in tree_cg.get_terminals() if c.name]
samples_cge = [c.name for c in tree_cge.get_terminals() if c.name]

samples_cg_mapped = [REFERENCE_MAPPING.get(s, s) for s in samples_cg]

common_samples = sorted(set(samples_cg_mapped) & set(samples_cge))
print(f"Số mẫu chung: {len(common_samples)}")

# ============================================================
# TÍNH TOÁN
# ============================================================

results = []

for sample in common_samples:
    if sample == 'NC_021726.1':
        cg_name = 'BJAB07104'
    else:
        cg_name = sample
    
    cg_info = get_bootstrap_info(tree_cg, cg_name)
    cge_info = get_bootstrap_info(tree_cge, sample)
    
    if cg_info is None or cge_info is None:
        continue
    
    # Diff cho clade nhỏ nhất
    diff_closest = cge_info['bs_closest'] - cg_info['bs_closest']
    # Diff cho clade cha
    diff_parent = cge_info['bs_parent'] - cg_info['bs_parent']
    
    results.append({
        'sample': sample,
        # Clade nhỏ nhất
        'bs_closest_cg': cg_info['bs_closest'],
        'n_closest_cg': cg_info['n_closest'],
        'bs_closest_cge': cge_info['bs_closest'],
        'n_closest_cge': cge_info['n_closest'],
        'diff_closest': diff_closest,
        # Clade cha
        'bs_parent_cg': cg_info['bs_parent'],
        'n_parent_cg': cg_info['n_parent'],
        'bs_parent_cge': cge_info['bs_parent'],
        'n_parent_cge': cge_info['n_parent'],
        'diff_parent': diff_parent,
    })

# Sắp xếp theo diff_closest giảm dần
results.sort(key=lambda x: x['diff_closest'], reverse=True)

# ============================================================
# IN KẾT QUẢ
# ============================================================

print(f"\n{'Mẫu':<18} "
      f"{'Closest (CG/CGE/Δ)':>30} "
      f"{'Parent (CG/CGE/Δ)':>30}")
print("-" * 85)

for r in results[:30]:
    closest_str = f"{r['bs_closest_cg']:>3.0f}/{r['bs_closest_cge']:>3.0f}/{r['diff_closest']:>+4.0f}"
    parent_str = f"{r['bs_parent_cg']:>3.0f}/{r['bs_parent_cge']:>3.0f}/{r['diff_parent']:>+4.0f}"
    print(f"{r['sample']:<18} {closest_str:>30} {parent_str:>30}")

print(f"\n... (tổng cộng {len(results)} mẫu)")

# ============================================================
# LƯU FILE TSV
# ============================================================

with open(OUTPUT_TSV, "w") as f:
    f.write("sample\t"
            "bs_closest_cg\tn_closest_cg\tbs_closest_cge\tn_closest_cge\tdiff_closest\t"
            "bs_parent_cg\tn_parent_cg\tbs_parent_cge\tn_parent_cge\tdiff_parent\n")
    for r in results:
        f.write(f"{r['sample']}\t"
                f"{r['bs_closest_cg']:.0f}\t{r['n_closest_cg']}\t"
                f"{r['bs_closest_cge']:.0f}\t{r['n_closest_cge']}\t{r['diff_closest']:+.0f}\t"
                f"{r['bs_parent_cg']:.0f}\t{r['n_parent_cg']}\t"
                f"{r['bs_parent_cge']:.0f}\t{r['n_parent_cge']}\t{r['diff_parent']:+.0f}\n")

print(f"\nĐã lưu kết quả vào: {OUTPUT_TSV}")

# ============================================================
# THỐNG KÊ
# ============================================================

n_improved_closest = sum(1 for r in results if r['diff_closest'] > 0)
n_improved_parent = sum(1 for r in results if r['diff_parent'] > 0)

print(f"\nThống kê:")
print(f"  Tổng số mẫu: {len(results)}")
print(f"  Số mẫu cải thiện (clade nhỏ nhất): {n_improved_closest}")
print(f"  Số mẫu cải thiện (clade cha): {n_improved_parent}")