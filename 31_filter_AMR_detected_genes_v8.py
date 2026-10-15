#!/usr/bin/env python3
"""
FILTER AMR DETECTION RESULTS - V2.67 (FINAL)
CHANGES:
1. ADC: GOM THEO NODE PATTERN CHIẾM ĐA SỐ, chọn gene tốt nhất
2. OXA: GOM THEO NODE PATTERN CHIẾM ĐA SỐ, chọn gene tốt nhất
3. APH, msr, mph, armA, qac, ade: KHÔNG GOM (giữ nguyên tên), threshold 80%
4. aadA1: chuẩn hóa thành ant(3'')-iia
5. APH: xử lý đúng aph(3'')-IIb và aph(3')-Ia là 2 gene khác nhau
6. APH: áp dụng threshold 80% (giống armA, msrE, mphE, qac, ade)
7. Trong cùng pattern: MỖI GENE LẤY CÁC CHỈ SỐ TỐT NHẤT
8. So sánh gene: ưu tiên match %, nếu coverage chênh lệch <= 5%, ưu tiên nodes
9. Loại bỏ ADC/OXA có RGI identity < 50%
"""

import csv
import re
import sys
import os
from collections import defaultdict
from typing import Dict, List, Set, Tuple
from functools import cmp_to_key

# ==================== CONFIGURATION ====================
DEFAULT_INPUT = "/Users/shinra/gene_detection_v2/validate_database/SAMN49741230/amr_detection_results/SAMN49741230_vs_amr_database.tsv"
DEFAULT_OUTPUT = "/Users/shinra/gene_detection_v2/validate_database/SAMN49741230/filtered_detected_AMR_gene/SAMN49741230_vs_amr_database_filtered_v3.tsv"

NODE_SIMILARITY_THRESHOLD = 90.0
MATCH_PERCENTAGE_TOLERANCE = 0.5  # Tolerance cho match %
COVERAGE_NODE_TOLERANCE = 5.0  # Nếu coverage chênh lệch <= 5%, ưu tiên nodes
DEFAULT_MIN_MATCH = 95.0
SPECIAL_MIN_MATCH = 80.0  # Cho APH, armA, msrE, mphE, qac, ade
MIN_RGI_IDENTITY = 50.0   # RGI identity tối thiểu để giữ lại

def normalize_gene_name(gene_name: str) -> str:
    """
    Normalize gene name
    - ADC: KHÔNG GOM THEO TÊN (sẽ gom theo node pattern sau)
    - OXA: KHÔNG GOM THEO TÊN (sẽ gom theo node pattern sau)
    - APH: KHÔNG GOM (giữ nguyên tên, phân biệt aph(3'')-IIb và aph(3')-Ia)
    - msrE / msr(E): KHÔNG GOM (giữ nguyên tên)
    - mphE / mph(E): KHÔNG GOM (giữ nguyên tên)
    - armA: KHÔNG GOM (giữ nguyên tên)
    - qac: KHÔNG GOM (giữ nguyên tên)
    - ade: KHÔNG GOM (giữ nguyên tên)
    - aadA1: chuẩn hóa thành ant(3'')-iia
    - TET: chỉ gom tetR và tetR(B) → 'tetR'
    - TEM: GOM THEO TÊN ('blaTEM')
    - Các gene khác: GOM THEO TÊN
    """
    gene_name_clean = gene_name.strip()
    gene_name_lower = gene_name_clean.lower()
    
    # ===== CÁC GENE KHÔNG GOM =====
    # ADC - KHÔNG GOM THEO TÊN (sẽ gom theo node pattern sau)
    if 'adc' in gene_name_lower:
        return gene_name_clean
    
    # OXA - KHÔNG GOM THEO TÊN (sẽ gom theo node pattern sau)
    if 'oxa' in gene_name_lower:
        return gene_name_clean
    
    # APH - KHÔNG GOM (giữ nguyên loại và suffix)
    if 'aph' in gene_name_lower:
        # Tìm pattern: aph(3'')-IIb, aph(3')-Ia, APH(3'')-IIb, ...
        aph_match = re.search(r'aph\s*\(([^)]+)\)\s*-?\s*([a-zA-Z0-9]+)', gene_name_lower, re.IGNORECASE)
        if aph_match:
            aph_type = aph_match.group(1).strip()  # 3'' hoặc 3'
            aph_suffix = aph_match.group(2).strip().lower()  # iib hoặc ia
            # Chuẩn hóa: giữ nguyên type và suffix
            return f"aph({aph_type})-{aph_suffix}"
        # Fallback: nếu không match được pattern chuẩn
        return gene_name_clean.lower().replace(' ', '')
    
    # msr - KHÔNG GOM
    if 'msr' in gene_name_lower:
        return gene_name_clean
    
    # mph - KHÔNG GOM
    if 'mph' in gene_name_lower:
        return gene_name_clean
    
    # armA - KHÔNG GOM
    if 'arm' in gene_name_lower:
        return gene_name_clean
    
    # qac - KHÔNG GOM
    if 'qac' in gene_name_lower:
        return gene_name_clean
    
    # ade - KHÔNG GOM
    if 'ade' in gene_name_lower:
        return gene_name_clean
    
    # aadA1 - CHỈ aadA1, chuẩn hóa thành ant(3'')-iia
    if 'aad' in gene_name_lower:
        if 'aad1' in gene_name_lower or 'aada1' in gene_name_lower:
            return "ant(3'')-iia"
        else:
            return gene_name_clean
    
    # ===== CÁC GENE GOM THEO TÊN =====
    # TET - chỉ gom tetR và tetR(B)
    if 'tet' in gene_name_lower:
        if 'tetr' in gene_name_lower:
            return 'tetR'
        else:
            return gene_name_clean
    
    # TEM
    if 'tem' in gene_name_lower:
        return 'blaTEM'
    
    # Các gene khác
    match = re.match(r'^([a-zA-Z]+)', gene_name_clean)
    if match:
        base = match.group(1).lower()
        # Các gene beta-lactamase
        if base in ['bla', 'oxa', 'adc', 'tet', 'aac', 'cat', 'sul', 'dfr', 'qnr']:
            if 'bla' in gene_name_lower:
                if 'adc' in gene_name_lower:
                    return 'blaADC'
                elif 'oxa' in gene_name_lower:
                    return 'blaOXA'
                elif 'kpc' in gene_name_lower:
                    return 'blaKPC'
                elif 'ndm' in gene_name_lower:
                    return 'blaNDM'
                elif 'vim' in gene_name_lower:
                    return 'blaVIM'
                elif 'imp' in gene_name_lower:
                    return 'blaIMP'
                else:
                    return base
            return base
    
    return gene_name_clean

def parse_blast_coverage(blast_validation: str) -> float:
    if not blast_validation or blast_validation == '':
        return 0.0
    
    match = re.search(r'Coverage:\s*([0-9.]+)%', blast_validation)
    if match:
        return float(match.group(1))
    
    match = re.search(r'Identity:\s*([0-9.]+)%', blast_validation)
    if match:
        return float(match.group(1))
    
    match = re.search(r'([0-9.]+)%', blast_validation)
    if match:
        return float(match.group(1))
    
    return 0.0

def has_blast_validation(blast_validation: str) -> bool:
    if not blast_validation or blast_validation == '':
        return False
    
    blast_keywords = ['Coverage:', 'Identity:', 'BLAST', 'Total_hits:', 'Hit']
    blast_lower = blast_validation.lower()
    
    for keyword in blast_keywords:
        if keyword.lower() in blast_lower:
            return True
    
    if len(blast_validation) > 50:
        return True
    
    return False

def is_refseq(original_contig: str) -> bool:
    if not original_contig:
        return False
    return original_contig.strip().startswith('NC_')

def calculate_node_similarity(nodes1: str, nodes2: str) -> float:
    if not nodes1 or not nodes2:
        return 0.0
    
    set1 = set(nodes1.split(','))
    set2 = set(nodes2.split(','))
    
    if not set1 or not set2:
        return 0.0
    
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    
    similarity = (intersection / union) * 100
    size_similarity = min(len(set1), len(set2)) / max(len(set1), len(set2)) * 100
    
    final_similarity = (similarity * 0.7 + size_similarity * 0.3)
    
    return final_similarity

def load_amr_results(input_file: str) -> List[Dict]:
    results = []
    
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            results.append(row)
    
    print(f"📊 Loaded {len(results)} results from {input_file}")
    return results

def get_original_rgi_identity(entry: Dict) -> float:
    """Lấy original RGI identity từ entry"""
    try:
        rgi_identity_str = entry.get('original_rgi_identity', '0')
        if rgi_identity_str and rgi_identity_str != 'N/A' and rgi_identity_str != '':
            return float(rgi_identity_str)
    except:
        pass
    return 0.0

def get_rgi_cutoff_score(entry: Dict) -> int:
    """Tính điểm cho RGI cutoff: Perfect=3, Strict=2, Loose=1, None=0"""
    cutoff = entry.get('rgi_cutoff', '').strip()
    if 'perfect' in cutoff.lower():
        return 3
    elif 'strict' in cutoff.lower():
        return 2
    elif 'loose' in cutoff.lower():
        return 1
    return 0

def get_dominant_node_pattern_no_threshold(nodes: str) -> str:
    """
    Xác định node pattern chiếm đa số (KHÔNG cần threshold)
    Trả về pattern xuất hiện nhiều nhất (3 số đầu)
    """
    if not nodes:
        return ""
    
    node_list = nodes.split(',')
    if not node_list:
        return ""
    
    # Đếm tần suất các pattern 3 số đầu
    pattern_counts = defaultdict(int)
    for node in node_list:
        node = node.strip()
        if len(node) >= 3:
            pattern = node[:3]
            pattern_counts[pattern] += 1
    
    if not pattern_counts:
        first_node = node_list[0].strip()
        return first_node[:3] if len(first_node) >= 3 else first_node
    
    # Tìm pattern xuất hiện nhiều nhất
    max_pattern = max(pattern_counts, key=pattern_counts.get)
    return max_pattern

def get_entry_score(entry: Dict) -> tuple:
    """
    Hệ thống tính điểm cho ADC và OXA
    Thứ tự ưu tiên:
    1. Match % (cao nhất)
    2. COMPLETE (không BLAST) 
    3. RGI identity (cao hơn = tốt hơn)
    4. RGI cutoff (Perfect > Strict > Loose)
    5. Coverage
    6. Node count
    7. RefSeq
    """
    match_pct = float(entry.get('match_percentage', '0%').replace('%', ''))
    blast_cov = parse_blast_coverage(entry.get('blast_validation', ''))
    db_nodes = int(entry.get('db_node_count', '0'))
    has_blast = has_blast_validation(entry.get('blast_validation', ''))
    is_ref = is_refseq(entry.get('original_contig', ''))
    rgi_identity = get_original_rgi_identity(entry)
    rgi_cutoff_score = get_rgi_cutoff_score(entry)
    
    # COMPLETE (không có BLAST) được ưu tiên
    is_complete = not has_blast
    
    return (match_pct, is_complete, rgi_identity, rgi_cutoff_score, blast_cov, db_nodes, is_ref)

def compare_scores(entry1: Dict, entry2: Dict) -> int:
    """
    So sánh 2 entry, trả về:
    1: entry1 > entry2
    -1: entry1 < entry2
    0: bằng nhau
    
    Logic:
    1. Match % (cao hơn = tốt hơn)
    2. COMPLETE (True > False)
    3. RGI identity (cao hơn = tốt hơn)
    4. RGI cutoff (Perfect > Strict > Loose)
    5. Coverage
    6. Node count
    7. RefSeq (True > False)
    """
    score1 = get_entry_score(entry1)
    score2 = get_entry_score(entry2)
    
    # 1. So sánh match %
    if score1[0] != score2[0]:
        return 1 if score1[0] > score2[0] else -1
    
    # 2. So sánh COMPLETE
    if score1[1] != score2[1]:
        return 1 if score1[1] > score2[1] else -1
    
    # 3. So sánh RGI identity
    if score1[2] != score2[2]:
        return 1 if score1[2] > score2[2] else -1
    
    # 4. So sánh RGI cutoff
    if score1[3] != score2[3]:
        return 1 if score1[3] > score2[3] else -1
    
    # 5. Coverage
    if score1[4] != score2[4]:
        return 1 if score1[4] > score2[4] else -1
    
    # 6. Node count
    if score1[5] != score2[5]:
        return 1 if score1[5] > score2[5] else -1
    
    # 7. RefSeq
    if score1[6] != score2[6]:
        return 1 if score1[6] > score2[6] else -1
    
    return 0  # Bằng nhau

def calculate_weighted_score(entry: Dict) -> float:
    """
    Tính điểm tổng hợp có trọng số cho entry
    Áp dụng sau khi đã lọc match %
    
    Trọng số:
    - COMPLETE 35%
    - Coverage 30%
    - Nodes 20%
    - RGI 8%
    - Cutoff 5%
    - RefSeq 2%
    """
    cov = parse_blast_coverage(entry.get('blast_validation', ''))
    nodes = int(entry.get('db_node_count', '0'))
    rgi_identity = get_original_rgi_identity(entry)
    rgi_cutoff_score = get_rgi_cutoff_score(entry)
    is_complete = 1 if not has_blast_validation(entry.get('blast_validation', '')) else 0
    is_ref = 1 if is_refseq(entry.get('original_contig', '')) else 0
    
    # Chuẩn hóa
    cov_score = cov / 100.0
    nodes_score = min(nodes / 300.0, 1.0)  # Giả sử max 300 nodes
    rgi_score = rgi_identity / 100.0
    cutoff_score = rgi_cutoff_score / 3.0  # Perfect=3, Strict=2, Loose=1
    
    # Trọng số
    WEIGHT_COMPLETE = 0.35
    WEIGHT_COVERAGE = 0.30
    WEIGHT_NODES = 0.20
    WEIGHT_RGI = 0.08
    WEIGHT_CUTOFF = 0.05
    WEIGHT_REFSEQ = 0.02
    
    total_score = (
        is_complete * WEIGHT_COMPLETE +
        cov_score * WEIGHT_COVERAGE +
        nodes_score * WEIGHT_NODES +
        rgi_score * WEIGHT_RGI +
        cutoff_score * WEIGHT_CUTOFF +
        is_ref * WEIGHT_REFSEQ
    )
    
    return total_score

def get_best_attributes_for_gene(entries: List[Dict]) -> Dict:
    """
    Lấy các chỉ số tốt nhất của 1 gene từ tất cả entries:
    - Match % cao nhất
    - Coverage cao nhất
    - Node count cao nhất
    - RGI identity cao nhất
    - RGI cutoff tốt nhất
    """
    best_match = 0
    best_cov = 0
    best_nodes = 0
    best_rgi = 0
    best_cutoff_score = 0
    best_entry = None
    
    for entry in entries:
        match = float(entry.get('match_percentage', '0%').replace('%', ''))
        cov = parse_blast_coverage(entry.get('blast_validation', ''))
        nodes = int(entry.get('db_node_count', '0'))
        rgi = get_original_rgi_identity(entry)
        cutoff_score = get_rgi_cutoff_score(entry)
        
        if match > best_match:
            best_match = match
        if cov > best_cov:
            best_cov = cov
        if nodes > best_nodes:
            best_nodes = nodes
        if rgi > best_rgi:
            best_rgi = rgi
        if cutoff_score > best_cutoff_score:
            best_cutoff_score = cutoff_score
        
        # Lưu entry có match cao nhất để làm đại diện
        if best_entry is None or match > float(best_entry.get('match_percentage', '0%').replace('%', '')):
            best_entry = entry
    
    return {
        'best_match': best_match,
        'best_cov': best_cov,
        'best_nodes': best_nodes,
        'best_rgi': best_rgi,
        'best_cutoff_score': best_cutoff_score,
        'best_entry': best_entry
    }

def filter_best_genes_by_pattern(entries: List[Dict], gene_family: str = "ADC") -> List[Dict]:
    """
    Lọc các gene ADC hoặc OXA:
    1. Loại bỏ các gene có RGI identity quá thấp (< 50)
    2. Tìm match % cao nhất trong group
    3. CHỈ GIỮ entries có match % >= max_match - 0.5%
    4. NHÓM THEO TÊN GENE
    5. Mỗi gene lấy các chỉ số tốt nhất (match, coverage, nodes, RGI)
    6. So sánh các gene: ưu tiên match %, sau đó coverage và nodes
    """
    if not entries:
        return []
    
    if len(entries) == 1:
        return entries
    
    # BƯỚC 1: Loại bỏ RGI identity quá thấp
    valid_entries = []
    removed_low_rgi = []
    
    for entry in entries:
        rgi_identity = get_original_rgi_identity(entry)
        if rgi_identity >= MIN_RGI_IDENTITY:
            valid_entries.append(entry)
        else:
            removed_low_rgi.append(entry)
    
    if removed_low_rgi:
        print(f"\n🗑️ REMOVED {gene_family} genes with low RGI identity (< {MIN_RGI_IDENTITY}%):")
        for gene in removed_low_rgi:
            gene_name = gene.get('gene_name', 'N/A')
            rgi = get_original_rgi_identity(gene)
            print(f"   - {gene_name}: RGI={rgi:.1f}%")
    
    if not valid_entries:
        print(f"   ❌ No valid entries with RGI >= {MIN_RGI_IDENTITY}%, removing entire group")
        return []
    
    # BƯỚC 2: Tìm match % cao nhất
    max_match = 0
    for entry in valid_entries:
        match = float(entry.get('match_percentage', '0%').replace('%', ''))
        if match > max_match:
            max_match = match
    
    print(f"\n📊 Max match in group: {max_match}%")
    
    # BƯỚC 3: CHỈ GIỮ entries có match % trong tolerance
    best_match_entries = []
    for entry in valid_entries:
        match = float(entry.get('match_percentage', '0%').replace('%', ''))
        if max_match - match <= MATCH_PERCENTAGE_TOLERANCE:
            best_match_entries.append(entry)
        else:
            gene_name = entry.get('gene_name', 'N/A')
            print(f"   ⚠️ REMOVED by match: {gene_name} (match={match}%, max={max_match}%)")
    
    if not best_match_entries:
        best_match_entries = valid_entries
    
    # BƯỚC 4: NHÓM THEO TÊN GENE
    gene_groups = defaultdict(list)
    for entry in best_match_entries:
        gene_name = entry.get('gene_name', '')
        gene_groups[gene_name].append(entry)
    
    print(f"\n📊 Grouped {len(best_match_entries)} entries into {len(gene_groups)} gene groups (match >= {max_match - MATCH_PERCENTAGE_TOLERANCE}%):")
    for gene_name, gene_entries in gene_groups.items():
        print(f"   {gene_name}: {len(gene_entries)} entries")
    
    # BƯỚC 5: Mỗi gene lấy các chỉ số tốt nhất
    gene_best_attributes = []
    
    for gene_name, gene_entries in gene_groups.items():
        attrs = get_best_attributes_for_gene(gene_entries)
        gene_best_attributes.append({
            'gene_name': gene_name,
            'best_match': attrs['best_match'],
            'best_cov': attrs['best_cov'],
            'best_nodes': attrs['best_nodes'],
            'best_rgi': attrs['best_rgi'],
            'best_cutoff_score': attrs['best_cutoff_score'],
            'entry': attrs['best_entry']
        })
        
        print(f"\n   📊 {gene_name}: best match={attrs['best_match']}%, best cov={attrs['best_cov']}%, best nodes={attrs['best_nodes']}")
    
    # BƯỚC 6: So sánh các gene
    # Ưu tiên 1: Match % cao nhất
    # Ưu tiên 2: Nếu coverage chênh lệch <= 5%, ưu tiên nodes cao hơn
    # Ưu tiên 3: Nếu coverage chênh lệch > 5%, ưu tiên coverage cao hơn
    
    def compare_genes(g1, g2):
        # So sánh match %
        if g1['best_match'] != g2['best_match']:
            return -1 if g1['best_match'] > g2['best_match'] else 1
        
        # Match % bằng nhau, so sánh coverage và nodes
        cov_diff = abs(g1['best_cov'] - g2['best_cov'])
        
        if cov_diff <= COVERAGE_NODE_TOLERANCE:
            # Coverage chênh lệch <= 5%, ưu tiên nodes cao hơn
            if g1['best_nodes'] != g2['best_nodes']:
                return -1 if g1['best_nodes'] > g2['best_nodes'] else 1
        else:
            # Coverage chênh lệch > 5%, ưu tiên coverage cao hơn
            if g1['best_cov'] != g2['best_cov']:
                return -1 if g1['best_cov'] > g2['best_cov'] else 1
        
        # So sánh RGI
        if g1['best_rgi'] != g2['best_rgi']:
            return -1 if g1['best_rgi'] > g2['best_rgi'] else 1
        
        # So sánh cutoff
        if g1['best_cutoff_score'] != g2['best_cutoff_score']:
            return -1 if g1['best_cutoff_score'] > g2['best_cutoff_score'] else 1
        
        return 0
    
    # Sắp xếp
    gene_best_attributes.sort(key=cmp_to_key(compare_genes))
    
    # Chọn gene tốt nhất
    best_gene = gene_best_attributes[0]
    best_entry = best_gene['entry']
    
    print(f"\n📊 Comparing {len(gene_best_attributes)} genes (priority: match % > coverage/nodes):")
    for i, item in enumerate(gene_best_attributes, 1):
        selected = " ✅ SELECTED" if i == 1 else ""
        print(f"   {i}. {item['gene_name']}: match={item['best_match']}%, cov={item['best_cov']}%, nodes={item['best_nodes']}{selected}")
    
    print(f"\n🏆 BEST {gene_family} in this pattern: {best_gene['gene_name']}")
    print(f"   Match: {best_gene['best_match']}%, Coverage: {best_gene['best_cov']}%, Nodes: {best_gene['best_nodes']}")
    
    # BƯỚC 7: Trả về entry của gene tốt nhất
    return [best_entry]

def group_similar_genes(results: List[Dict], threshold: float = 90.0) -> List[List[Dict]]:
    """
    Group genes:
    - ADC: gom theo node pattern chiếm đa số
    - OXA: gom theo node pattern chiếm đa số
    - Các gene khác: gom theo tên hoặc node similarity
    """
    groups = []
    
    # Đầu tiên, xử lý ADC và OXA: gom theo node pattern chiếm đa số
    adc_entries = []
    oxa_entries = []
    other_entries = []
    
    for result in results:
        gene_name = result.get('gene_name', '').lower()
        if 'adc' in gene_name:
            adc_entries.append(result)
        elif 'oxa' in gene_name:
            oxa_entries.append(result)
        else:
            other_entries.append(result)
    
    print(f"📊 Found {len(adc_entries)} ADC entries, {len(oxa_entries)} OXA entries, {len(other_entries)} other entries")
    
    # Gom ADC theo node pattern chiếm đa số
    adc_groups = defaultdict(list)
    for entry in adc_entries:
        node_pattern = get_dominant_node_pattern_no_threshold(entry.get('sample_node_sequence', ''))
        if node_pattern:
            adc_groups[node_pattern].append(entry)
        else:
            adc_groups[f"unique_{entry.get('gene_name', '')}"].append(entry)
    
    print(f"📊 Grouped ADC into {len(adc_groups)} groups by dominant node pattern")
    for pattern, entries in adc_groups.items():
        gene_names = [e.get('gene_name', '') for e in entries]
        print(f"   Pattern {pattern}: {len(entries)} entries ({', '.join(gene_names)})")
    
    # Gom OXA theo node pattern chiếm đa số
    oxa_groups = defaultdict(list)
    for entry in oxa_entries:
        node_pattern = get_dominant_node_pattern_no_threshold(entry.get('sample_node_sequence', ''))
        if node_pattern:
            oxa_groups[node_pattern].append(entry)
        else:
            oxa_groups[f"unique_{entry.get('gene_name', '')}"].append(entry)
    
    print(f"📊 Grouped OXA into {len(oxa_groups)} groups by dominant node pattern")
    for pattern, entries in oxa_groups.items():
        gene_names = [e.get('gene_name', '') for e in entries]
        print(f"   Pattern {pattern}: {', '.join(gene_names)}")
    
    # Xử lý other entries với grouping thông thường
    other_groups = []
    used_other = set()
    
    for i, result1 in enumerate(other_entries):
        if i in used_other:
            continue
        
        group = [result1]
        used_other.add(i)
        
        for j, result2 in enumerate(other_entries):
            if j in used_other:
                continue
            
            nodes1 = result1.get('sample_node_sequence', '')
            nodes2 = result2.get('sample_node_sequence', '')
            
            similarity = calculate_node_similarity(nodes1, nodes2)
            
            gene1_norm = normalize_gene_name(result1.get('gene_name', ''))
            gene2_norm = normalize_gene_name(result2.get('gene_name', ''))
            
            if similarity >= threshold or gene1_norm == gene2_norm:
                group.append(result2)
                used_other.add(j)
        
        other_groups.append(group)
    
    # Gom tất cả groups lại
    all_groups = []
    
    # Thêm ADC groups (đã được filter chọn gene tốt nhất)
    for pattern, entries in adc_groups.items():
        filtered_entries = filter_best_genes_by_pattern(entries, "ADC")
        if filtered_entries:
            all_groups.append(filtered_entries)
            print(f"   ✅ ADC group {pattern} kept: {filtered_entries[0].get('gene_name', 'N/A')}")
        else:
            print(f"   ❌ ADC group {pattern} removed (no valid entries)")
    
    # Thêm OXA groups (đã được filter chọn gene tốt nhất)
    for pattern, entries in oxa_groups.items():
        filtered_entries = filter_best_genes_by_pattern(entries, "OXA")
        if filtered_entries:
            all_groups.append(filtered_entries)
            print(f"   ✅ OXA group {pattern} kept: {filtered_entries[0].get('gene_name', 'N/A')}")
        else:
            print(f"   ❌ OXA group {pattern} removed (no valid entries)")
    
    # Thêm other groups
    all_groups.extend(other_groups)
    
    print(f"\n📊 Total groups: {len(all_groups)}")
    
    return all_groups

def select_best_genes_from_group(group: List[Dict]) -> List[Dict]:
    """
    Chọn gene tốt nhất từ group:
    1. Gom tất cả entries vào 1 list
    2. So sánh trực tiếp tất cả entries với nhau, chọn entry tốt nhất
    """
    if len(group) == 1:
        return [group[0]]
    
    # BƯỚC 1: Gom tất cả entries vào 1 list (không phân biệt gene_name)
    all_entries = group[:]
    
    # BƯỚC 2: So sánh tất cả entries với nhau, chọn entry tốt nhất
    best_entry = None
    best_score = None
    
    for entry in all_entries:
        if best_entry is None or compare_scores(entry, best_entry) == 1:
            best_entry = entry
            best_score = get_entry_score(entry)
    
    # Debug
    print(f"\n🔍 Group with {len(all_entries)} entries (comparing all):")
    
    # Gom theo gene_name để hiển thị
    gene_entries_map = defaultdict(list)
    for gene in all_entries:
        gene_full_name = gene.get('gene_name', '')
        gene_entries_map[gene_full_name].append(gene)
    
    for gene_name, entries in gene_entries_map.items():
        print(f"   {gene_name}: {len(entries)} entries")
        for gene in entries:
            match_pct = gene.get('match_percentage', 'N/A')
            has_blast = has_blast_validation(gene.get('blast_validation', ''))
            status = "BLAST" if has_blast else "COMPLETE"
            blast_cov = parse_blast_coverage(gene.get('blast_validation', ''))
            db_nodes = int(gene.get('db_node_count', '0'))
            rgi_identity = get_original_rgi_identity(gene)
            rgi_cutoff = gene.get('rgi_cutoff', 'N/A')
            is_selected = "✅ SELECTED" if gene == best_entry else ""
            print(f"      - match={match_pct} RGI={rgi_identity:.1f} cutoff={rgi_cutoff} {status} cov={blast_cov:.1f}% nodes={db_nodes} {is_selected}")
    
    match_pct = best_entry.get('match_percentage', 'N/A')
    has_blast = has_blast_validation(best_entry.get('blast_validation', ''))
    status = "COMPLETE" if not has_blast else "BLAST"
    blast_cov = parse_blast_coverage(best_entry.get('blast_validation', ''))
    db_nodes = int(best_entry.get('db_node_count', '0'))
    rgi_identity = get_original_rgi_identity(best_entry)
    rgi_cutoff = best_entry.get('rgi_cutoff', 'N/A')
    
    print(f"\n   ✅ FINAL SELECTED: {best_entry.get('gene_name', 'N/A')} "
          f"(match={match_pct}, RGI={rgi_identity:.1f}, cutoff={rgi_cutoff}, {status}, cov={blast_cov:.1f}%, nodes={db_nodes})")
    
    return [best_entry]

def is_special_threshold_gene(gene_name: str) -> bool:
    """Kiểm tra gene có threshold đặc biệt 80% không (APH, armA, msrE, mphE, qac, ade)"""
    name_lower = gene_name.lower()
    return ('aph' in name_lower or      # 🆕 APH được áp dụng threshold 80%
            'msr' in name_lower or 
            'mph' in name_lower or 
            'arm' in name_lower or
            'qac' in name_lower or
            'ade' in name_lower)

def filter_by_match_percentage(unique_genes: List[Dict]) -> List[Dict]:
    """Filter genes by match percentage - special genes 80%, others 95%"""
    filtered = []
    removed = []
    
    for gene in unique_genes:
        gene_name = gene.get('gene_name', '')
        match_pct = float(gene.get('match_percentage', '0%').replace('%', ''))
        
        # Xác định threshold
        if is_special_threshold_gene(gene_name):
            threshold = SPECIAL_MIN_MATCH
        else:
            threshold = DEFAULT_MIN_MATCH
        
        if match_pct >= threshold:
            filtered.append(gene)
        else:
            removed.append((gene, threshold))
    
    if removed:
        print(f"\n🗑️ REMOVED genes with match below threshold:")
        for gene, threshold in removed:
            gene_name = gene.get('gene_name', 'N/A')
            match = gene.get('match_percentage', 'N/A')
            print(f"   - {gene_name:<35} match={match} (threshold: {threshold}%)")
    
    return filtered

def write_filtered_results(unique_genes: List[Dict], output_file: str):
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"📁 Created output directory: {output_dir}")
    
    with open(output_file, 'w', encoding='utf-8', newline='') as f:
        if not unique_genes:
            print("⚠️ No unique genes to write")
            return
        
        fieldnames = list(unique_genes[0].keys())
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        
        for gene in unique_genes:
            writer.writerow(gene)
    
    print(f"✅ Filtered results saved to: {output_file}")

def print_summary(original_count: int, unique_genes: List[Dict], groups: List[List[Dict]]):
    print("\n" + "=" * 80)
    print("FILTERING SUMMARY V2.67")
    print("=" * 80)
    print(f"Original results: {original_count} entries")
    print(f"Unique gene groups: {len(groups)}")
    print(f"Final unique genes: {len(unique_genes)}")
    print(f"Threshold: 95% (APH/armA/msrE/mphE/qac/ade: 80%)")
    print(f"Min RGI identity: {MIN_RGI_IDENTITY}%")
    print(f"Match tolerance: ±{MATCH_PERCENTAGE_TOLERANCE}%")
    print(f"Coverage-Nodes tolerance: ±{COVERAGE_NODE_TOLERANCE}% (ưu tiên nodes nếu coverage chênh lệch <= 5%)")
    
    print("\n📊 FINAL UNIQUE GENES BY FAMILY:")
    print("-" * 80)
    families = defaultdict(int)
    for gene in unique_genes:
        name = gene.get('gene_name', '')
        if 'OXA' in name:
            family = 'OXA'
        elif 'APH' in name or 'aph' in name:
            family = 'APH'
        elif 'ADC' in name:
            family = 'ADC'
        elif 'TET' in name or 'tet' in name:
            family = 'TET'
        elif 'TEM' in name:
            family = 'TEM'
        elif 'msr' in name.lower():
            family = 'MSR'
        elif 'mph' in name.lower():
            family = 'MPH'
        elif 'arm' in name.lower():
            family = 'ARM'
        elif 'qac' in name.lower():
            family = 'QAC'
        elif 'ade' in name.lower():
            family = 'ADE'
        elif 'ant' in name.lower():
            family = 'ANT'
        else:
            family = 'Other'
        families[family] += 1
    
    for family, count in sorted(families.items()):
        print(f"  {family}: {count} genes")
    
    print("\n🔍 CHECKING IMPORTANT GENES:")
    print("-" * 80)
    
    genes_to_check = [
        'OXA-23', 'OXA-66', 'OXA-98', 'OXA-91', 'OXA-259',
        'ADC-73', 'ADC-82', 'ADC-30', 'ADC-52', 'ADC-198',
        "APH(3')-Ia", "APH(3'')-Ib", "APH(3'')-IIb", 'APH(6)-Id',
        'tetR', 'tet(B)', 'tetR(B)',
        'msrE', 'mphE', 'armA',
        'sul1', 'qacEdelta1',
        'adeA', 'adeB', 'adeC', 'adeF', 'adeG', 'adeH', 'adeI', 'adeJ', 'adeK', 'adeL', 'adeN', 'adeR', 'adeS',
        "ant(3'')-iia"
    ]
    
    for gene_name in genes_to_check:
        found = [g for g in unique_genes if gene_name in g.get('gene_name', '')]
        status = '✅ FOUND' if found else '❌ NOT FOUND'
        if found:
            match = found[0].get('match_percentage', 'N/A')
            has_blast = has_blast_validation(found[0].get('blast_validation', ''))
            gene_status = "BLAST" if has_blast else "COMPLETE"
            db_nodes = found[0].get('db_node_count', 'N/A')
            blast_cov = parse_blast_coverage(found[0].get('blast_validation', ''))
            rgi_identity = get_original_rgi_identity(found[0])
            rgi_cutoff = found[0].get('rgi_cutoff', 'N/A')
            print(f"  {gene_name}: {status} (match={match}, RGI={rgi_identity:.1f}, cutoff={rgi_cutoff}, cov={blast_cov:.1f}%, nodes={db_nodes})")
        else:
            print(f"  {gene_name}: {status}")

def main():
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    else:
        input_file = DEFAULT_INPUT
    
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    else:
        output_file = DEFAULT_OUTPUT
    
    print("=" * 80)
    print("AMR RESULTS FILTERING - V2.67 (FINAL)")
    print("=" * 80)
    print(f"Input: {input_file}")
    print(f"Output: {output_file}")
    print(f"Node similarity threshold: {NODE_SIMILARITY_THRESHOLD}%")
    print(f"Threshold: 95% (APH/armA/msrE/mphE/qac/ade: 80%)")
    print(f"Min RGI identity: {MIN_RGI_IDENTITY}%")
    print(f"Match tolerance: ±{MATCH_PERCENTAGE_TOLERANCE}%")
    print(f"Coverage-Nodes tolerance: ±{COVERAGE_NODE_TOLERANCE}% (ưu tiên nodes nếu coverage chênh lệch <= 5%)")
    print("\n📝 NGUYÊN LÝ:")
    print("  1️⃣ ADC: GOM THEO NODE PATTERN CHIẾM ĐA SỐ")
    print("  2️⃣ OXA: GOM THEO NODE PATTERN CHIẾM ĐA SỐ")
    print("  3️⃣ Trong cùng pattern: MỖI GENE LẤY CÁC CHỈ SỐ TỐT NHẤT")
    print("  4️⃣ Match %: ưu tiên cao nhất (tolerance 0.5%)")
    print("  5️⃣ Nếu coverage chênh lệch <= 5%, ưu tiên nodes cao hơn")
    print("  6️⃣ Nếu coverage chênh lệch > 5%, ưu tiên coverage cao hơn")
    print("  7️⃣ APH: phân biệt aph(3'')-IIb và aph(3')-Ia (không gom chung)")
    print("  8️⃣ APH, msr, mph, armA, qac, ade: KHÔNG GOM, threshold 80%")
    print("  9️⃣ aadA1: chuẩn hóa thành ant(3'')-iia")
    print("  🔟 Loại bỏ ADC/OXA có RGI identity < 50%")
    print("=" * 80)
    print()
    
    if not os.path.exists(input_file):
        print(f"❌ Input file not found: {input_file}")
        sys.exit(1)
    
    results = load_amr_results(input_file)
    
    if not results:
        print("❌ No results loaded")
        sys.exit(1)
    
    groups = group_similar_genes(results, NODE_SIMILARITY_THRESHOLD)
    
    unique_genes = []
    for group in groups:
        if not group:
            continue
        best_genes = select_best_genes_from_group(group)
        unique_genes.extend(best_genes)
    
    unique_genes = filter_by_match_percentage(unique_genes)
    
    write_filtered_results(unique_genes, output_file)
    print_summary(len(results), unique_genes, groups)
    
    print("\n✅ Filtering completed successfully!")

if __name__ == "__main__":
    main()
