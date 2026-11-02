#!/usr/bin/env python3
"""
FILTER AMR DETECTION RESULTS - V2.68 (FINAL)
CHANGES:
1. ADC: GOM THEO NODE PATTERN CHIẾM ĐA SỐ, chọn gene tốt nhất
2. OXA: GOM THEO NODE PATTERN CHIẾM ĐA SỐ, chọn gene tốt nhất
3. APH, msr, mph, armA, qac, ade, TEM: KHÔNG GOM (giữ nguyên tên), threshold 80%
4. Ưu tiên match % cao nhất TRƯỚC (chỉ giữ entries có match % cao nhất)
5. Sau đó tính điểm: Coverage (40%) > Complete (30%) > Nodes (15%) > RGI (8%) > Cutoff (5%) > RefSeq (2%)
6. Loại bỏ ADC/OXA có RGI identity < 50%
"""

import csv
import re
import sys
import os
from collections import defaultdict
from typing import Dict, List, Set, Tuple

# ==================== CONFIGURATION ====================
DEFAULT_INPUT = "/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur/amr_detection_results_filtered/AB_Precipitated_Purple_20_vs_amr_database.tsv"
DEFAULT_OUTPUT = "/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur/filtered_detected_AMR_gene/AB_Precipitated_Purple_20_vs_amr_database_filtered.tsv"

NODE_SIMILARITY_THRESHOLD = 90.0
MATCH_PERCENTAGE_TOLERANCE = 0.5
DEFAULT_MIN_MATCH = 95.0
SPECIAL_MIN_MATCH = 80.0
MIN_RGI_IDENTITY = 50.0

def normalize_gene_name(gene_name: str) -> str:
    """Normalize gene name"""
    gene_name_clean = gene_name.strip()
    gene_name_lower = gene_name_clean.lower()
    
    # KHÔNG GOM
    if 'adc' in gene_name_lower:
        return gene_name_clean
    if 'oxa' in gene_name_lower:
        return gene_name_clean
    if 'aph' in gene_name_lower:
        return gene_name_clean
    if 'msr' in gene_name_lower:
        return gene_name_clean
    if 'mph' in gene_name_lower:
        return gene_name_clean
    if 'arm' in gene_name_lower:
        return gene_name_clean
    if 'qac' in gene_name_lower:
        return gene_name_clean
    if 'ade' in gene_name_lower:
        return gene_name_clean
    if 'tem' in gene_name_lower:
        return gene_name_clean
    
    # GOM THEO TÊN
    if 'tet' in gene_name_lower:
        if 'tetr' in gene_name_lower:
            return 'tetR'
        else:
            return gene_name_clean
    
    match = re.match(r'^([a-zA-Z]+)', gene_name_clean)
    if match:
        base = match.group(1).lower()
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
    try:
        rgi_identity_str = entry.get('original_rgi_identity', '0')
        if rgi_identity_str and rgi_identity_str != 'N/A' and rgi_identity_str != '':
            return float(rgi_identity_str)
    except:
        pass
    return 0.0

def get_rgi_cutoff_score(entry: Dict) -> int:
    cutoff = entry.get('rgi_cutoff', '').strip()
    if 'perfect' in cutoff.lower():
        return 3
    elif 'strict' in cutoff.lower():
        return 2
    elif 'loose' in cutoff.lower():
        return 1
    return 0

def get_dominant_node_pattern_no_threshold(nodes: str) -> str:
    if not nodes:
        return ""
    
    node_list = nodes.split(',')
    if not node_list:
        return ""
    
    pattern_counts = defaultdict(int)
    for node in node_list:
        node = node.strip()
        if len(node) >= 3:
            pattern = node[:3]
            pattern_counts[pattern] += 1
    
    if not pattern_counts:
        first_node = node_list[0].strip()
        return first_node[:3] if len(first_node) >= 3 else first_node
    
    max_pattern = max(pattern_counts, key=pattern_counts.get)
    return max_pattern

def calculate_weighted_score(entry: Dict) -> float:
    """
    Tính điểm tổng hợp có trọng số cho entry
    1. Coverage (40%) - QUAN TRỌNG NHẤT sau match %
    2. COMPLETE (30%)
    3. Nodes (15%)
    4. RGI identity (8%)
    5. RGI cutoff (5%)
    6. RefSeq (2%)
    """
    cov = parse_blast_coverage(entry.get('blast_validation', ''))
    nodes = int(entry.get('db_node_count', '0'))
    rgi_identity = get_original_rgi_identity(entry)
    rgi_cutoff_score = get_rgi_cutoff_score(entry)
    is_complete = 1 if not has_blast_validation(entry.get('blast_validation', '')) else 0
    is_ref = 1 if is_refseq(entry.get('original_contig', '')) else 0
    
    # Chuẩn hóa
    cov_score = cov / 100.0
    nodes_score = min(nodes / 300.0, 1.0)
    rgi_score = rgi_identity / 100.0
    cutoff_score = rgi_cutoff_score / 3.0
    
    # Trọng số
    WEIGHT_COVERAGE = 0.40   # Tăng lên 40%
    WEIGHT_COMPLETE = 0.30
    WEIGHT_NODES = 0.15      # Giảm xuống 15%
    WEIGHT_RGI = 0.08
    WEIGHT_CUTOFF = 0.05
    WEIGHT_REFSEQ = 0.02
    
    total_score = (
        cov_score * WEIGHT_COVERAGE +
        is_complete * WEIGHT_COMPLETE +
        nodes_score * WEIGHT_NODES +
        rgi_score * WEIGHT_RGI +
        cutoff_score * WEIGHT_CUTOFF +
        is_ref * WEIGHT_REFSEQ
    )
    
    return total_score

def filter_best_genes_by_pattern(entries: List[Dict], gene_family: str = "ADC") -> List[Dict]:
    """
    Lọc các gene ADC hoặc OXA:
    1. Loại bỏ RGI identity < 50
    2. Tìm match % cao nhất
    3. Chỉ giữ entries có match % cao nhất
    4. Tính điểm các tiêu chí khác
    5. Chọn entry có điểm cao nhất
    """
    if not entries:
        return []
    
    if len(entries) == 1:
        return entries
    
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
    
    # Tìm match % cao nhất
    max_match = 0
    for entry in valid_entries:
        match = float(entry.get('match_percentage', '0%').replace('%', ''))
        if match > max_match:
            max_match = match
    
    print(f"\n📊 Max match in group: {max_match}%")
    
    # Chỉ giữ entries có match % cao nhất
    best_match_entries = []
    for entry in valid_entries:
        match = float(entry.get('match_percentage', '0%').replace('%', ''))
        if match == max_match:
            best_match_entries.append(entry)
        else:
            gene_name = entry.get('gene_name', 'N/A')
            print(f"   ⚠️ REMOVED by match: {gene_name} (match={match}%, max={max_match}%)")
    
    if not best_match_entries:
        best_match_entries = valid_entries
    
    # Tính điểm
    scored_entries = []
    for entry in best_match_entries:
        score = calculate_weighted_score(entry)
        scored_entries.append({
            'entry': entry,
            'score': score,
            'match': float(entry.get('match_percentage', '0%').replace('%', '')),
            'cov': parse_blast_coverage(entry.get('blast_validation', '')),
            'nodes': int(entry.get('db_node_count', '0')),
            'rgi': get_original_rgi_identity(entry),
            'cutoff': entry.get('rgi_cutoff', 'N/A'),
            'is_complete': not has_blast_validation(entry.get('blast_validation', ''))
        })
    
    # Sắp xếp theo điểm
    scored_entries.sort(key=lambda x: x['score'], reverse=True)
    
    best_entry = scored_entries[0]['entry']
    best_score = scored_entries[0]['score']
    
    # Debug
    best_match = scored_entries[0]['match']
    best_cov = scored_entries[0]['cov']
    best_nodes = scored_entries[0]['nodes']
    best_rgi = scored_entries[0]['rgi']
    best_cutoff = scored_entries[0]['cutoff']
    best_complete = "COMPLETE" if scored_entries[0]['is_complete'] else "BLAST"
    best_node_pattern = get_dominant_node_pattern_no_threshold(best_entry.get('sample_node_sequence', ''))
    
    print(f"\n📊 Scored {len(scored_entries)} entries with max match ({max_match}%)")
    print(f"\n🏆 BEST {gene_family}: {best_entry.get('gene_name', 'N/A')}")
    print(f"   Score: {best_score:.4f}")
    print(f"   Match: {best_match}%, {best_complete}, Coverage: {best_cov}%, Nodes: {best_nodes}")
    print(f"   RGI: {best_rgi:.1f}, Cutoff: {best_cutoff}, Pattern: {best_node_pattern}")
    
    print(f"\n📊 Top 5 entries:")
    for i, item in enumerate(scored_entries[:5], 1):
        gene_name = item['entry'].get('gene_name', 'N/A')
        complete_status = "COMPLETE" if item['is_complete'] else "BLAST"
        print(f"   {i}. {gene_name}: score={item['score']:.4f}, match={item['match']}%, {complete_status}, cov={item['cov']}%, nodes={item['nodes']}")
    
    kept_genes = [best_entry]
    removed_genes = []
    
    for entry in valid_entries:
        if entry == best_entry:
            continue
        removed_genes.append((entry, f"not selected"))
    
    if removed_genes:
        print(f"\n🗑️ REMOVED {gene_family} genes (only best kept):")
        for gene, reason in removed_genes:
            gene_name = gene.get('gene_name', 'N/A')
            match = gene.get('match_percentage', 'N/A')
            rgi = get_original_rgi_identity(gene)
            pattern = get_dominant_node_pattern_no_threshold(gene.get('sample_node_sequence', ''))
            print(f"   - {gene_name}: match={match}, RGI={rgi:.1f}, pattern={pattern} ({reason})")
    
    return kept_genes

def group_similar_genes(results: List[Dict], threshold: float = 90.0) -> List[List[Dict]]:
    groups = []
    
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
        print(f"   Pattern {pattern}: {', '.join(gene_names)}")
    
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
    
    all_groups = []
    
    for pattern, entries in adc_groups.items():
        filtered_entries = filter_best_genes_by_pattern(entries, "ADC")
        if filtered_entries:
            all_groups.append(filtered_entries)
    
    for pattern, entries in oxa_groups.items():
        filtered_entries = filter_best_genes_by_pattern(entries, "OXA")
        if filtered_entries:
            all_groups.append(filtered_entries)
    
    all_groups.extend(other_groups)
    
    print(f"📊 Total groups: {len(all_groups)}")
    
    return all_groups

def select_best_genes_from_group(group: List[Dict]) -> List[Dict]:
    if len(group) == 1:
        return [group[0]]
    
    all_entries = group[:]
    
    # Tìm match % cao nhất
    max_match = 0
    for entry in all_entries:
        match = float(entry.get('match_percentage', '0%').replace('%', ''))
        if match > max_match:
            max_match = match
    
    # Chỉ giữ entries có match % cao nhất
    best_match_entries = []
    for entry in all_entries:
        match = float(entry.get('match_percentage', '0%').replace('%', ''))
        if match == max_match:
            best_match_entries.append(entry)
    
    if not best_match_entries:
        best_match_entries = all_entries
    
    # Tính điểm
    best_entry = None
    best_score = -1
    
    for entry in best_match_entries:
        score = calculate_weighted_score(entry)
        if score > best_score:
            best_score = score
            best_entry = entry
    
    print(f"\n🔍 Group with {len(all_entries)} entries (keeping only max match = {max_match}%):")
    
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
            score = calculate_weighted_score(gene)
            is_selected = "✅ SELECTED" if gene == best_entry else ""
            in_tolerance = "✅" if (max_match - float(match_pct.replace('%', ''))) == 0 else "❌"
            print(f"      - match={match_pct} {in_tolerance} score={score:.4f} {status} cov={blast_cov:.1f}% nodes={db_nodes} {is_selected}")
    
    match_pct = best_entry.get('match_percentage', 'N/A')
    has_blast = has_blast_validation(best_entry.get('blast_validation', ''))
    status = "COMPLETE" if not has_blast else "BLAST"
    blast_cov = parse_blast_coverage(best_entry.get('blast_validation', ''))
    db_nodes = int(best_entry.get('db_node_count', '0'))
    
    print(f"\n   ✅ FINAL SELECTED: {best_entry.get('gene_name', 'N/A')} "
          f"(match={match_pct}, score={best_score:.4f}, {status}, cov={blast_cov:.1f}%, nodes={db_nodes})")
    
    return [best_entry]

def is_special_threshold_gene(gene_name: str) -> bool:
    name_lower = gene_name.lower()
    return ('msr' in name_lower or 
            'mph' in name_lower or 
            'arm' in name_lower or
            'qac' in name_lower or
            'ade' in name_lower or
            'tem' in name_lower)

def filter_by_match_percentage(unique_genes: List[Dict]) -> List[Dict]:
    filtered = []
    removed = []
    
    for gene in unique_genes:
        gene_name = gene.get('gene_name', '')
        match_pct = float(gene.get('match_percentage', '0%').replace('%', ''))
        
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
    print("FILTERING SUMMARY V2.68")
    print("=" * 80)
    print(f"Original results: {original_count} entries")
    print(f"Unique gene groups: {len(groups)}")
    print(f"Final unique genes: {len(unique_genes)}")
    print(f"Threshold: 95% (armA/msrE/mphE/qac/ade/tem: 80%)")
    print(f"Min RGI identity: {MIN_RGI_IDENTITY}%")
    print(f"Weights: Coverage 40%, Complete 30%, Nodes 15%, RGI 8%, Cutoff 5%, RefSeq 2%")
    print(f"Priority: Match % > Coverage > Complete > Nodes > RGI > Cutoff > RefSeq")
    
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
        "APH(3')-Ia", "APH(3'')-Ib", 'APH(6)-Id',
        'tetR', 'tet(B)', 'tetR(B)',
        'msrE', 'mphE', 'armA',
        'sul1', 'qacEdelta1',
        'adeA', 'adeB', 'adeC', 'adeF', 'adeG', 'adeH', 'adeI', 'adeJ', 'adeK', 'adeL', 'adeN', 'adeR', 'adeS',
        'TEM-1', 'TEM-12'
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
            score = calculate_weighted_score(found[0])
            print(f"  {gene_name}: {status} (score={score:.4f}, match={match}, RGI={rgi_identity:.1f}, cutoff={rgi_cutoff}, cov={blast_cov:.1f}%, nodes={db_nodes})")
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
    print("AMR RESULTS FILTERING - V2.68")
    print("=" * 80)
    print(f"Input: {input_file}")
    print(f"Output: {output_file}")
    print(f"Node similarity threshold: {NODE_SIMILARITY_THRESHOLD}%")
    print(f"Threshold: 95% (armA/msrE/mphE/qac/ade/tem: 80%)")
    print(f"Min RGI identity: {MIN_RGI_IDENTITY}%")
    print(f"Weights: Coverage 40%, Complete 30%, Nodes 15%, RGI 8%, Cutoff 5%, RefSeq 2%")
    print(f"Priority: Match % > Coverage > Complete > Nodes > RGI > Cutoff > RefSeq")
    print("\n📝 NGUYÊN LÝ:")
    print("  1️⃣ ADC: GOM THEO NODE PATTERN CHIẾM ĐA SỐ, chọn gene tốt nhất")
    print("  2️⃣ OXA: GOM THEO NODE PATTERN CHIẾM ĐA SỐ, chọn gene tốt nhất")
    print("  3️⃣ APH, msr, mph, armA, qac, ade, TEM: KHÔNG GOM (giữ nguyên tên), threshold 80%")
    print("  4️⃣ Ưu tiên match % CAO NHẤT (loại bỏ các entry có match % thấp hơn)")
    print("  5️⃣ Sau đó tính điểm: Coverage 40% > Complete 30% > Nodes 15% > RGI 8% > Cutoff 5% > RefSeq 2%")
    print("  6️⃣ Loại bỏ ADC/OXA có RGI identity < 50%")
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