#!/usr/bin/env python3
"""
AMR GENE DETECTION FROM GRAPH MAPPING - OPTIMIZED FOR LARGE FILES
Xử lý file GAM lớn bằng cách đọc streaming, không load toàn bộ vào memory
"""

import json
import csv
import subprocess
import os
import sys
import gzip
from collections import defaultdict
from typing import Dict, List, Set

# ==================== HARDCODED CONFIGURATION ====================
BASE_DIR = "/Users/shinra"

# Database paths
AMR_DB = f"{BASE_DIR}/Bioinformatics/WGS_AMR/data_public/Graph_annotation/amr_gene_nodes_database_v10_final.tsv"
NON_AMR_DB = f"{BASE_DIR}/Bioinformatics/WGS_AMR/data_public/Graph_annotation/non_amr_reference_database_v1.tsv"

# Test sample - THAY ĐỔI TẠI ĐÂY
TEST_SAMPLE_ID = "SAMN49741230"  # <--- CHANGE THIS

# Path to GAM file
GAM_FILE = f"{BASE_DIR}/gene_detection_v2/validate_database/{TEST_SAMPLE_ID}/reads_mapped_to_graph/{TEST_SAMPLE_ID}.gam"

# Threshold
NODE_MATCH_THRESHOLD = 75.0

# Output directory
OUTPUT_DIR = f"{BASE_DIR}/gene_detection_v2/validate_database/{TEST_SAMPLE_ID}/amr_detection_results"
OUTPUT_AMR_FILE = f"{OUTPUT_DIR}/{TEST_SAMPLE_ID}_vs_amr_database.tsv"
OUTPUT_NON_AMR_FILE = f"{OUTPUT_DIR}/{TEST_SAMPLE_ID}_vs_non_amr_database.tsv"
OUTPUT_BLAST_FILE = f"{OUTPUT_DIR}/{TEST_SAMPLE_ID}_blast_details.txt"

VERBOSE = True
CHUNK_SIZE = 100000  # Số dòng xử lý mỗi lần

# ==================== FUNCTIONS ====================

def log(message: str):
    if VERBOSE:
        print(f"[INFO] {message}")

def log_error(message: str):
    print(f"[ERROR] {message}", file=sys.stderr)

def load_amr_database(db_path: str) -> Dict:
    """Load AMR database"""
    database_entries = {}
    
    if not os.path.exists(db_path):
        log_error(f"AMR database not found: {db_path}")
        return {}
    
    print(f"📊 Loading AMR database from: {db_path}")
    
    with open(db_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        for row in reader:
            sample_id = row.get('sample_id', '').strip()
            gene_id = row.get('gene_id', '').strip()
            
            if not sample_id or not gene_id:
                continue
            
            node_sequence = row.get('node_sequence', '').strip()
            if not node_sequence or node_sequence == 'N/A':
                continue
            
            nodes = [n.strip() for n in node_sequence.split(',') if n.strip()]
            if not nodes:
                continue
            
            key = (sample_id, gene_id)
            
            database_entries[key] = {
                'sample_id': sample_id,
                'gene_id': gene_id,
                'gene_name': row.get('gene_name', ''),
                'original_contig': row.get('original_contig', ''),
                'gene_start': row.get('gene_start', ''),
                'gene_end': row.get('gene_end', ''),
                'gene_length': row.get('gene_length', ''),
                'node_count': row.get('node_count', ''),
                'node_path': nodes,
                'original_rgi_identity': row.get('original_rgi_identity', ''),
                'blast_validation': row.get('blast_validation', ''),
                'blast_validation_detailed': row.get('blast_validation_detailed', ''),
                'abricate_gene': row.get('Abricate_Gene', ''),
                'abricate_identity': row.get('Abricate_Identity', ''),
                'abricate_coverage': row.get('Abricate_Coverage', ''),
                'rgi_gene': row.get('RGI_Gene', ''),
                'rgi_identity': row.get('RGI_Identity', ''),
                'rgi_cutoff': row.get('RGI_Cutoff', ''),
                'validation_status': row.get('Validation_Status', '')
            }
    
    unique_samples = len(set(s[0] for s in database_entries.keys()))
    print(f"✅ Loaded {len(database_entries):,} AMR entries from {unique_samples} samples")
    return database_entries

def load_non_amr_database(db_path: str) -> Dict:
    """Load non-AMR reference database"""
    database_entries = {}
    
    if not os.path.exists(db_path):
        log_error(f"Non-AMR database not found: {db_path}")
        return {}
    
    print(f"📊 Loading non-AMR database from: {db_path}")
    
    try:
        with open(db_path, 'r', encoding='utf-8') as f:
            first_line = f.readline().strip()
            f.seek(0)
            
            if '\t' in first_line and ('gene_id' in first_line.lower() or 'sample_id' in first_line.lower()):
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    gene_id = row.get('gene_id', '').strip()
                    if not gene_id:
                        continue
                    
                    node_sequence = row.get('node_sequence', '')
                    if not node_sequence:
                        node_sequence = row.get('nodes', '')
                    if not node_sequence:
                        node_sequence = row.get('node_path', '')
                    
                    node_sequence = node_sequence.strip()
                    if not node_sequence or node_sequence == 'N/A':
                        continue
                    
                    nodes = [n.strip() for n in node_sequence.split(',') if n.strip()]
                    if not nodes:
                        continue
                    
                    key = ("non_amr_ref", gene_id)
                    database_entries[key] = {
                        'sample_id': row.get('sample_id', 'non_amr_reference'),
                        'gene_id': gene_id,
                        'gene_name': row.get('gene_name', row.get('name', gene_id)),
                        'original_contig': row.get('original_contig', ''),
                        'gene_start': row.get('gene_start', ''),
                        'gene_end': row.get('gene_end', ''),
                        'gene_length': row.get('gene_length', ''),
                        'node_count': str(len(nodes)),
                        'node_path': nodes,
                        'reference': row.get('reference', 'non_amr')
                    }
            else:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        gene_id = parts[0].strip()
                        node_sequence = parts[1].strip()
                        gene_name = parts[2].strip() if len(parts) > 2 else gene_id
                        original_contig = parts[3].strip() if len(parts) > 3 else ''
                        gene_start = parts[4].strip() if len(parts) > 4 else ''
                        gene_end = parts[5].strip() if len(parts) > 5 else ''
                        
                        nodes = [n.strip() for n in node_sequence.split(',') if n.strip()]
                        if nodes:
                            key = ("non_amr_ref", gene_id)
                            database_entries[key] = {
                                'sample_id': 'non_amr_reference',
                                'gene_id': gene_id,
                                'gene_name': gene_name,
                                'original_contig': original_contig,
                                'gene_start': gene_start,
                                'gene_end': gene_end,
                                'gene_length': str(int(gene_end) - int(gene_start) + 1) if gene_start and gene_end else '',
                                'node_count': str(len(nodes)),
                                'node_path': nodes,
                                'reference': 'non_amr'
                            }
        
        print(f"✅ Loaded {len(database_entries):,} non-AMR entries")
        
    except Exception as e:
        log_error(f"Error loading non-AMR database: {e}")
    
    return database_entries

def find_gam_file(sample_id: str) -> str:
    """Find GAM file for test sample"""
    possible_paths = [
        f"{BASE_DIR}/gene_detection_v2/validate_database/{sample_id}/reads_mapped_to_graph/{sample_id}.gam",
        f"{BASE_DIR}/{sample_id}/reads_mapped_to_graph/{sample_id}.gam",
        f"{BASE_DIR}/Bioinformatics/WGS_AMR/data_public/GAM_files/{sample_id}.gam",
        f"{BASE_DIR}/Bioinformatics/WGS_AMR/data_public/GAM_files/{sample_id}/{sample_id}.gam",
    ]
    
    for path in possible_paths:
        if os.path.exists(path):
            file_size = os.path.getsize(path) / (1024**3)
            print(f"📁 Found GAM file: {path}")
            print(f"   Size: {file_size:.2f} GB")
            return path
    
    log_error(f"GAM file not found for sample {sample_id}")
    return None

def extract_all_nodes_from_gam_optimized(gam_file: str) -> Set[str]:
    """
    Extract ALL nodes from GAM file - OPTIMIZED for large files
    Đọc streaming, không load toàn bộ file vào memory
    """
    print(f"🔍 Extracting nodes from GAM file (optimized for large files)...")
    
    nodes = set()
    line_count = 0
    error_count = 0
    
    try:
        # Thử đọc với gzip với error handling
        with gzip.open(gam_file, 'rt', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if not line.strip():
                    continue
                
                line_count += 1
                
                # Parse JSON
                try:
                    data = json.loads(line)
                    
                    # Lấy node từ path/mapping
                    if 'path' in data and 'mapping' in data['path']:
                        for mapping in data['path']['mapping']:
                            if isinstance(mapping, dict) and 'position' in mapping:
                                position = mapping['position']
                                if 'node_id' in position:
                                    nodes.add(str(position['node_id']))
                    
                    # Lấy node trực tiếp
                    elif 'node' in data:
                        nodes.add(str(data['node']))
                        
                except json.JSONDecodeError:
                    error_count += 1
                    if error_count <= 10:
                        log(f"   ⚠️ JSON decode error at line {line_count}")
                    continue
                
                # Progress report
                if line_count % 500000 == 0:
                    print(f"   Processed {line_count:,} lines, found {len(nodes):,} nodes...")
                    
                    # Giới hạn memory: nếu nodes quá lớn, không cần lưu tất cả
                    # Vì chúng ta chỉ cần so sánh với database
                    if len(nodes) > 1000000:
                        print(f"   ⚠️ Nodes count > 1,000,000. Continuing...")
            
            print(f"✅ Extracted {len(nodes):,} nodes from {line_count:,} lines")
            if error_count > 0:
                print(f"   ⚠️ {error_count} JSON errors encountered")
            return nodes
            
    except Exception as e:
        log_error(f"Failed to read GAM file: {e}")
        
        # Fallback: vg view với streaming
        print("   Trying vg view with streaming...")
        try:
            cmd = f"vg view -a '{gam_file}' 2>/dev/null"
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            
            line_count = 0
            for line in process.stdout:
                if not line.strip():
                    continue
                
                line_count += 1
                try:
                    data = json.loads(line)
                    if 'path' in data and 'mapping' in data['path']:
                        for mapping in data['path']['mapping']:
                            if 'position' in mapping and 'node_id' in mapping['position']:
                                nodes.add(str(mapping['position']['node_id']))
                except:
                    continue
                
                if line_count % 500000 == 0:
                    print(f"   Processed {line_count:,} lines, found {len(nodes):,} nodes...")
            
            process.wait()
            print(f"✅ Extracted {len(nodes):,} nodes using vg view")
            return nodes
                
        except Exception as e2:
            log_error(f"vg view also failed: {e2}")
    
    return nodes

def calculate_match_percentage(db_nodes: List[str], test_nodes: Set[str]) -> float:
    """Calculate percentage match"""
    if not db_nodes:
        return 0.0
    
    db_nodes_set = set(db_nodes)
    matching = db_nodes_set.intersection(test_nodes)
    return (len(matching) / len(db_nodes)) * 100

def get_matching_nodes(db_nodes: List[str], test_nodes: Set[str]) -> str:
    """Get comma-separated list of matching nodes"""
    db_nodes_set = set(db_nodes)
    matching = db_nodes_set.intersection(test_nodes)
    return ','.join(sorted(matching))

def debug_specific_genes(amr_database: Dict, test_nodes: Set[str], gene_names: List[str] = None):
    """
    Debug các gene cụ thể để xem match % của chúng
    """
    if gene_names is None:
        gene_names = ['msrE', 'mphE', 'armA']
    
    print(f"\n🔍 DEBUG: Checking specific genes: {', '.join(gene_names)}")
    print("-" * 60)
    
    found_genes = []
    
    for (db_sample, db_gene), info in amr_database.items():
        gene_name = info.get('gene_name', '')
        gene_name_lower = gene_name.lower()
        
        # Kiểm tra nếu gene_name chứa bất kỳ gene nào trong danh sách
        for target in gene_names:
            if target.lower() in gene_name_lower:
                match_pct = calculate_match_percentage(info['node_path'], test_nodes)
                found_genes.append({
                    'gene_name': gene_name,
                    'sample_id': db_sample,
                    'gene_id': db_gene,
                    'node_count': len(info['node_path']),
                    'match_percentage': round(match_pct, 2),
                    'is_above_threshold': match_pct >= NODE_MATCH_THRESHOLD
                })
                break
    
    if found_genes:
        print(f"\n{'Gene':<30} {'Sample':<15} {'Nodes':<8} {'Match %':<10} {'Status'}")
        print("-" * 80)
        for g in found_genes:
            status = "✅ KEPT" if g['is_above_threshold'] else "❌ BELOW THRESHOLD"
            print(f"{g['gene_name']:<30} {g['sample_id']:<15} {g['node_count']:<8} {g['match_percentage']:<10}% {status}")
    else:
        print(f"❌ No entries found for genes: {', '.join(gene_names)}")
    
    return found_genes

def write_amr_tsv(results: List[Dict], output_file: str):
    """Write AMR results"""
    
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    with open(output_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        
        writer.writerow([
            'test_sample_id',
            'database_sample_id',
            'gene_id',
            'gene_name',
            'original_contig',
            'gene_start',
            'gene_end',
            'gene_length',
            'sample_node_sequence',
            'db_node_count',
            'match_percentage',
            'original_rgi_identity',
            'blast_validation',
            'blast_validation_detailed',
            'abricate_gene',
            'abricate_identity',
            'abricate_coverage',
            'rgi_gene',
            'rgi_identity',
            'rgi_cutoff',
            'validation_status'
        ])
        
        for r in results:
            writer.writerow([
                r['test_sample_id'],
                r['database_sample_id'],
                r['gene_id'],
                r['gene_name'],
                r['original_contig'],
                r['gene_start'],
                r['gene_end'],
                r['gene_length'],
                r['sample_node_sequence'],
                r['db_node_count'],
                f"{r['match_percentage']}%",
                r.get('original_rgi_identity', ''),
                r.get('blast_validation', ''),
                r.get('blast_validation_detailed', ''),
                r.get('abricate_gene', ''),
                r.get('abricate_identity', ''),
                r.get('abricate_coverage', ''),
                r.get('rgi_gene', ''),
                r.get('rgi_identity', ''),
                r.get('rgi_cutoff', ''),
                r.get('validation_status', '')
            ])

def write_non_amr_tsv(results: List[Dict], output_file: str):
    """Write non-AMR results"""
    
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    with open(output_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        
        writer.writerow([
            'test_sample_id',
            'gene_id',
            'gene_name',
            'original_contig',
            'gene_start',
            'gene_end',
            'gene_length',
            'sample_node_sequence',
            'db_node_count',
            'match_percentage',
            'reference'
        ])
        
        for r in results:
            writer.writerow([
                r['test_sample_id'],
                r['gene_id'],
                r['gene_name'],
                r.get('original_contig', ''),
                r.get('gene_start', ''),
                r.get('gene_end', ''),
                r.get('gene_length', ''),
                r['sample_node_sequence'],
                r['db_node_count'],
                f"{r['match_percentage']}%",
                r.get('reference', 'non_amr')
            ])

def write_blast_details_file(results: List[Dict], output_file: str):
    """Write BLAST details"""
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("=" * 80 + "\n")
        f.write("BLAST VALIDATION DETAILS\n")
        f.write(f"Test Sample: {TEST_SAMPLE_ID}\n")
        f.write("=" * 80 + "\n\n")
        
        blast_count = 0
        for i, r in enumerate(results, 1):
            if r.get('blast_validation_detailed') and r['blast_validation_detailed'].strip():
                blast_count += 1
                f.write(f"[{i}] {r['database_sample_id']} - {r['gene_id']} ({r['gene_name']})\n")
                f.write(f"Match percentage: {r['match_percentage']}%\n")
                f.write(f"Validation status: {r.get('validation_status', 'N/A')}\n")
                f.write("-" * 40 + "\n")
                f.write(r['blast_validation_detailed'])
                f.write("\n" + "=" * 80 + "\n\n")
        
        if blast_count == 0:
            f.write("No BLAST validation details available.\n")

def main():
    print("=" * 80)
    print("AMR GENE DETECTION - OPTIMIZED FOR LARGE FILES")
    print("=" * 80)
    print(f"Test Sample: {TEST_SAMPLE_ID}")
    print(f"Threshold: >= {NODE_MATCH_THRESHOLD}% node match")
    print(f"Output AMR: {OUTPUT_AMR_FILE}")
    print(f"Output Non-AMR: {OUTPUT_NON_AMR_FILE}")
    print("=" * 80)
    print()
    
    # Find GAM file
    gam_file = find_gam_file(TEST_SAMPLE_ID)
    if not gam_file:
        sys.exit(1)
    
    # Extract nodes - OPTIMIZED
    print("\n📌 Step 1: Extracting nodes from GAM file (optimized)...")
    test_nodes = extract_all_nodes_from_gam_optimized(gam_file)
    
    if not test_nodes:
        log_error("No nodes extracted")
        sys.exit(1)
    
    # Show sample nodes
    node_list = list(test_nodes)[:20]
    print(f"\n   Sample nodes (first 20): {', '.join(node_list)}...")
    print(f"   Total unique nodes: {len(test_nodes):,}")
    
    # Load AMR database
    print("\n📌 Step 2: Loading AMR database...")
    amr_database = load_amr_database(AMR_DB)
    
    # Debug specific genes
    debug_specific_genes(amr_database, test_nodes, ['msrE', 'mphE', 'armA'])
    
    # Load non-AMR database
    print("\n📌 Step 3: Loading non-AMR database...")
    non_amr_database = load_non_amr_database(NON_AMR_DB)
    
    # Match AMR
    print("\n📌 Step 4: Matching against AMR database...")
    amr_results = []
    
    total_entries = len(amr_database)
    processed = 0
    
    for (db_sample, db_gene), info in amr_database.items():
        processed += 1
        if processed % 1000 == 0:
            print(f"   Processing {processed}/{total_entries} entries...")
        
        match_pct = calculate_match_percentage(info['node_path'], test_nodes)
        
        if match_pct >= NODE_MATCH_THRESHOLD:
            sample_node_sequence = get_matching_nodes(info['node_path'], test_nodes)
            
            amr_results.append({
                'test_sample_id': TEST_SAMPLE_ID,
                'database_sample_id': db_sample,
                'gene_id': db_gene,
                'gene_name': info['gene_name'],
                'original_contig': info['original_contig'],
                'gene_start': info['gene_start'],
                'gene_end': info['gene_end'],
                'gene_length': info['gene_length'],
                'sample_node_sequence': sample_node_sequence,
                'db_node_count': info['node_count'],
                'match_percentage': round(match_pct, 2),
                'original_rgi_identity': info.get('original_rgi_identity', ''),
                'blast_validation': info.get('blast_validation', ''),
                'blast_validation_detailed': info.get('blast_validation_detailed', ''),
                'abricate_gene': info.get('abricate_gene', ''),
                'abricate_identity': info.get('abricate_identity', ''),
                'abricate_coverage': info.get('abricate_coverage', ''),
                'rgi_gene': info.get('rgi_gene', ''),
                'rgi_identity': info.get('rgi_identity', ''),
                'rgi_cutoff': info.get('rgi_cutoff', ''),
                'validation_status': info.get('validation_status', '')
            })
    
    amr_results.sort(key=lambda x: x['match_percentage'], reverse=True)
    print(f"   ✅ Found {len(amr_results)} AMR matches")
    
    # Match non-AMR
    print("\n📌 Step 5: Matching against non-AMR database...")
    non_amr_results = []
    
    total_entries = len(non_amr_database)
    processed = 0
    
    for (_, gene_id), info in non_amr_database.items():
        processed += 1
        if processed % 1000 == 0:
            print(f"   Processing {processed}/{total_entries} entries...")
        
        match_pct = calculate_match_percentage(info['node_path'], test_nodes)
        
        if match_pct >= NODE_MATCH_THRESHOLD:
            sample_node_sequence = get_matching_nodes(info['node_path'], test_nodes)
            
            non_amr_results.append({
                'test_sample_id': TEST_SAMPLE_ID,
                'gene_id': gene_id,
                'gene_name': info['gene_name'],
                'original_contig': info.get('original_contig', ''),
                'gene_start': info.get('gene_start', ''),
                'gene_end': info.get('gene_end', ''),
                'gene_length': info.get('gene_length', ''),
                'sample_node_sequence': sample_node_sequence,
                'db_node_count': info['node_count'],
                'match_percentage': round(match_pct, 2),
                'reference': info.get('reference', 'non_amr')
            })
    
    non_amr_results.sort(key=lambda x: x['match_percentage'], reverse=True)
    print(f"   ✅ Found {len(non_amr_results)} non-AMR matches")
    
    # Write outputs
    print("\n📌 Step 6: Writing output files...")
    write_amr_tsv(amr_results, OUTPUT_AMR_FILE)
    write_non_amr_tsv(non_amr_results, OUTPUT_NON_AMR_FILE)
    write_blast_details_file(amr_results, OUTPUT_BLAST_FILE)
    
    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Test Sample: {TEST_SAMPLE_ID}")
    print(f"Nodes extracted: {len(test_nodes):,}")
    print(f"\n🔬 AMR matches: {len(amr_results)}")
    print(f"📚 Non-AMR matches: {len(non_amr_results)}")
    
    if amr_results:
        print(f"\n   Top 5 AMR matches:")
        for r in amr_results[:5]:
            print(f"      {r['database_sample_id']} - {r['gene_id']} ({r['match_percentage']}%)")
    
    if non_amr_results:
        print(f"\n   Top 5 Non-AMR matches:")
        for r in non_amr_results[:5]:
            print(f"      {r['gene_id']} - {r['gene_name']} ({r['match_percentage']}%)")
    
    print(f"\n📁 Output files:")
    print(f"   AMR: {OUTPUT_AMR_FILE}")
    print(f"   Non-AMR: {OUTPUT_NON_AMR_FILE}")
    print(f"   BLAST details: {OUTPUT_BLAST_FILE}")
    print("=" * 80)

if __name__ == "__main__":
    main()
