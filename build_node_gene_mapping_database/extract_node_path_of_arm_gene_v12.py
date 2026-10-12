#!/usr/bin/env python3
"""
AMR GENE NODE EXTRACTION DATABASE - VERSION 12
- Giữ tên gốc từ bakta database
- Ưu tiên dữ liệu QC từ Master Validation Database
- Xử lý trường hợp RGI reassembled bị lỗi identity
"""

import json
import csv
import subprocess
import os
import sys
import re
import gzip
from collections import defaultdict

# ========== CONFIG - CHỈ CẦN THAY ĐỔI SAMPLE_ID Ở ĐÂY ==========
SAMPLE_ID = "SAMEA8218225"  # <--- THAY ĐỔI SAMPLE_ID TẠI ĐÂY

# ========== CÁC ĐƯỜNG DẪN CỐ ĐỊNH ==========
BASE_DIR = "/Users/shinra"
OUTPUT_DB = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/Graph_annotation/amr_gene_nodes_database_v10.tsv"

# ========== DATABASE BAKTA_RGI_INTEGRATION ==========
BAKTA_RGI_DB = f"{BASE_DIR}/Bioinformatics/WGS_AMR/data_public/Graph_annotation/bakta_rgi_integration_107_samples_v5.tsv"

# ========== MASTER VALIDATION DATABASE ==========
MASTER_VALIDATION_DB = f"{BASE_DIR}/node_gene_mapping_database_validation/MASTER_VALIDATION_DATABASE_v9.tsv"

# ========== ĐƯỜNG DẪN FILE RGI RESULTS GỐC ==========
ORIGINAL_RGI_FILE = f"{BASE_DIR}/Bioinformatics/WGS_AMR/data_public/Rgi_results/{SAMPLE_ID}/rgi_results.txt"

# ========== HÀM ĐỌC DỮ LIỆU TỪ VALIDATION DATABASE ==========
def load_validation_data(sample_id):
    """Load QC data from master validation database"""
    validation_data = {}
    
    if not os.path.exists(MASTER_VALIDATION_DB):
        print(f"  ⚠️ Master validation database not found: {MASTER_VALIDATION_DB}")
        return validation_data
    
    try:
        with open(MASTER_VALIDATION_DB, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if row.get('Sample_ID') == sample_id:
                    gene_id = row.get('Gene_ID')
                    if gene_id:
                        validation_data[gene_id] = {
                            'original_identity': row.get('Original_RGI_Identity', '0'),
                            'validation_identity': row.get('Validation_Identity', '0'),
                            'qc_status': row.get('QC_Status', 'UNKNOWN'),
                            'gene_name': row.get('Gene_Name', ''),
                            'match_status': row.get('Gene_Match_Status', '')
                        }
        print(f"  📊 Loaded validation data for {len(validation_data)} genes")
    except Exception as e:
        print(f"  ⚠️ Error loading validation database: {e}")
    
    return validation_data

# ========== CÁC HÀM KHÁC GIỮ NGUYÊN ==========
def detect_gene_prefix():
    """Tự động phát hiện gene prefix từ các directories"""
    sample_dir = f"{BASE_DIR}/{SAMPLE_ID}"
    amr_gene_dir = f"{sample_dir}/amr_gene_assemblies"
    
    directories_to_check = [
        amr_gene_dir,
        sample_dir,
        f"{BASE_DIR}/Bioinformatics/WGS_AMR/data_public/GAM_files/{SAMPLE_ID}"
    ]
    
    gene_prefixes = set()
    
    for check_dir in directories_to_check:
        if not os.path.exists(check_dir):
            continue
        
        try:
            items = os.listdir(check_dir)
            for item in items:
                item_path = os.path.join(check_dir, item)
                if os.path.isdir(item_path):
                    match = re.match(r'^([A-Z]{2,})_(\d+)$', item)
                    if match:
                        prefix = match.group(1)
                        gene_prefixes.add(prefix)
        except Exception:
            continue
    
    if gene_prefixes:
        prefix_counts = {}
        for prefix in gene_prefixes:
            count = sum(1 for check_dir in directories_to_check 
                       if os.path.exists(check_dir) and 
                       any(item.startswith(prefix + '_') for item in os.listdir(check_dir)))
            prefix_counts[prefix] = count
        
        best_prefix = max(prefix_counts.items(), key=lambda x: x[1])[0]
        print(f"🔍 Auto-detected gene prefix: {best_prefix}")
        return best_prefix
    else:
        if os.path.exists(OUTPUT_DB):
            try:
                with open(OUTPUT_DB, 'r') as f:
                    for line in f:
                        if line.startswith('gene_id\t') or not line.strip():
                            continue
                        parts = line.split('\t')
                        if len(parts) > 0 and '_' in parts[0]:
                            prefix = parts[0].split('_')[0]
                            if prefix and prefix.isalpha() and prefix.isupper():
                                print(f"📁 Found existing prefix in database: {prefix}")
                                return prefix
            except:
                pass
        
        print("⚠️  Could not auto-detect gene prefix, using fallback: check all")
        return None

def get_original_contig_from_bakta_db():
    """Lấy Original_Contig từ file Bakta RGI integration database"""
    print(f"\n📊 Đang lấy Original_Contig từ Bakta DB...")
    print(f"   File: {BAKTA_RGI_DB}")
    
    if not os.path.exists(BAKTA_RGI_DB):
        print(f"❌ File Bakta DB không tồn tại: {BAKTA_RGI_DB}")
        return {}
    
    gene_to_contig = {}
    
    try:
        with open(BAKTA_RGI_DB, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        header_line = None
        for i, line in enumerate(lines):
            if 'Sample_ID' in line and 'Original_Contig' in line and 'GeneID' in line:
                header_line = line.strip()
                break
        
        if not header_line:
            print("❌ Không tìm thấy header trong Bakta DB")
            return {}
        
        headers = header_line.split('\t')
        sample_idx = headers.index('Sample_ID') if 'Sample_ID' in headers else -1
        gene_id_idx = headers.index('GeneID') if 'GeneID' in headers else -1
        contig_idx = headers.index('Original_Contig') if 'Original_Contig' in headers else -1
        gene_name_idx = headers.index('Gene_Name') if 'Gene_Name' in headers else -1
        
        if sample_idx == -1 or gene_id_idx == -1 or contig_idx == -1:
            print("❌ Thiếu cột cần thiết trong Bakta DB")
            return {}
        
        for line in lines[i+1:]:
            if not line.strip():
                continue
            
            parts = line.strip().split('\t')
            if len(parts) <= max(sample_idx, gene_id_idx, contig_idx):
                continue
            
            current_sample = parts[sample_idx].strip()
            if current_sample != SAMPLE_ID:
                continue
            
            gene_id = parts[gene_id_idx].strip()
            original_contig = parts[contig_idx].strip()
            gene_name = parts[gene_name_idx].strip() if gene_name_idx < len(parts) else ""
            
            if gene_id and original_contig and gene_id not in gene_to_contig:
                gene_to_contig[gene_id] = {
                    'contig': original_contig,
                    'gene_name': gene_name
                }
        
        print(f"✅ Đã lấy Original_Contig cho {len(gene_to_contig)} genes từ Bakta DB")
        
        if gene_to_contig:
            print("   Examples:")
            for i, (gene_id, info) in enumerate(list(gene_to_contig.items())[:3]):
                print(f"     {gene_id}: {info['gene_name']} ({info['contig']})")
        
    except Exception as e:
        print(f"❌ Lỗi đọc Bakta DB: {e}")
        import traceback
        traceback.print_exc()
        return {}
    
    return gene_to_contig

def extract_original_rgi_identity(gene_prefix, bakta_gene_info):
    """Lấy Original_RGI_Identity từ file RGI gốc bằng mapping ARO ID"""
    print(f"\n📊 Đang lấy Original_RGI_Identity từ file RGI gốc...")
    print(f"   File: {ORIGINAL_RGI_FILE}")
    
    if not os.path.exists(ORIGINAL_RGI_FILE):
        print(f"❌ File RGI results gốc không tồn tại: {ORIGINAL_RGI_FILE}")
        return {}
    
    gene_mapping = {}
    aro_to_gene = {}
    
    try:
        with open(BAKTA_RGI_DB, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        header_line = None
        for i, line in enumerate(lines):
            if 'Sample_ID' in line and 'GeneID' in line:
                header_line = line.strip()
                break
        
        if not header_line:
            print("❌ Không tìm thấy header trong Bakta DB")
            return {}
        
        headers = header_line.split('\t')
        sample_idx = headers.index('Sample_ID') if 'Sample_ID' in headers else -1
        gene_id_idx = headers.index('GeneID') if 'GeneID' in headers else -1
        contig_idx = headers.index('Original_Contig') if 'Original_Contig' in headers else -1
        aro_idx = headers.index('AMR_Gene') if 'AMR_Gene' in headers else -1
        gene_name_idx = headers.index('Gene_Name') if 'Gene_Name' in headers else -1
        
        if sample_idx == -1 or gene_id_idx == -1:
            print("❌ Thiếu cột cần thiết trong Bakta DB")
            return {}
        
        for line in lines[i+1:]:
            if not line.strip():
                continue
            
            parts = line.strip().split('\t')
            if len(parts) <= max(sample_idx, gene_id_idx, contig_idx, aro_idx):
                continue
            
            current_sample = parts[sample_idx].strip()
            if current_sample != SAMPLE_ID:
                continue
            
            gene_id = parts[gene_id_idx].strip()
            
            if gene_prefix:
                if not gene_id.startswith(gene_prefix + '_'):
                    continue
            else:
                if not re.match(r'^[A-Z]{2,}_\d+$', gene_id):
                    continue
            
            contig = parts[contig_idx].strip() if contig_idx < len(parts) else ""
            aro_id = parts[aro_idx].strip() if aro_idx < len(parts) else ""
            gene_name = parts[gene_name_idx].strip() if gene_name_idx < len(parts) else ""
            
            if contig and aro_id:
                gene_mapping[gene_id] = {
                    'contig': contig,
                    'aro_id': aro_id,
                    'gene_name': gene_name
                }
                aro_to_gene[aro_id] = gene_id
        
        print(f"   ✅ Đã tìm thấy {len(gene_mapping)} genes trong Bakta DB")
        
    except Exception as e:
        print(f"❌ Lỗi đọc Bakta DB: {e}")
        import traceback
        traceback.print_exc()
        return {}
    
    rgi_identities = {}
    
    try:
        with open(ORIGINAL_RGI_FILE, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        header_line = None
        header_idx = -1
        for i, line in enumerate(lines):
            if 'ORF_ID' in line and 'Contig' in line and 'Best_Hit_ARO' in line:
                header_line = line.strip()
                header_idx = i
                break
        
        if not header_line:
            print("❌ Không tìm thấy header trong RGI file")
            return {}
        
        headers = header_line.split('\t')
        
        try:
            contig_idx = headers.index('Contig')
            gene_name_idx = headers.index('Best_Hit_ARO')
            identity_idx = headers.index('Best_Identities')
            aro_idx = headers.index('ARO') if 'ARO' in headers else -1
        except ValueError:
            contig_idx = 1
            gene_name_idx = 8
            identity_idx = 9
            aro_idx = 10
        
        found_count = 0
        
        for line_num, line in enumerate(lines[header_idx+1:], start=header_idx+2):
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('\t')
            if len(parts) <= max(contig_idx, gene_name_idx, identity_idx, aro_idx):
                continue
            
            rgi_contig = parts[contig_idx].strip()
            rgi_gene_name = parts[gene_name_idx].strip()
            identity_str = parts[identity_idx].strip().replace('%', '')
            rgi_aro_id = parts[aro_idx].strip() if aro_idx != -1 and aro_idx < len(parts) else ""
            
            try:
                identity = float(identity_str)
            except ValueError:
                continue
            
            matched_gene_id = None
            
            if rgi_aro_id and rgi_aro_id in aro_to_gene:
                matched_gene_id = aro_to_gene[rgi_aro_id]
            
            if not matched_gene_id:
                for gene_id, info in gene_mapping.items():
                    if rgi_contig.startswith(info['contig']):
                        matched_gene_id = gene_id
                        break
            
            if matched_gene_id:
                if matched_gene_id not in rgi_identities or identity > rgi_identities[matched_gene_id]:
                    rgi_identities[matched_gene_id] = identity
                    found_count += 1
        
        print(f"✅ Đã lấy identity cho {found_count} genes từ RGI file gốc")
        
    except Exception as e:
        print(f"❌ Lỗi đọc RGI file gốc: {e}")
        import traceback
        traceback.print_exc()
        return {}
    
    return rgi_identities

def remove_sample_from_database():
    """Xóa tất cả entries của sample hiện tại khỏi database"""
    print(f"\n🗑️  Đang xóa entries cũ của sample {SAMPLE_ID} từ database...")
    
    if not os.path.exists(OUTPUT_DB):
        print(f"✅ Database không tồn tại, không cần xóa")
        return True
    
    try:
        with open(OUTPUT_DB, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        if not lines:
            print(f"✅ Database trống")
            return True
        
        header_idx = 0
        for i, line in enumerate(lines):
            if 'gene_id' in line.lower() and 'sample_id' in line.lower():
                header_idx = i
                break
        
        header = lines[header_idx]
        kept_entries = [header]
        removed_count = 0
        
        for i in range(header_idx + 1, len(lines)):
            line = lines[i].strip()
            if not line:
                kept_entries.append(line + '\n')
                continue
            
            parts = line.split('\t')
            if len(parts) >= 3:
                sample_in_line = parts[2]
                if sample_in_line == SAMPLE_ID:
                    removed_count += 1
                    continue
            
            kept_entries.append(line + '\n')
        
        with open(OUTPUT_DB, 'w', encoding='utf-8') as f:
            f.writelines(kept_entries)
        
        print(f"✅ Đã xóa {removed_count} entries của sample {SAMPLE_ID}")
        return True
        
    except Exception as e:
        print(f"❌ Lỗi khi xóa sample khỏi database: {e}")
        return False

def calculate_qc_status(original_identity, validation_identity):
    """Tính QC status dựa trên chênh lệch identity"""
    if original_identity == 0.0:
        return "NO_DATA"
    
    diff = abs(float(original_identity) - float(validation_identity))
    
    if diff <= 5.0:
        return "PASSED"
    else:
        return f"FAILED ({diff:.1f}% diff)"

def extract_nodes_fixed(gene_id, gam_base_dir):
    """Trích xuất nodes từ GAM file"""
    paths = [
        f"{gam_base_dir}/arm_contigs_mapped_to_graph/{gene_id}/{gene_id}_mapped.gam",
        f"{gam_base_dir}/amr_gene_contig_mapped_to_graph/{gene_id}/{gene_id}_mapped.gam",
        f"{gam_base_dir}/mapped_gams/{gene_id}_mapped.gam",
        f"{BASE_DIR}/{SAMPLE_ID}/{gene_id}/{gene_id}_mapped.gam",
        f"/Users/shinra/Bioinformatics/WGS_AMR/data_public/GAM_files/{SAMPLE_ID}/{gene_id}_mapped.gam",
    ]
    
    for path in paths:
        if os.path.exists(path):
            print(f"   📁 Tìm thấy GAM file: {path}")
            
            try:
                cmd = f"vg view -a '{path}' -j 2>&1"
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)
                
                if result.returncode != 0:
                    continue
                
                if not result.stdout.strip():
                    continue
                
                nodes = []
                lines = result.stdout.strip().split('\n')
                
                for i, line in enumerate(lines):
                    if not line.strip():
                        continue
                    
                    try:
                        data = json.loads(line.strip())
                        
                        if 'path' in data and 'mapping' in data['path']:
                            for mapping in data['path']['mapping']:
                                if isinstance(mapping, dict) and 'position' in mapping:
                                    position = mapping['position']
                                    if 'node_id' in position:
                                        node_id = str(position['node_id'])
                                        nodes.append(node_id)
                    
                    except json.JSONDecodeError:
                        continue
                    except Exception:
                        continue
                
                if nodes:
                    seen = set()
                    unique_nodes = []
                    for node in nodes:
                        if node not in seen:
                            seen.add(node)
                            unique_nodes.append(node)
                    
                    print(f"   ✅ Tìm thấy {len(nodes)} nodes ({len(unique_nodes)} unique)")
                    return unique_nodes
                    
            except subprocess.TimeoutExpired:
                continue
            except Exception:
                continue
    
    print(f"   ❌ Không tìm thấy GAM file hợp lệ cho {gene_id}")
    return []

def extract_nodes_bgzf(gene_id, gam_base_dir):
    """Alternative method: đọc BGZF GAM file trực tiếp"""
    paths = [
        f"{gam_base_dir}/arm_contigs_mapped_to_graph/{gene_id}/{gene_id}_mapped.gam",
        f"{gam_base_dir}/amr_gene_contig_mapped_to_graph/{gene_id}/{gene_id}_mapped.gam",
    ]
    
    for path in paths:
        if os.path.exists(path):
            try:
                nodes = []
                with gzip.open(path, 'rt', encoding='utf-8') as f:
                    for line in f:
                        if not line.strip():
                            continue
                        
                        try:
                            data = json.loads(line.strip())
                            
                            if 'path' in data and 'mapping' in data['path']:
                                for mapping in data['path']['mapping']:
                                    if ('position' in mapping and 
                                        'node_id' in mapping['position']):
                                        nodes.append(str(mapping['position']['node_id']))
                                        
                        except json.JSONDecodeError:
                            continue
                
                if nodes:
                    unique_nodes = list(set(nodes))
                    return unique_nodes
                    
            except Exception:
                continue
    
    return []

def select_best_gene(gene_id, bakta_gene_name, ab, rgi):
    """CHỌN GENE - LUÔN GIỮ TÊN GỐC TỪ BAKTA DB"""
    print(f"\n  📍 {gene_id}")
    print(f"     Bakta DB (KEPT): {bakta_gene_name}")
    if rgi:
        print(f"     RGI suggests: {rgi['gene_name']} ({rgi.get('identity', 'N/A')}%)")
    if ab:
        print(f"     Abricate suggests: {ab['gene_name']} ({ab.get('identity', 'N/A')}%)")
    
    print(f"     ✅ Keeping original name from Bakta: {bakta_gene_name}")
    
    if rgi:
        rgi['gene_name'] = bakta_gene_name
        return rgi
    elif ab:
        ab['gene_name'] = bakta_gene_name
        return ab
    else:
        return None

def main():
    print("=" * 80)
    print(f"AMR GENE NODE EXTRACTION DATABASE - VERSION 12")
    print(f"SAMPLE: {SAMPLE_ID}")
    print("Version: 12 - ƯU TIÊN DỮ LIỆU TỪ VALIDATION DATABASE")
    print("=" * 80)
    
    # ========== TỰ ĐỘNG PHÁT HIỆN GENE PREFIX ==========
    gene_prefix = detect_gene_prefix()
    
    # ========== XÓA SAMPLE CŨ KHỎI DATABASE ==========
    if not remove_sample_from_database():
        print("\n⚠️  Không thể xóa sample cũ, dừng script!")
        sys.exit(1)
    
    # ========== LẤY DỮ LIỆU TỪ VALIDATION DATABASE ==========
    validation_data = load_validation_data(SAMPLE_ID)
    
    # ========== LẤY ORIGINAL CONTIG VÀ GENE NAME TỪ BAKTA DB ==========
    bakta_gene_info = get_original_contig_from_bakta_db()
    
    # ========== LẤY ORIGINAL RGI IDENTITY ==========
    original_identities = extract_original_rgi_identity(gene_prefix, bakta_gene_info)
    
    # ========== ĐƯỜNG DẪN FILES ==========
    abricate_file = f"{BASE_DIR}/{SAMPLE_ID}/amr_gene_contig_abricate_results/{SAMPLE_ID}/{SAMPLE_ID}_abricate_results.tab"
    rgi_file = f"{BASE_DIR}/{SAMPLE_ID}/amr_gene_contig_rgi_results/{SAMPLE_ID}/{SAMPLE_ID}_rgi_results.tab"
    gam_base_dir = f"{BASE_DIR}/{SAMPLE_ID}"
    
    # ========== KIỂM TRA FILES TỒN TẠI ==========
    print("\n🔍 Kiểm tra file tồn tại...")
    for name, path in [("Abricate", abricate_file), ("RGI reassembled", rgi_file)]:
        if not os.path.exists(path):
            print(f"❌ File không tồn tại: {path}")
            sys.exit(1)
        print(f"✅ {name}: {os.path.basename(path)}")
    
    # ========== PARSE ABRICATE ==========
    print(f"\n📄 Đang đọc Abricate results...")
    abricate_results = []
    
    try:
        with open(abricate_file, 'r', encoding='utf-8', errors='ignore') as f:
            lines = [line.strip() for line in f if line.strip()]
        
        header = None
        for line in lines:
            if line.startswith('GENE_ID') and '\t' in line:
                header = line
                break
        
        if not header:
            print("❌ Không tìm thấy header trong Abricate file")
            sys.exit(1)
        
        columns = header.split('\t')
        data_start = lines.index(header) + 1
        
        for line in lines[data_start:]:
            if not line or line.startswith('#'):
                continue
            
            fields = line.split('\t')
            if len(fields) < 2 or fields[1] == '#FILE':
                continue
            
            row = {}
            for i, col in enumerate(columns):
                if i < len(fields):
                    row[col.lower()] = fields[i]
            
            gene_id = row.get('gene_id', '')
            if not gene_id:
                continue
            
            try:
                identity_str = row.get('%identity', '0').replace('%', '').strip()
                entry = {
                    'gene_id': gene_id,
                    'gene_name': row.get('gene', ''),
                    'contig': row.get('sequence', f"{gene_id}_contig00001"),
                    'start': int(row.get('start', 0)),
                    'end': int(row.get('end', 0)),
                    'identity': float(identity_str),
                    'source': 'abricate'
                }
                abricate_results.append(entry)
            except:
                continue
        
        print(f"✅ Đọc được {len(abricate_results)} genes từ Abricate")
        
    except Exception as e:
        print(f"❌ Lỗi đọc Abricate file: {e}")
        sys.exit(1)
    
    # ========== PARSE RGI (reassembled) ==========
    print(f"\n📄 Đang đọc RGI results (reassembled)...")
    rgi_results = []
    
    try:
        with open(rgi_file, 'r', encoding='utf-8', errors='ignore') as f:
            lines = [line.strip() for line in f if line.strip()]
        
        header = None
        for line in lines:
            if 'GENE_ID' in line and 'Best_Hit_ARO' in line:
                header = line
                break
        
        if not header:
            print("❌ Không tìm thấy header trong RGI file")
            sys.exit(1)
        
        columns = header.split('\t')
        data_start = lines.index(header) + 1
        
        gene_entries = defaultdict(list)
        
        for line in lines[data_start:]:
            if not line:
                continue
            
            fields = line.split('\t')
            if len(fields) < len(columns):
                continue
            
            row = {}
            for i, col in enumerate(columns):
                if i < len(fields):
                    norm_col = col.lower().replace('(%)', '').replace(' ', '_').replace('(', '').replace(')', '')
                    row[norm_col] = fields[i]
            
            gene_id = row.get('gene_id', '')
            if not gene_id:
                continue
            
            try:
                identity_str = row.get('best_identities', '0').replace('%', '').strip()
                entry = {
                    'gene_id': gene_id,
                    'gene_name': row.get('best_hit_aro', ''),
                    'contig': row.get('contig', f"{gene_id}_contig00001"),
                    'start': int(row.get('start', 0)),
                    'end': int(row.get('stop', 0)),
                    'identity': float(identity_str),
                    'cut_off': row.get('cut_off', 'Loose'),
                    'source': 'rgi'
                }
                gene_entries[gene_id].append(entry)
            except:
                continue
        
        for gene_id, entries in gene_entries.items():
            if len(entries) == 1:
                rgi_results.append(entries[0])
            else:
                sorted_entries = sorted(entries,
                    key=lambda x: (
                        {'Perfect': 3, 'Strict': 2, 'Loose': 1}.get(x.get('cut_off', 'Loose'), 0),
                        -x['identity']
                    ),
                    reverse=True
                )
                rgi_results.append(sorted_entries[0])
        
        print(f"✅ Đọc được {len(rgi_results)} genes từ RGI reassembled")
        
    except Exception as e:
        print(f"❌ Lỗi đọc RGI file: {e}")
        sys.exit(1)
    
    # ========== CHUYỂN THÀNH DICT ==========
    ab_dict = {r['gene_id']: r for r in abricate_results}
    rgi_dict = {r['gene_id']: r for r in rgi_results}
    
    # ========== CHỈ LẤY GENES CÓ TRONG BAKTA DB ==========
    all_genes = set(bakta_gene_info.keys())
    selected_genes = {}
    
    print(f"\n🎯 Đang chọn genes...")
    
    for gene_id in sorted(all_genes):
        ab = ab_dict.get(gene_id)
        rgi = rgi_dict.get(gene_id)
        
        # Lấy gene name từ Bakta DB
        bakta_gene_name = bakta_gene_info[gene_id]['gene_name'] if gene_id in bakta_gene_info else ""
        
        selected = select_best_gene(gene_id, bakta_gene_name, ab, rgi)
        if selected:
            selected_genes[gene_id] = selected
    
    print(f"✅ Đã chọn {len(selected_genes)} genes")
    
    # ========== EXTRACT NODES VÀ TẠO ENTRIES ==========
    print(f"\n🗺️  Đang trích xuất nodes từ GAM files...")
    
    entries = []
    genes_with_nodes = 0
    total_nodes = 0
    qc_stats = {'PASSED': 0, 'FAILED': 0, 'NO_DATA': 0}
    
    for i, (gene_id, gene) in enumerate(selected_genes.items()):
        if (i + 1) % 10 == 0:
            print(f"  Đang xử lý gene {i+1}/{len(selected_genes)}...")
        
        # Thử phương pháp chính
        nodes = extract_nodes_fixed(gene_id, gam_base_dir)
        
        # Nếu không được, thử phương pháp BGZF
        if not nodes:
            nodes = extract_nodes_bgzf(gene_id, gam_base_dir)
        
        node_count = len(nodes)
        
        if nodes:
            genes_with_nodes += 1
            total_nodes += node_count
        
        # Lấy Original_Contig từ mapping
        original_info = bakta_gene_info.get(gene_id, {})
        original_contig = original_info.get('contig', '') if isinstance(original_info, dict) else original_info
        original_gene_name = original_info.get('gene_name', '') if isinstance(original_info, dict) else ''
        
        if not original_contig:
            original_contig = gene['contig']
        
        # Lấy Original_RGI_Identity từ mapping
        original_identity = original_identities.get(gene_id, 0.0)
        
        # ========== QUAN TRỌNG: ƯU TIÊN DỮ LIỆU TỪ VALIDATION DATABASE ==========
        val_data = validation_data.get(gene_id, {})
        
        if val_data.get('qc_status') == 'PASSED':
            # Dùng identity từ validation thay vì từ RGI reassembled
            validation_identity = float(val_data.get('validation_identity', 0))
            qc_status = val_data.get('qc_status', 'PASSED')
            print(f"     ✅ Using validation data: validation_identity={validation_identity}%")
        else:
            # Fallback: tính từ RGI reassembled
            validation_identity = gene['identity']
            qc_status = calculate_qc_status(original_identity, validation_identity)
        
        qc_key = qc_status.split()[0] if ' ' in qc_status else qc_status
        qc_stats[qc_key] = qc_stats.get(qc_key, 0) + 1
        
        # Tạo entry
        source = gene['source']
        
        entry = {
            'gene_id': gene_id,
            'gene_name': original_gene_name if original_gene_name else gene['gene_name'],
            'sample_id': SAMPLE_ID,
            'original_contig': original_contig,
            'gene_start': gene['start'],
            'gene_end': gene['end'],
            'gene_length': gene['end'] - gene['start'] + 1,
            'node_sequence': ','.join(nodes) if nodes else 'N/A',
            'node_count': node_count,
            'original_rgi_identity': f"{original_identity:.2f}",
            'validation_identity': f"{validation_identity:.2f}",
            'qc_status': qc_status,
            'source': source,
            'cut_off': gene.get('cut_off', 'NA') if source == 'rgi' else 'NA',
            'note': 'HAS_NODES' if nodes else 'NO_NODES'
        }
        
        entries.append(entry)
    
    print(f"✅ Đã trích xuất nodes cho {genes_with_nodes}/{len(selected_genes)} genes")
    print(f"✅ Tổng số nodes: {total_nodes}")
    
    # ========== GHI VÀO DATABASE ==========
    print(f"\n📁 Đang ghi vào database: {OUTPUT_DB}")
    
    os.makedirs(os.path.dirname(OUTPUT_DB), exist_ok=True)
    
    fields = [
        'gene_id', 'gene_name', 'sample_id', 'original_contig',
        'gene_start', 'gene_end', 'gene_length', 'node_sequence',
        'node_count', 'original_rgi_identity', 'validation_identity',
        'qc_status', 'source', 'cut_off', 'note'
    ]
    
    file_exists = os.path.exists(OUTPUT_DB)
    
    with open(OUTPUT_DB, 'a' if file_exists else 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        
        if not file_exists:
            writer.writeheader()
        
        for entry in entries:
            writer.writerow(entry)
    
    print(f"✅ Đã thêm {len(entries)} entries vào database")
    
    # ========== HIỂN THỊ KẾT QUẢ ==========
    print("\n" + "=" * 80)
    print("HOÀN THÀNH! ✅")
    print("=" * 80)
    print(f"SAMPLE: {SAMPLE_ID}")
    print(f"Gene Prefix: {gene_prefix}")
    print(f"Tổng genes: {len(selected_genes)}")
    print(f"Genes có nodes: {genes_with_nodes}")
    print(f"Tổng nodes: {total_nodes}")
    print(f"Database: {OUTPUT_DB}")
    
    # Kiểm tra gene bị QC FAILED
    print("\n🔍 KIỂM TRA QC STATUS:")
    failed_genes = [e for e in entries if 'FAILED' in e['qc_status']]
    if failed_genes:
        print(f"  ⚠️ Có {len(failed_genes)} genes bị QC FAILED:")
        for gene in failed_genes[:5]:
            print(f"     - {gene['gene_id']}: {gene['qc_status']}")
    else:
        print(f"  ✅ Không có gene nào bị QC FAILED")
    
    print("\n" + "=" * 80)

if __name__ == "__main__":
    main()