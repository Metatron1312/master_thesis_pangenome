#!/usr/bin/env python3
"""
Script để xử lý và tổng hợp kết quả phát hiện gene kháng
Cho pipeline pangenome so sánh với RGI CARD
"""

import os
import re
import pandas as pd
import subprocess
from pathlib import Path
import json
import sys
import shutil
from datetime import datetime
from collections import defaultdict

class AMRGeneValidator:
    def __init__(self, base_dir):
        self.base_dir = Path(base_dir)
        self.results = {}
        self.RGI_THRESHOLD = 85
        self.block_info_cache = {}
        self.debug_mode = True
        self.aac_warnings = {}  # Lưu warning theo gene_name
        
    def normalize_gene_name(self, gene_name, source='pipeline'):
        """Chuẩn hóa tên gene để so sánh giữa RGI và pipeline"""
        if not gene_name:
            return ''
            
        gene_clean = gene_name.strip()
        gene_lower = gene_clean.lower()
        
        # ===== XỬ LÝ ĐẶC BIỆT - ƯU TIÊN CAO NHẤT =====
        
        # 1. gyrA / gyrB - xử lý đầu tiên
        if 'gyr' in gene_lower:
            # Loại bỏ prefix
            gene_lower = re.sub(r'acinetobacter\s+baumannii\s+', '', gene_lower)
            gene_lower = re.sub(r'acinetobacter\s+', '', gene_lower)
            gene_lower = re.sub(r'baumannii\s+', '', gene_lower)
            # Giữ nguyên tên
            gyr_match = re.search(r'gyr([ab])', gene_lower)
            if gyr_match:
                return f"gyr{gyr_match.group(1)}"
            return 'gyr'
        
        # 2. abaQ - xử lý sớm
        if 'abaq' in gene_lower:
            return 'abaq'
        
        # 3. msrE / msr(E) - chuẩn hóa về 'msr'
        if 'msr' in gene_lower:
            return 'msr'
        
        # 4. mphE / mph(E) - chuẩn hóa về 'mph'
        if 'mph' in gene_lower:
            return 'mph'
        
        # 5. tet genes - phân biệt tetR và tet(B)
        if 'tet' in gene_lower:
            if 'tetr' in gene_lower:
                return 'tetr'
            elif 'tet(b)' in gene_lower or 'tetb' in gene_lower:
                return 'tetb'
            else:
                return 'tet'
        
        # 6. APH genes - giữ nguyên loại và suffix
        if 'aph' in gene_lower:
            aph_match = re.search(r'aph\s*\(([^)]+)\)\s*-?\s*([a-zA-Z]+)', gene_lower, re.IGNORECASE)
            if aph_match:
                aph_type = aph_match.group(1).strip()
                aph_suffix = aph_match.group(2).strip()
                return f"aph({aph_type})-{aph_suffix.lower()}"
            return gene_lower.replace(' ', '')
        
        # 7. AAC genes - GOM CHUNG Ib và Ib8, có cảnh báo
        if 'aac' in gene_lower:
            # Gom tất cả AAC(6')-Ib và biến thể thành 'aac(6')-ib'
            if "aac(6')" in gene_lower or 'aac(6)-' in gene_lower:
                # Nếu là biến thể Ib8, ghi nhận để cảnh báo
                if 'ib8' in gene_lower:
                    self._add_aac_warning(gene_clean, 'Ib8')
                return 'aac(6)-ib'
            # Các loại AAC khác giữ nguyên
            aac_match = re.search(r'aac\s*\(([^)]+)\)\s*-?\s*([a-zA-Z0-9]+)', gene_lower, re.IGNORECASE)
            if aac_match:
                aac_type = aac_match.group(1).strip()
                aac_suffix = aac_match.group(2).strip()
                return f"aac({aac_type})-{aac_suffix.lower()}"
            return gene_lower.replace(' ', '')
        
        # 8. ANT genes - giữ nguyên loại và suffix
        if 'ant' in gene_lower:
            ant_match = re.search(r'ant\s*\(([^)]+)\)\s*-?\s*([a-zA-Z0-9]+)', gene_lower, re.IGNORECASE)
            if ant_match:
                ant_type = ant_match.group(1).strip()
                ant_suffix = ant_match.group(2).strip()
                return f"ant({ant_type})-{ant_suffix.lower()}"
            return gene_lower.replace(' ', '')
        
        # 9. OXA genes
        if 'oxa' in gene_lower:
            oxa_match = re.search(r'oxa[-_]?(\d+)', gene_lower, re.IGNORECASE)
            if oxa_match:
                oxa_num = oxa_match.group(1)
                return f"blaoxa-{oxa_num}"
            return gene_lower.replace(' ', '')
        
        # 10. ADC genes
        if 'adc' in gene_lower:
            adc_match = re.search(r'adc[-_]?(\d+)', gene_lower, re.IGNORECASE)
            if adc_match:
                adc_num = adc_match.group(1)
                return f"blaadc-{adc_num}"
            return gene_lower.replace(' ', '')
        
        # 11. TEM genes
        if 'tem' in gene_lower:
            tem_match = re.search(r'tem[-_]?(\d+)', gene_lower, re.IGNORECASE)
            if tem_match:
                tem_num = tem_match.group(1)
                return f"blatem-{tem_num}"
            return 'blatem-1'
        
        # 12. Các gene đặc biệt - bỏ prefix
        special_genes = {
            'acinetobacter baumannii aba': 'aba',
            'acinetobacter baumannii abe': 'abe',
            'acinetobacter baumannii ade': 'ade',
            'acinetobacter baumannii amv': 'amv',
            'acinetobacter baumannii par': 'par',
            'acinetobacter baumannii lps': 'lps',
            'acinetobacter baumannii rsm': 'rsm',
            'acinetobacter baumannii arm': 'arm',
        }
        
        for prefix, replacement in special_genes.items():
            if prefix in gene_lower:
                gene_lower = gene_lower.replace(prefix, replacement)
        
        # Xóa các cụm thừa
        remove_phrases = [
            'conferring resistance to fluoroquinolones',
            'acinetobacter baumannii ',
            'acinetobacter ',
            'baumannii ',
            'resistance to ',
        ]
        for phrase in remove_phrases:
            gene_lower = gene_lower.replace(phrase, '')
        
        # Xóa khoảng trắng và ký tự đặc biệt
        gene_lower = re.sub(r'[^a-z0-9()\-]', '', gene_lower)
        
        return gene_lower.strip()
    
    def _add_aac_warning(self, gene_name, variant):
        """Ghi nhận cảnh báo cho biến thể AAC"""
        warning = f"AAC variant: {gene_name} is variant {variant} of AAC(6')-Ib"
        self.aac_warnings[gene_name] = warning
        print(f"\n⚠️ {warning}")

    def parse_rgi_results(self, rgi_file):
        """Đọc file RGI và lọc gene với identity > threshold"""
        rgi_genes = []
        try:
            with open(rgi_file, 'r') as f:
                lines = f.readlines()
            
            if len(lines) < 2:
                print("❌ File RGI rỗng hoặc không có dữ liệu")
                return rgi_genes
            
            header = lines[0].strip().split('\t')
            
            best_identity_col = None
            gene_name_col = None
            
            for i, col in enumerate(header):
                if col == 'Best_Identities':
                    best_identity_col = i
                elif col == 'Best_Hit_ARO':
                    gene_name_col = i
            
            if best_identity_col is None:
                print("❌ Không tìm thấy cột 'Best_Identities' trong file RGI")
                return rgi_genes
            
            for line in lines[1:]:
                if not line.strip():
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) <= best_identity_col:
                    continue
                
                identity_str = parts[best_identity_col].strip()
                try:
                    identity = float(identity_str)
                except ValueError:
                    continue
                
                if gene_name_col is not None and len(parts) > gene_name_col:
                    gene_name = parts[gene_name_col].strip()
                else:
                    gene_name = parts[0].strip().split('#')[0].strip()
                
                if not gene_name:
                    continue
                
                if identity > self.RGI_THRESHOLD:
                    normalized = self.normalize_gene_name(gene_name, source='rgi')
                    rgi_genes.append({
                        'gene_name': gene_name,
                        'identity': identity,
                        'normalized': normalized
                    })
                        
        except Exception as e:
            print(f"❌ Lỗi khi đọc file RGI: {e}")
            
        return rgi_genes

    def parse_pipeline_results(self, pipeline_file):
        """Đọc file kết quả từ pipeline"""
        pipeline_genes = []
        try:
            df = pd.read_csv(pipeline_file, sep='\t')
            
            for _, row in df.iterrows():
                gene_name = row['gene_name']
                normalized = self.normalize_gene_name(gene_name, source='pipeline')
                pipeline_genes.append({
                    'gene_name': gene_name,
                    'gene_id': row.get('gene_id', ''),
                    'normalized': normalized,
                    'match_percentage': row.get('match_percentage', None),
                    'contig': row.get('original_contig', None),
                    'sample_node_sequence': row.get('sample_node_sequence', '')
                })
                
        except Exception as e:
            print(f"❌ Lỗi khi đọc file pipeline: {e}")
            
        return pipeline_genes

    def compare_results(self, rgi_file, pipeline_file):
        """So sánh kết quả giữa RGI và pipeline"""
        print("\n📖 Đang đọc file RGI...")
        self.aac_warnings = {}  # Reset warnings
        rgi_genes = self.parse_rgi_results(rgi_file)
        print(f"   ✅ Đọc được {len(rgi_genes)} gene từ RGI (identity > {self.RGI_THRESHOLD}%)")
        
        print("\n📖 Đang đọc file pipeline...")
        pipeline_genes = self.parse_pipeline_results(pipeline_file)
        print(f"   ✅ Đọc được {len(pipeline_genes)} gene từ pipeline")
        
        rgi_dict = {g['normalized']: g for g in rgi_genes}
        pipeline_dict = {g['normalized']: g for g in pipeline_genes}
        
        # In cảnh báo AAC variants nếu có
        if self.aac_warnings:
            print("\n📋 AAC VARIANT WARNINGS:")
            for gene_name, warning in self.aac_warnings.items():
                print(f"   ⚠️ {warning}")
        
        # DEBUG
        if self.debug_mode:
            print("\n🔍 DEBUG - NORMALIZED NAMES:")
            print("-" * 60)
            
            # In tất cả normalized names để kiểm tra
            print("\n📋 RGI normalized names:")
            for norm, info in sorted(rgi_dict.items()):
                print(f"   '{norm}' ← {info['gene_name']}")
            
            print("\n📋 Pipeline normalized names:")
            for norm, info in sorted(pipeline_dict.items()):
                print(f"   '{norm}' ← {info['gene_name']}")
            
            # Tìm gene RGI không có trong pipeline
            rgi_only = set(rgi_dict.keys()) - set(pipeline_dict.keys())
            if rgi_only:
                print("\n⚠️ RGI-ONLY GENES (pipeline KHÔNG phát hiện):")
                for norm in sorted(rgi_only):
                    info = rgi_dict[norm]
                    print(f"   - '{norm}' → {info['gene_name']} ({info['identity']:.2f}%)")
            
            # Tìm gene pipeline không có trong RGI
            pipeline_only = set(pipeline_dict.keys()) - set(rgi_dict.keys())
            if pipeline_only:
                print("\n⚠️ PIPELINE-ONLY GENES (RGI KHÔNG phát hiện):")
                for norm in sorted(pipeline_only):
                    info = pipeline_dict[norm]
                    print(f"   - '{norm}' → {info['gene_name']} ({info.get('match_percentage', 'N/A')})")
        
        common_genes = set(rgi_dict.keys()) & set(pipeline_dict.keys())
        pipeline_only = set(pipeline_dict.keys()) - set(rgi_dict.keys())
        rgi_only = set(rgi_dict.keys()) - set(pipeline_dict.keys())
        
        return {
            'rgi_genes': rgi_genes,
            'pipeline_genes': pipeline_genes,
            'common_genes': common_genes,
            'pipeline_only': pipeline_only,
            'rgi_only': rgi_only,
            'rgi_count': len(rgi_genes),
            'pipeline_count': len(pipeline_genes),
            'common_count': len(common_genes)
        }

    def parse_blast_report(self, blast_report_file):
        """Parse file báo cáo BLAST"""
        blast_info = {
            'identity': None,
            'coverage': None,
            'evalue': None,
            'contig': None,
            'found': False,
            'gene_name': None,
            'num_blocks': 1,
            'blocks': []
        }
        
        try:
            with open(blast_report_file, 'r') as f:
                content = f.read()
                
            gene_name_match = re.search(r'/([^/]+)_report\.txt$', str(blast_report_file))
            if gene_name_match:
                blast_info['gene_name'] = gene_name_match.group(1)
            
            identity_match = re.search(r'Identity:\s*([\d.]+)%', content)
            if identity_match:
                blast_info['identity'] = float(identity_match.group(1))
                
            coverage_match = re.search(r'Coverage:\s*([\d.]+)%', content)
            if coverage_match:
                blast_info['coverage'] = float(coverage_match.group(1))
                
            evalue_match = re.search(r'E-value:\s*([\d.eE+-]+)', content)
            if evalue_match:
                blast_info['evalue'] = evalue_match.group(1)
                
            contig_match = re.search(r'Contig:\s*(\S+)', content)
            if contig_match:
                blast_info['contig'] = contig_match.group(1)
            
            total_hits_match = re.search(r'Tổng số hits:\s*(\d+)', content)
            if total_hits_match:
                blast_info['num_blocks'] = int(total_hits_match.group(1))
            
            position_matches = re.findall(r'(\d+)\.\s+(\S+):\s+([\d.]+)% identity,\s+([\d.]+)% coverage,\s+position\s+(\d+)-(\d+)\s+\(([+-])\)', content)
            if position_matches:
                for match in position_matches:
                    hit_num, contig_name, identity_val, coverage_val, start, end, strand = match
                    blast_info['blocks'].append({
                        'block_num': int(hit_num),
                        'contig': contig_name,
                        'identity': float(identity_val),
                        'coverage': float(coverage_val),
                        'start': int(start),
                        'end': int(end),
                        'strand': strand
                    })
            
            if '✅' in content and 'phát hiện' in content:
                blast_info['found'] = True
                
        except Exception as e:
            print(f"Lỗi khi parse BLAST report {blast_report_file}: {e}")
            
        return blast_info

    def get_average_depth(self, bam_file, contig_name, region_start=None, region_end=None):
        """Lấy độ phủ trung bình của contig hoặc region từ file BAM"""
        try:
            if region_start is not None and region_end is not None:
                region = f"{contig_name}:{region_start}-{region_end}"
            else:
                region = contig_name
                
            cmd = f"samtools depth -r '{region}' {bam_file} 2>/dev/null | awk '{{sum+=$3; count++}} END {{if(count>0) print sum/count; else print 0}}'"
            
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode == 0 and result.stdout.strip():
                depth = float(result.stdout.strip())
                return depth
            else:
                return 0.0
                
        except Exception as e:
            print(f"❌ Lỗi khi lấy depth: {e}")
            
        return 0.0

    def get_average_depth_from_gam(self, sample_gam, node_sequence, gene_name):
        """
        Lấy depth trung bình từ GAM file cho gene không có trong reference
        Sử dụng node sequence từ pipeline result
        """
        print(f"\n   📊 Tính depth từ GAM file (gene không có trong reference)...")
        
        if not node_sequence:
            print(f"      ⚠️ Không có node sequence cho {gene_name}")
            return 0.0
        
        nodes = [n.strip() for n in node_sequence.split(',') if n.strip()]
        if not nodes:
            print(f"      ⚠️ Node sequence rỗng cho {gene_name}")
            return 0.0
        
        print(f"      Số nodes: {len(nodes)}")
        print(f"      Node đầu: {nodes[0]}")
        print(f"      Node cuối: {nodes[-1]}")
        
        node_set = set(nodes)
        node_counts = defaultdict(int)
        processed_lines = 0
        total_reads = 0
        
        try:
            # Dùng vg view để đọc GAM
            cmd = f"vg view -a '{sample_gam}' -j 2>/dev/null"
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            
            for line in process.stdout:
                if not line.strip():
                    continue
                processed_lines += 1
                if processed_lines % 100000 == 0:
                    print(f"      Đã xử lý {processed_lines:,} reads...")
                
                try:
                    data = json.loads(line)
                    if 'path' in data and 'mapping' in data['path']:
                        for mapping in data['path']['mapping']:
                            if 'position' in mapping and 'node_id' in mapping['position']:
                                node_id = str(mapping['position']['node_id'])
                                if node_id in node_set:
                                    node_counts[node_id] += 1
                                    total_reads += 1
                except:
                    continue
            
            process.terminate()
            
        except Exception as e:
            print(f"      ❌ Lỗi khi đọc GAM: {e}")
            return 0.0
        
        # Tính depth trung bình
        num_nodes = len(nodes)
        if num_nodes == 0:
            return 0.0
        
        avg_depth = total_reads / num_nodes
        
        print(f"\n      📊 Kết quả depth từ GAM:")
        print(f"         Tổng số reads: {total_reads}")
        print(f"         Số nodes: {num_nodes}")
        print(f"         Depth trung bình: {avg_depth:.2f} reads/node")
        print(f"         Số node có reads: {len(node_counts)}/{num_nodes}")
        
        # Hiển thị phân phối depth
        if node_counts:
            depths = [node_counts.get(node, 0) for node in nodes]
            max_depth = max(depths) if depths else 0
            min_depth = min([d for d in depths if d > 0]) if any(d > 0 for d in depths) else 0
            
            print(f"         Depth cao nhất: {max_depth} reads")
            print(f"         Depth thấp nhất (có reads): {min_depth} reads")
            
            # Phân phối
            ranges = {'1-5 reads': 0, '6-10 reads': 0, '11-20 reads': 0, '>20 reads': 0}
            for d in depths:
                if d == 0:
                    continue
                elif d <= 5:
                    ranges['1-5 reads'] += 1
                elif d <= 10:
                    ranges['6-10 reads'] += 1
                elif d <= 20:
                    ranges['11-20 reads'] += 1
                else:
                    ranges['>20 reads'] += 1
            
            if any(ranges.values()):
                print(f"         Phân phối depth:")
                for range_name, count in ranges.items():
                    if count > 0:
                        print(f"            {range_name}: {count} nodes")
        
        return avg_depth

    def get_node_sequence_for_gene(self, gene_name, pipeline_info, pipeline_only_list, pipeline_file):
        """Lấy node sequence cho gene từ nhiều nguồn"""
        node_sequence = ''
        
        # 1. Từ pipeline_info
        if pipeline_info and pipeline_info.get('sample_node_sequence'):
            node_seq = pipeline_info.get('sample_node_sequence', '')
            if node_seq and node_seq != '' and node_seq != 'N/A':
                node_sequence = node_seq
                print(f"      📍 Lấy node sequence từ pipeline_info")
                return node_sequence
        
        # 2. Từ pipeline_only_list
        if pipeline_only_list:
            for item in pipeline_only_list:
                if item.get('gene_name') == gene_name:
                    node_seq = item.get('sample_node_sequence', '')
                    if node_seq and node_seq != '' and node_seq != 'N/A':
                        node_sequence = node_seq
                        print(f"      📍 Lấy node sequence từ pipeline_only_list")
                        return node_sequence
        
        # 3. Từ pipeline file filtered (TSV)
        if pipeline_file and os.path.exists(pipeline_file):
            try:
                df = pd.read_csv(pipeline_file, sep='\t')
                gene_rows = df[df['gene_name'] == gene_name]
                if not gene_rows.empty:
                    for _, row in gene_rows.iterrows():
                        for col in ['sample_node_sequence', 'node_sequence', 'sample_node_seq']:
                            if col in row.index:
                                node_seq = str(row[col]) if row[col] is not None else ''
                                if node_seq and node_seq != '' and node_seq != 'nan' and node_seq != 'N/A':
                                    node_sequence = node_seq
                                    print(f"      📍 Lấy node sequence từ pipeline file (cột: {col})")
                                    return node_sequence
            except Exception as e:
                print(f"      ⚠️ Lỗi đọc pipeline file: {e}")
        
        # 4. Từ AMR detection raw result
        if pipeline_file:
            try:
                sample_id = Path(pipeline_file).parent.parent.name
                amr_file = Path(pipeline_file).parent.parent / 'amr_detection_results' / f"{sample_id}_vs_amr_database.tsv"
                
                if os.path.exists(amr_file):
                    cmd = f"grep -m 1 '{gene_name}' {amr_file} | cut -f9"
                    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                    if result.stdout:
                        node_seq = result.stdout.strip()
                        if node_seq and node_seq != '' and node_seq != 'N/A':
                            node_sequence = node_seq
                            print(f"      📍 Lấy node sequence từ AMR result")
                            return node_sequence
            except Exception as e:
                print(f"      ⚠️ Lỗi đọc AMR result: {e}")
        
        # 5. Từ database
        db_path = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/Graph_annotation/amr_gene_nodes_database_v10_final.tsv"
        try:
            cmd = f"grep -m 1 '{gene_name}' {db_path} | cut -f8"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if result.stdout:
                node_seq = result.stdout.strip()
                if node_seq and node_seq != '' and node_seq != 'N/A':
                    node_sequence = node_seq
                    print(f"      📍 Lấy node sequence từ database")
                    return node_sequence
        except Exception as e:
            print(f"      ⚠️ Lỗi đọc database: {e}")
        
        print(f"      ⚠️ Không tìm thấy node sequence cho {gene_name}")
        return node_sequence

    def find_blast_folders(self, pipeline_only_dir):
        """Tìm tất cả các thư mục có dạng *_blast trong pipeline_only_dir"""
        blast_folders = []
        
        if not pipeline_only_dir or not os.path.exists(pipeline_only_dir):
            print(f"⚠️ Thư mục không tồn tại: {pipeline_only_dir}")
            return blast_folders
            
        try:
            for item in os.listdir(pipeline_only_dir):
                item_path = os.path.join(pipeline_only_dir, item)
                if os.path.isdir(item_path) and item.endswith('_blast'):
                    blast_folders.append({
                        'folder_name': item,
                        'folder_path': item_path,
                        'gene_name': item.replace('_blast', '')
                    })
        except Exception as e:
            print(f"Lỗi khi tìm thư mục blast: {e}")
            
        return blast_folders

    def clear_sample_output(self, sample_output_dir):
        """Xóa toàn bộ output cũ của sample trước khi ghi mới"""
        if os.path.exists(sample_output_dir):
            print(f"   🗑️ Xóa kết quả cũ: {sample_output_dir}")
            shutil.rmtree(sample_output_dir)
        os.makedirs(sample_output_dir, exist_ok=True)

    def get_block_info_from_user(self, gene_name, blast_info, pipeline_info=None, pipeline_only_list=None, pipeline_file=None):
        """
        Lấy thông tin block từ người dùng (CHỈ CHO PIPELINE-ONLY GENES)
        RGI-only genes sẽ không gọi hàm này
        """
        if gene_name in self.block_info_cache:
            print(f"\n   📦 Sử dụng thông tin block đã nhập cho gene {gene_name}")
            return self.block_info_cache[gene_name]
        
        print(f"\n   📋 THÔNG TIN TỪ BLAST:")
        print(f"      Identity: {blast_info.get('identity', 'N/A')}%")
        print(f"      Coverage: {blast_info.get('coverage', 'N/A')}%")
        print(f"      Contig: {blast_info.get('contig', 'N/A')}")
        
        print(f"\n   ❓ Gene này có bao nhiêu block?")
        print(f"      (Enter để dùng số block từ BLAST: {blast_info.get('num_blocks', 1)})")
        print(f"      (Nhập 0 nếu gene không có trong reference - sẽ tính depth từ GAM)")
        num_blocks_input = input(f"   Số block: ").strip()
        
        if num_blocks_input:
            try:
                num_blocks = int(num_blocks_input)
            except:
                print(f"      ⚠️ Số không hợp lệ, dùng số từ BLAST")
                num_blocks = blast_info.get('num_blocks', 1)
        else:
            num_blocks = blast_info.get('num_blocks', 1)
        
        # Nếu num_blocks = 0, gene không có trong reference
        if num_blocks == 0:
            print(f"\n   📌 Gene {gene_name} không có trong reference")
            print(f"   Sẽ tính depth từ GAM file (dựa trên node sequence)")
            
            # Lấy node sequence từ nhiều nguồn
            node_sequence = self.get_node_sequence_for_gene(
                gene_name, 
                pipeline_info, 
                pipeline_only_list, 
                pipeline_file
            )
            
            self.block_info_cache[gene_name] = {
                'num_blocks': 0,
                'blocks': [],
                'use_gam': True,
                'node_sequence': node_sequence
            }
            
            return self.block_info_cache[gene_name]
        
        print(f"   ✅ Số block: {num_blocks}")
        
        blocks = []
        print(f"\n   📝 Nhập thông tin cho {num_blocks} block(s):")
        print(f"   (Định dạng contig: tên_contig, Position: start-end)")
        
        for i in range(1, num_blocks + 1):
            print(f"\n   🔹 Block {i}:")
            print(f"      Nhập tên contig (ví dụ: NC_021726.1 hoặc block_xxx):")
            contig = input(f"      Contig block {i}: ").strip()
            
            print(f"      Nhập vị trí trên contig (start-end, ví dụ: 1847888-1848470):")
            position_input = input(f"      Position block {i}: ").strip()
            
            start = None
            end = None
            if position_input and '-' in position_input:
                try:
                    start_str, end_str = position_input.split('-')
                    start = int(start_str.strip())
                    end = int(end_str.strip())
                except:
                    print(f"      ⚠️ Format không hợp lệ, để trống để lấy toàn bộ contig")
            
            blocks.append({
                'block_num': i,
                'contig': contig,
                'start': start,
                'end': end
            })
        
        self.block_info_cache[gene_name] = {
            'num_blocks': num_blocks,
            'blocks': blocks,
            'use_gam': False,
            'node_sequence': ''
        }
        
        return self.block_info_cache[gene_name]

    def process_sample(self, sample_id, rgi_file, pipeline_file, pipeline_only_dir, bam_file, gam_file, output_dir):
        """Xử lý một mẫu hoàn chỉnh"""
        print(f"\n{'='*60}")
        print(f"Xử lý mẫu: {sample_id}")
        print(f"{'='*60}")
        
        self.block_info_cache = {}
        
        sample_output_dir = Path(output_dir) / sample_id
        self.clear_sample_output(sample_output_dir)
        
        print("\n[1] SO SÁNH KẾT QUẢ RGI vs PIPELINE")
        print("-" * 40)
        
        comparison = self.compare_results(rgi_file, pipeline_file)
        
        print(f"\n📊 Kết quả so sánh:")
        print(f"   Số lượng gene RGI (identity > {self.RGI_THRESHOLD}%): {comparison['rgi_count']}")
        print(f"   Số lượng gene pipeline: {comparison['pipeline_count']}")
        print(f"   Số lượng gene chung: {comparison['common_count']}")
        print(f"   Số lượng gene chỉ pipeline: {len(comparison['pipeline_only'])}")
        print(f"   Số lượng gene chỉ RGI (pipeline KHÔNG phát hiện): {len(comparison['rgi_only'])}")
        
        # Lưu common genes
        common_gene_list = []
        for gene_norm in comparison['common_genes']:
            rgi_info = next((g for g in comparison['rgi_genes'] if g['normalized'] == gene_norm), {})
            pipeline_info = next((g for g in comparison['pipeline_genes'] if g['normalized'] == gene_norm), {})
            
            gene_name = rgi_info.get('gene_name', '')
            warning = self.aac_warnings.get(gene_name, '')
            
            common_gene_list.append({
                'gene_name': gene_name,
                'pipeline_name': pipeline_info.get('gene_name', ''),
                'rgi_identity': rgi_info.get('identity', 0),
                'normalized': gene_norm,
                'warning': warning
            })
        
        if common_gene_list:
            df_common = pd.DataFrame(common_gene_list)
            df_common.to_csv(sample_output_dir / f'{sample_id}_common_genes.csv', index=False)
            print(f"\n📁 Đã lưu danh sách common genes: {sample_output_dir}/{sample_id}_common_genes.csv")
        
        # Lưu RGI-only genes (chưa có BLAST)
        rgi_only_list = []
        if comparison['rgi_only']:
            for gene_norm in comparison['rgi_only']:
                rgi_info = next((g for g in comparison['rgi_genes'] if g['normalized'] == gene_norm), {})
                gene_name = rgi_info.get('gene_name', '')
                warning = self.aac_warnings.get(gene_name, '')
                
                rgi_only_list.append({
                    'gene_name': gene_name,
                    'normalized': gene_norm,
                    'rgi_identity': rgi_info.get('identity', 0),
                    'warning': warning
                })
            df_rgi_only = pd.DataFrame(rgi_only_list)
            df_rgi_only.to_csv(sample_output_dir / f'{sample_id}_rgi_only_genes.csv', index=False)
            print(f"\n📁 Đã lưu danh sách RGI-only genes: {sample_output_dir}/{sample_id}_rgi_only_genes.csv")
            
            print(f"\n⚠️ DANH SÁCH GENE PIPELINE KHÔNG PHÁT HIỆN (RGI-ONLY):")
            for item in rgi_only_list:
                warning_str = f" ⚠️ {item['warning']}" if item['warning'] else ""
                print(f"   - {item['gene_name']} (identity: {item['rgi_identity']:.2f}%){warning_str}")
        
        # Lưu pipeline-only genes
        pipeline_only_list = []
        for gene_norm in comparison['pipeline_only']:
            pipeline_info = next((g for g in comparison['pipeline_genes'] if g['normalized'] == gene_norm), {})
            pipeline_only_list.append({
                'gene_name': pipeline_info.get('gene_name', ''),
                'normalized': gene_norm,
                'match_percentage': pipeline_info.get('match_percentage', ''),
                'sample_node_sequence': pipeline_info.get('sample_node_sequence', '')
            })
        
        if pipeline_only_list:
            df_pipeline_only = pd.DataFrame(pipeline_only_list)
            df_pipeline_only.to_csv(sample_output_dir / f'{sample_id}_pipeline_only_genes.csv', index=False)
            print(f"\n📁 Đã lưu danh sách pipeline-only genes: {sample_output_dir}/{sample_id}_pipeline_only_genes.csv")
        
        # ===== XỬ LÝ PIPELINE-ONLY GENES (CÓ BLAST + DEPTH) =====
        print(f"\n[2] GENE CHỈ PIPELINE TÌM ĐƯỢC - LẤY THÔNG TIN BLAST & DEPTH")
        print("-" * 40)
        
        pipeline_only_genes = []
        
        # Lấy danh sách normalized names của pipeline-only genes
        pipeline_only_norm_set = set(comparison['pipeline_only'])
        
        blast_folders = self.find_blast_folders(pipeline_only_dir)
        print(f"Tìm thấy {len(blast_folders)} thư mục *_blast")
        
        # Debug
        if self.debug_mode:
            print("\n🔍 PIPELINE-ONLY GENES từ comparison:")
            for gene_norm in sorted(comparison['pipeline_only']):
                pipeline_info = next((g for g in comparison['pipeline_genes'] if g['normalized'] == gene_norm), {})
                print(f"   {gene_norm} → {pipeline_info.get('gene_name', '')}")
            
            print("\n🔍 BLAST FOLDERS tìm thấy:")
            for bf in blast_folders:
                print(f"   {bf['gene_name']}")
        
        for blast_folder in blast_folders:
            gene_name = blast_folder['gene_name']
            folder_path = blast_folder['folder_path']
            
            # Sử dụng normalized name để kiểm tra
            gene_norm = self.normalize_gene_name(gene_name, source='pipeline')
            is_pipeline_only = gene_norm in pipeline_only_norm_set
            
            if self.debug_mode:
                print(f"\n🔍 Kiểm tra {gene_name}: normalized='{gene_norm}', is_pipeline_only={is_pipeline_only}")
            
            if not is_pipeline_only:
                print(f"   ⏭️ Bỏ qua {gene_name} (không phải pipeline-only, normalized='{gene_norm}')")
                continue
            
            report_files = list(Path(folder_path).glob('*_report.txt'))
            
            if report_files:
                blast_file = report_files[0]
                print(f"\n🔬 Pipeline-only Gene: {gene_name}")
                print(f"   📁 {blast_folder['folder_name']}")
                
                blast_info = self.parse_blast_report(str(blast_file))
                
                print(f"   📊 BLAST results:")
                print(f"      Identity: {blast_info.get('identity', 'N/A')}%")
                print(f"      Coverage: {blast_info.get('coverage', 'N/A')}%")
                print(f"      E-value: {blast_info.get('evalue', 'N/A')}")
                
                # Tìm pipeline info để lấy node sequence
                pipeline_info = next((g for g in pipeline_only_list if g['gene_name'] == gene_name), {})
                
                block_info = self.get_block_info_from_user(
                    gene_name, 
                    blast_info, 
                    pipeline_info,
                    pipeline_only_list,
                    pipeline_file
                )
                num_blocks = block_info['num_blocks']
                use_gam = block_info.get('use_gam', False)
                
                depth_info = {
                    'num_blocks': num_blocks,
                    'block_depths': [],
                    'average_depth': 0.0,
                    'total_depth': 0.0,
                    'use_gam': use_gam
                }
                
                if use_gam:
                    # Tính depth từ GAM
                    print(f"\n   📊 Gene không có trong reference, tính depth từ GAM...")
                    
                    node_sequence = block_info.get('node_sequence', '')
                    if gam_file and os.path.exists(gam_file) and node_sequence:
                        depth = self.get_average_depth_from_gam(gam_file, node_sequence, gene_name)
                        depth_info['average_depth'] = depth
                        depth_info['total_depth'] = depth
                        depth_info['block_depths'] = [depth]
                        print(f"\n      ✅ Depth trung bình từ GAM: {depth:.2f} reads/node")
                    else:
                        if not gam_file or not os.path.exists(gam_file):
                            print(f"      ⚠️ Không tìm thấy GAM file: {gam_file}")
                        if not node_sequence:
                            print(f"      ⚠️ Không có node sequence cho {gene_name}")
                        depth_info['average_depth'] = 0.0
                        depth_info['total_depth'] = 0.0
                        depth_info['block_depths'] = [0.0]
                    
                elif bam_file and os.path.exists(bam_file):
                    # Tính depth từ BAM (gene có trong reference)
                    print(f"\n   📊 Lấy depth cho {num_blocks} block(s) từ BAM:")
                    print(f"   BAM file: {bam_file}")
                    
                    blocks = block_info['blocks']
                    for i, block in enumerate(blocks, 1):
                        contig = block.get('contig', '')
                        start = block.get('start')
                        end = block.get('end')
                        
                        if contig and start and end:
                            d = self.get_average_depth(bam_file, contig, start, end)
                            depth_info['block_depths'].append(d)
                            depth_info['total_depth'] += d
                            print(f"      Block {i}: {contig} {start}-{end} → depth={d:.2f}")
                        else:
                            print(f"      Block {i}: ⚠️ Không có thông tin đầy đủ")
                            depth_info['block_depths'].append(0.0)
                    
                    depth_info['average_depth'] = depth_info['total_depth'] / num_blocks if num_blocks > 0 else 0.0
                    depth = depth_info['average_depth']
                    
                    print(f"\n   📊 Tổng hợp depth từ BAM:")
                    for i, d in enumerate(depth_info['block_depths'], 1):
                        print(f"      Block {i}: {d:.2f}")
                    print(f"      Average depth: {depth_info['average_depth']:.2f}")
                    print(f"      Total depth: {depth_info['total_depth']:.2f}")
                else:
                    print(f"   ⚠️ Không có BAM hoặc GAM file để lấy depth")
                    depth_info['average_depth'] = 0.0
                    depth_info['total_depth'] = 0.0
                    depth_info['block_depths'] = [0.0]
                
                pipeline_only_genes.append({
                    'gene_name': gene_name,
                    'folder_name': blast_folder['folder_name'],
                    'identity': blast_info.get('identity'),
                    'coverage': blast_info.get('coverage'),
                    'evalue': blast_info.get('evalue'),
                    'contig': blast_info.get('contig'),
                    'num_blocks': num_blocks,
                    'block_depths': depth_info['block_depths'],
                    'depth': depth_info['average_depth'],
                    'total_depth': depth_info['total_depth'],
                    'found': blast_info.get('found', False),
                    'use_gam': use_gam,
                    'node_count': len(block_info.get('node_sequence', '').split(',')) if block_info.get('node_sequence') else 0
                })
                
                if blast_info.get('found'):
                    print(f"   ✅ Phát hiện thành công")
                else:
                    print(f"   ⚠️ Không tìm thấy trong BLAST report")
            else:
                print(f"\n⚠️ Không tìm thấy file report trong: {folder_path}")
        
        if pipeline_only_genes:
            df_only = pd.DataFrame(pipeline_only_genes)
            df_only.to_csv(sample_output_dir / f'{sample_id}_pipeline_only_genes_with_blast.csv', index=False)
        
        # ===== XỬ LÝ RGI-ONLY GENES (CHỈ BLAST, KHÔNG DEPTH) =====
        print(f"\n[3] GENE CHỈ RGI PHÁT HIỆN - LẤY THÔNG TIN BLAST (KHÔNG DEPTH)")
        print("-" * 40)
        
        rgi_only_blast_genes = []
        
        if comparison['rgi_only']:
            # Tìm thư mục blast cho RGI-only genes
            blast_folders = self.find_blast_folders(pipeline_only_dir)
            
            # Lọc các blast folder thuộc RGI-only
            rgi_only_blast_folders = []
            for blast_folder in blast_folders:
                gene_name = blast_folder['gene_name']
                # Sử dụng normalized name để so sánh
                blast_norm = self.normalize_gene_name(gene_name, source='pipeline')
                for rgi_gene in rgi_only_list:
                    if blast_norm == rgi_gene['normalized']:
                        rgi_only_blast_folders.append({
                            'folder_name': blast_folder['folder_name'],
                            'folder_path': blast_folder['folder_path'],
                            'gene_name': gene_name,
                            'rgi_info': rgi_gene
                        })
                        break
            
            print(f"Tìm thấy {len(rgi_only_blast_folders)} thư mục *_blast cho RGI-only genes")
            
            for blast_folder in rgi_only_blast_folders:
                gene_name = blast_folder['gene_name']
                rgi_info = blast_folder['rgi_info']
                folder_path = blast_folder['folder_path']
                
                report_files = list(Path(folder_path).glob('*_report.txt'))
                
                if report_files:
                    blast_file = report_files[0]
                    print(f"\n🔬 RGI-only Gene: {gene_name}")
                    print(f"   📁 {blast_folder['folder_name']}")
                    print(f"   📋 RGI Info: identity={rgi_info['rgi_identity']:.2f}%")
                    
                    blast_info = self.parse_blast_report(str(blast_file))
                    
                    print(f"   📊 BLAST results (chỉ thông tin, không có depth):")
                    print(f"      Identity: {blast_info.get('identity', 'N/A')}%")
                    print(f"      Coverage: {blast_info.get('coverage', 'N/A')}%")
                    print(f"      E-value: {blast_info.get('evalue', 'N/A')}")
                    print(f"      Contig: {blast_info.get('contig', 'N/A')}")
                    
                    rgi_only_blast_genes.append({
                        'gene_name': gene_name,
                        'rgi_identity': rgi_info['rgi_identity'],
                        'folder_name': blast_folder['folder_name'],
                        'identity': blast_info.get('identity'),
                        'coverage': blast_info.get('coverage'),
                        'evalue': blast_info.get('evalue'),
                        'contig': blast_info.get('contig'),
                        'found': blast_info.get('found', False),
                        'warning': rgi_info.get('warning', '')
                    })
                    
                    if blast_info.get('found'):
                        print(f"   ✅ Phát hiện thành công")
                    else:
                        print(f"   ⚠️ Không tìm thấy trong BLAST report")
                else:
                    print(f"\n⚠️ Không tìm thấy file report cho gene {gene_name} trong: {folder_path}")
            
            # Lưu RGI-only genes với BLAST info
            if rgi_only_blast_genes:
                df_rgi_only_blast = pd.DataFrame(rgi_only_blast_genes)
                df_rgi_only_blast.to_csv(sample_output_dir / f'{sample_id}_rgi_only_genes_with_blast.csv', index=False)
                print(f"\n📁 Đã lưu danh sách RGI-only genes với BLAST info: {sample_output_dir}/{sample_id}_rgi_only_genes_with_blast.csv")
            
            # Tìm RGI-only genes không có BLAST
            rgi_only_no_blast = []
            for rgi_gene in rgi_only_list:
                found = False
                for blast_gene in rgi_only_blast_genes:
                    if blast_gene['gene_name'] == rgi_gene['gene_name']:
                        found = True
                        break
                if not found:
                    rgi_only_no_blast.append(rgi_gene)
            
            if rgi_only_no_blast:
                df_rgi_only_no_blast = pd.DataFrame(rgi_only_no_blast)
                df_rgi_only_no_blast.to_csv(sample_output_dir / f'{sample_id}_rgi_only_genes_no_blast.csv', index=False)
                print(f"\n📁 Đã lưu danh sách RGI-only genes KHÔNG có BLAST: {sample_output_dir}/{sample_id}_rgi_only_genes_no_blast.csv")
        else:
            print("   ✅ Không có RGI-only genes")
        
        # ===== XỬ LÝ COMMON GENES (CÓ BLAST) =====
        print(f"\n[4] GENE CHUNG (COMMON) - LẤY THÔNG TIN BLAST (KHÔNG DEPTH)")
        print("-" * 40)
        
        common_blast_genes = []
        
        if comparison['common_genes']:
            blast_folders = self.find_blast_folders(pipeline_only_dir)
            
            # Lọc các blast folder thuộc common genes
            common_blast_folders = []
            for blast_folder in blast_folders:
                gene_name = blast_folder['gene_name']
                blast_norm = self.normalize_gene_name(gene_name, source='pipeline')
                if blast_norm in comparison['common_genes']:
                    common_blast_folders.append({
                        'folder_name': blast_folder['folder_name'],
                        'folder_path': blast_folder['folder_path'],
                        'gene_name': gene_name,
                        'normalized': blast_norm
                    })
            
            print(f"Tìm thấy {len(common_blast_folders)} thư mục *_blast cho common genes")
            
            for blast_folder in common_blast_folders:
                gene_name = blast_folder['gene_name']
                folder_path = blast_folder['folder_path']
                
                report_files = list(Path(folder_path).glob('*_report.txt'))
                
                if report_files:
                    blast_file = report_files[0]
                    print(f"\n🔬 Common Gene với BLAST: {gene_name}")
                    print(f"   📁 {blast_folder['folder_name']}")
                    
                    blast_info = self.parse_blast_report(str(blast_file))
                    
                    print(f"   📊 BLAST results:")
                    print(f"      Identity: {blast_info.get('identity', 'N/A')}%")
                    print(f"      Coverage: {blast_info.get('coverage', 'N/A')}%")
                    print(f"      E-value: {blast_info.get('evalue', 'N/A')}")
                    print(f"      Contig: {blast_info.get('contig', 'N/A')}")
                    
                    common_blast_genes.append({
                        'gene_name': gene_name,
                        'folder_name': blast_folder['folder_name'],
                        'identity': blast_info.get('identity'),
                        'coverage': blast_info.get('coverage'),
                        'evalue': blast_info.get('evalue'),
                        'contig': blast_info.get('contig'),
                        'found': blast_info.get('found', False),
                        'num_blocks': blast_info.get('num_blocks', 1)
                    })
                    
                    if blast_info.get('found'):
                        print(f"   ✅ Phát hiện thành công")
                    else:
                        print(f"   ⚠️ Không tìm thấy trong BLAST report")
                else:
                    print(f"\n⚠️ Không tìm thấy file report cho gene {gene_name} trong: {folder_path}")
            
            # Lưu common genes với BLAST info
            if common_blast_genes:
                df_common_blast = pd.DataFrame(common_blast_genes)
                df_common_blast.to_csv(sample_output_dir / f'{sample_id}_common_genes_with_blast.csv', index=False)
                print(f"\n📁 Đã lưu danh sách common genes với BLAST info: {sample_output_dir}/{sample_id}_common_genes_with_blast.csv")
        else:
            print("   ✅ Không có common genes")
        
        # Lưu kết quả tổng hợp
        result = {
            'sample_id': sample_id,
            'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'rgi_threshold': self.RGI_THRESHOLD,
            'statistics': {
                'rgi_count': comparison['rgi_count'],
                'pipeline_count': comparison['pipeline_count'],
                'common_count': comparison['common_count'],
                'pipeline_only_count': len(comparison['pipeline_only']),
                'rgi_only_count': len(comparison['rgi_only']),
                'rgi_only_with_blast_count': len(rgi_only_blast_genes),
                'rgi_only_no_blast_count': len(rgi_only_no_blast) if 'rgi_only_no_blast' in locals() else 0,
                'common_with_blast_count': len(common_blast_genes)
            },
            'common_genes': common_gene_list,
            'common_genes_with_blast': common_blast_genes,
            'pipeline_only_genes': pipeline_only_genes,
            'rgi_only_genes': rgi_only_list,
            'rgi_only_genes_with_blast': rgi_only_blast_genes,
            'rgi_only_genes_no_blast': rgi_only_no_blast if 'rgi_only_no_blast' in locals() else [],
            'aac_warnings': list(self.aac_warnings.values())
        }
        
        self.results[sample_id] = result
        
        json_file = sample_output_dir / f'{sample_id}_validation_results.json'
        with open(json_file, 'w') as f:
            json.dump(result, f, indent=2, default=str)
        
        print(f"\n✅ Đã lưu kết quả cho {sample_id} vào: {sample_output_dir}")
        print(f"   - {sample_id}_validation_results.json")
        print(f"   - {sample_id}_common_genes.csv")
        if common_blast_genes:
            print(f"   - {sample_id}_common_genes_with_blast.csv")
        if pipeline_only_genes:
            print(f"   - {sample_id}_pipeline_only_genes_with_blast.csv")
        if comparison['rgi_only']:
            print(f"   - {sample_id}_rgi_only_genes.csv")
            if rgi_only_blast_genes:
                print(f"   - {sample_id}_rgi_only_genes_with_blast.csv")
            if rgi_only_no_blast:
                print(f"   - {sample_id}_rgi_only_genes_no_blast.csv")
        
        return result

def main():
    """Hàm chính"""
    BASE_DIR = "/Users/shinra/gene_detection_v2/validate_database"
    OUTPUT_DIR = f"{BASE_DIR}/validation_results"
    
    samples = [
        {
            'id': 'SAMN48847114',
            'rgi_file': f'{BASE_DIR}/SAMN48847114/rgi_result/rgi_strict.txt',
            'pipeline_file': f'{BASE_DIR}/SAMN48847114/filtered_detected_AMR_gene/SAMN48847114_vs_amr_database_filtered_v3.tsv',
            'pipeline_only_dir': f'{BASE_DIR}/SAMN48847114',
            'bam_file': f'{BASE_DIR}/SAMN48847114/reads_mapped_to_graph/SAMN48847114.all.bam',
            'gam_file': f'{BASE_DIR}/SAMN48847114/reads_mapped_to_graph/SAMN48847114.gam'
        },
    ]
    
    validator = AMRGeneValidator(BASE_DIR)
    
    if len(sys.argv) > 1:
        try:
            validator.RGI_THRESHOLD = float(sys.argv[1])
            print(f"✅ Sử dụng threshold identity = {validator.RGI_THRESHOLD}%")
        except:
            pass
    
    for sample in samples:
        if not os.path.exists(sample['rgi_file']):
            print(f"⚠️ Không tìm thấy file RGI: {sample['rgi_file']}")
            continue
            
        if not os.path.exists(sample['pipeline_file']):
            print(f"⚠️ Không tìm thấy file pipeline: {sample['pipeline_file']}")
            continue
        
        if not os.path.exists(sample['bam_file']):
            print(f"⚠️ Không tìm thấy BAM file: {sample['bam_file']}")
            print(f"   Vẫn tiếp tục nhưng sẽ không lấy được depth từ BAM")
        
        if not os.path.exists(sample['gam_file']):
            print(f"⚠️ Không tìm thấy GAM file: {sample['gam_file']}")
            print(f"   Vẫn tiếp tục nhưng sẽ không lấy được depth từ GAM")
        
        validator.process_sample(
            sample_id=sample['id'],
            rgi_file=sample['rgi_file'],
            pipeline_file=sample['pipeline_file'],
            pipeline_only_dir=sample.get('pipeline_only_dir'),
            bam_file=sample.get('bam_file'),
            gam_file=sample.get('gam_file'),
            output_dir=OUTPUT_DIR
        )
    
    print("\n" + "="*60)
    print("✅ HOÀN THÀNH XỬ LÝ TẤT CẢ MẪU")
    print("="*60)
    print(f"\n📁 Kết quả lưu tại: {OUTPUT_DIR}")
    print("   Mỗi sample có thư mục riêng với kết quả chi tiết")

if __name__ == "__main__":
    main()
