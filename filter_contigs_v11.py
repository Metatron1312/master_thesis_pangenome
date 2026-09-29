import os
import csv
from collections import defaultdict
from datetime import datetime
import shutil

# ==================== CẤU HÌNH ĐƯỜNG DẪN ====================
THU_MUC_GOC = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/Platon_results_annotated"
THU_MUC_KRAKEN2 = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/Kraken2_results_v2"
DUONG_DAN_REPORT = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/Kraken2_results_v2/combined_reports_v2.csv"
DUONG_DAN_BLAST = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/Kraken2_results_v2/blast_results_filtered.csv"
# =============================================================

def tao_thu_muc_filter_contigs():
    """Tạo thư mục output"""
    thu_muc_output = os.path.join(os.getcwd(), "filter_contigs")
    cac_thu_muc_con = [
        'filtered',  # Luôn chứa AB + nhiễm (để kiểm tra sau)
        'contaminated_chromosome',  # Chứa nhiễm từ chromosome (để kiểm tra kháng sinh)
        'contaminated_plasmid',  # Chứa nhiễm từ plasmid (nếu ≥1%)
        'reports'
    ]
    
    os.makedirs(thu_muc_output, exist_ok=True)
    for thu_muc in cac_thu_muc_con:
        os.makedirs(os.path.join(thu_muc_output, thu_muc), exist_ok=True)
    
    return thu_muc_output

def phan_tich_file_report(duong_dan_report):
    """Đọc file report và loại bỏ Homo sapiens"""
    thong_tin = []
    taxid_mapping = {}
    human_taxids = set()
    
    try:
        with open(duong_dan_report, 'r', encoding='utf-8') as f:
            delimiter = '\t' if '\t' in f.readline() else ','
            f.seek(0)
            reader = csv.DictReader(f, delimiter=delimiter)
            
            current_main_taxid = None
            for row in reader:
                rank_code = row.get('rank_code', '').strip().upper()
                if rank_code in ['S', 'S1', 'S2']:
                    scientific_name = row['scientific_name'].strip().lower()
                    
                    # Loại bỏ Homo sapiens hoàn toàn
                    if 'homo sapiens' in scientific_name:
                        if rank_code == 'S':
                            human_taxids.add(row['taxonomy_id'])
                        continue
                    
                    entry = {
                        'biosample_id': row['biosample_id'],
                        'file_type': row['file_type'],
                        'clade_fragments': int(row['clade_fragments']),
                        'taxonomy_id': row['taxonomy_id'],
                        'scientific_name': scientific_name,
                        'rank_code': rank_code
                    }
                    thong_tin.append(entry)
                    
                    if rank_code == 'S':
                        current_main_taxid = row['taxonomy_id']
                    elif rank_code in ['S1', 'S2'] and current_main_taxid:
                        taxid_mapping[row['taxonomy_id']] = current_main_taxid
            
            print(f"Đã đọc {len(thong_tin)} mục rank S/S1/S2 (đã loại bỏ Homo sapiens)")
            return thong_tin, taxid_mapping, human_taxids
            
    except Exception as e:
        print(f"Lỗi khi đọc file report: {str(e)}")
        return [], {}, set()

def tinh_toan_so_lieu(thong_tin_report):
    """Tính toán tỷ lệ từ các entry rank S"""
    ab_contigs = 0
    other_contigs = defaultdict(int)
    ab_taxids = set()
    
    for entry in thong_tin_report:
        if entry['rank_code'] == 'S':
            if 'acinetobacter baumannii' in entry['scientific_name']:
                ab_contigs += entry['clade_fragments']
                ab_taxids.add(entry['taxonomy_id'])
            else:
                other_contigs[entry['scientific_name']] += entry['clade_fragments']
    
    return {
        'total_contigs': ab_contigs + sum(other_contigs.values()),
        'ab_contigs': ab_contigs,
        'other_contigs': other_contigs,
        'ab_taxids': ab_taxids
    }

def doc_contig_ids_tu_output(duong_dan_output, taxid_mapping, ab_taxids, human_taxids, thong_tin_report):
    """Đọc file output Kraken2 và loại bỏ Homo sapiens"""
    contig_ids_ab = set()
    contig_ids_khac = {}  # {contig_id: {'taxid', 'main_taxid', 'name'}}
    contig_ids_unclassified = set()
    
    taxid_to_info = {}
    for entry in thong_tin_report:
        if entry['rank_code'] == 'S':
            taxid_to_info[entry['taxonomy_id']] = {
                'name': entry['scientific_name'],
                'main_taxid': entry['taxonomy_id']
            }
        elif entry['taxonomy_id'] in taxid_mapping:
            taxid_to_info[entry['taxonomy_id']] = {
                'name': entry['scientific_name'],
                'main_taxid': taxid_mapping[entry['taxonomy_id']]
            }
    
    try:
        with open(duong_dan_output, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    status = parts[0]  # C/U
                    contig_id = parts[1]
                    taxid = parts[2]
                    
                    if status == 'C':
                        if taxid in taxid_to_info:
                            info = taxid_to_info[taxid]
                            main_taxid = info['main_taxid']
                            
                            # Bỏ qua Homo sapiens
                            if main_taxid in human_taxids:
                                continue
                                
                            if main_taxid in ab_taxids:
                                contig_ids_ab.add(contig_id)
                            else:
                                contig_ids_khac[contig_id] = {
                                    'taxid': taxid,
                                    'main_taxid': main_taxid,
                                    'name': info['name']
                                }
                    elif status == 'U':
                        contig_ids_unclassified.add(contig_id)
                        
        print(f"Phát hiện: {len(contig_ids_ab)} contigs AB, {len(contig_ids_khac)} contigs khác, {len(contig_ids_unclassified)} unclassified")
        return contig_ids_ab, contig_ids_khac, contig_ids_unclassified
    
    except Exception as e:
        print(f"LỖI khi đọc file output: {str(e)}")
        return set(), {}, set()

def trich_xuat_contigs(duong_dan_fasta, contig_ids, duong_dan_luu):
    """Trích xuất contigs từ file fasta"""
    os.makedirs(os.path.dirname(duong_dan_luu), exist_ok=True)
    
    try:
        contigs_found = 0
        with open(duong_dan_fasta, 'r') as f_in, open(duong_dan_luu, 'w') as f_out:
            current_contig = None
            write_contig = False
            
            for line in f_in:
                if line.startswith('>'):
                    current_contig = line.split()[0][1:]
                    write_contig = current_contig in contig_ids
                    if write_contig:
                        contigs_found += 1
                        f_out.write(line)
                elif write_contig:
                    f_out.write(line)
        
        print(f"ĐÃ TRÍCH XUẤT {contigs_found}/{len(contig_ids)} CONTIGS")
        return contigs_found > 0
        
    except Exception as e:
        print(f"LỖI khi trích xuất contigs: {str(e)}")
        return False

def xu_ly_unclassified(contig_ids_unclassified):
    """Xử lý contigs unclassified bằng BLAST"""
    contig_ids_ab = set()
    contig_ids_khac = set()
    
    if DUONG_DAN_BLAST and os.path.exists(DUONG_DAN_BLAST):
        try:
            with open(DUONG_DAN_BLAST, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    contig_id = row.get('Contig_ID', '')
                    if contig_id in contig_ids_unclassified:
                        if 'acinetobacter baumannii' in row.get('Description', '').lower():
                            contig_ids_ab.add(contig_id)
                        else:
                            contig_ids_khac.add(contig_id)
            
            print(f"Phân loại unclassified: {len(contig_ids_ab)} AB, {len(contig_ids_khac)} không phải AB")
            return contig_ids_ab, contig_ids_khac
            
        except Exception as e:
            print(f"Lỗi khi đọc file blast: {str(e)}")
    
    return set(), set()

def tinh_toan_ty_le_tu_contigs(contig_ids_ab, contig_ids_khac, ab_taxids):
    """Tính toán tỷ lệ nhiễm từ contigs"""
    total_contigs = len(contig_ids_ab) + len(contig_ids_khac)
    contaminated_contigs = len(contig_ids_khac)
    ty_le_nhiem = (contaminated_contigs / total_contigs) * 100 if total_contigs > 0 else 0
    
    return {
        'total_contigs': total_contigs,
        'ab_contigs': len(contig_ids_ab),
        'contaminated_contigs': contaminated_contigs,
        'ty_le_nhiem': ty_le_nhiem
    }

def xu_ly_chromosome(thong_tin_mau, taxid_mapping, duong_dan_output, duong_dan_fasta, thu_muc_output, report_chi_tiet, human_taxids):
    """Xử lý chromosome: Loại bỏ nhiễm ≥10% trước, sau đó xử lý phần còn lại như case 1-3"""
    so_lieu = tinh_toan_so_lieu([e for e in thong_tin_mau if e['rank_code'] == 'S'])
    ab_taxids = so_lieu['ab_taxids']
    
    # Bước 1: Đọc contigs từ file Kraken2
    contig_ids_ab, contig_ids_khac, contig_ids_unclassified = doc_contig_ids_tu_output(
        duong_dan_output, taxid_mapping, ab_taxids, human_taxids, thong_tin_mau)
    
    # Bước 2: Xử lý unclassified bằng BLAST
    contig_ids_ab_blast, contig_ids_nhiem_blast = xu_ly_unclassified(contig_ids_unclassified)
    contig_ids_ab.update(contig_ids_ab_blast)
    for cid in contig_ids_nhiem_blast:
        contig_ids_khac[cid] = {'taxid': 'unclassified', 'main_taxid': 'unclassified', 'name': 'unclassified (blast)'}
    
    # Bước 3: Tính tỷ lệ ban đầu
    ty_le_ban_dau = tinh_toan_ty_le_tu_contigs(contig_ids_ab, contig_ids_khac, ab_taxids)
    print(f"Tỷ lệ nhiễm ban đầu: {ty_le_ban_dau['ty_le_nhiem']:.2f}%")

    # Bước 4: Xác định contigs cần xử lý
    contigs_giu_lai = contig_ids_ab.copy()
    contigs_nhiem = set()
    action = ""

    # Case 4: Nếu nhiễm ≥10% - LỌC CÁC TAXON NHIỄM ≥10% TRƯỚC
    if ty_le_ban_dau['ty_le_nhiem'] >= 10:
        # Tính toán main_taxid_counts
        main_taxid_counts = defaultdict(int)
        for cid, info in contig_ids_khac.items():
            main_taxid_counts[info['main_taxid']] += 1
        
        # Xác định các taxid có tỷ lệ ≥10%
        main_taxid_contaminated = set()
        for taxid, count in main_taxid_counts.items():
            if taxid not in ab_taxids:
                taxon_rate = (count / ty_le_ban_dau['total_contigs']) * 100
                if taxon_rate >= 10:
                    main_taxid_contaminated.add(taxid)
        
        # Loại bỏ các contigs thuộc taxon nhiễm ≥10%
        contigs_nhiem_geq10 = {cid for cid, info in contig_ids_khac.items() 
                              if info['main_taxid'] in main_taxid_contaminated}
        
        # Giữ lại các contigs không bị loại bỏ (AB + nhiễm <10%)
        contig_ids_khac = {cid: info for cid, info in contig_ids_khac.items() 
                          if cid not in contigs_nhiem_geq10}
        
        action = f"Filtered taxa with ≥10% contamination (removed {len(contigs_nhiem_geq10)} contigs)"
        
        # Tính toán lại tỷ lệ sau khi loại bỏ
        ty_le = tinh_toan_ty_le_tu_contigs(contig_ids_ab, contig_ids_khac, ab_taxids)
        print(f"Tỷ lệ nhiễm sau khi lọc ≥10%: {ty_le['ty_le_nhiem']:.2f}%")
    else:
        ty_le = ty_le_ban_dau

    # Bước 5: Xử lý phần còn lại theo case 1-3
    ten_file_base = os.path.basename(duong_dan_fasta).replace('.fasta', '')
    duong_dan_filtered = os.path.join(thu_muc_output, 'filtered', f"{ten_file_base}.filtered.fasta")
    
    # LUÔN giữ lại AB + nhiễm còn lại trong file filtered
    contigs_giu_lai.update(contig_ids_khac.keys())
    
    # Tạo file filtered
    if os.path.exists(duong_dan_filtered):
        os.remove(duong_dan_filtered)
    
    if not trich_xuat_contigs(duong_dan_fasta, contigs_giu_lai, duong_dan_filtered):
        print("KHÔNG THỂ TẠO FILE FILTERED")
        return

    # Xác định contigs nhiễm cần copy vào contaminated (tùy case)
    if ty_le['ty_le_nhiem'] < 1:
        action += "; Keep all remaining contigs (contamination <1%)"
        contigs_nhiem = set()  # Không lưu nhiễm
    
    elif 1 <= ty_le['ty_le_nhiem'] < 5:
        action += "; Filtered non-Acinetobacter from remaining contigs (1-5% contamination)"
        contigs_nhiem = {cid for cid, info in contig_ids_khac.items() 
                        if 'acinetobacter' not in info['name']}
    
    elif ty_le['ty_le_nhiem'] >= 5:
        action += "; Filtered all remaining contaminants (≥5% contamination)"
        contigs_nhiem = set(contig_ids_khac.keys())  # Tất cả nhiễm còn lại

    # Lưu contigs nhiễm vào contaminated (nếu có)
    if contigs_nhiem:
        ten_file_nhiem = f"{ten_file_base}.contaminated.fasta"
        duong_dan_nhiem = os.path.join(thu_muc_output, 'contaminated_chromosome', ten_file_nhiem)
        trich_xuat_contigs(duong_dan_fasta, contigs_nhiem, duong_dan_nhiem)

    # Báo cáo
    report_chi_tiet.append({
        'Biosample': os.path.basename(duong_dan_fasta).split('.')[0],
        'File_Type': 'chromosome',
        'Total_Contigs': ty_le_ban_dau['total_contigs'],
        'AB_Contigs': ty_le['ab_contigs'],
        'Contaminated_Contigs': ty_le['contaminated_contigs'],
        'Contamination_Rate': f"{ty_le['ty_le_nhiem']:.2f}%",
        'Action': action,
        'Output_File': duong_dan_filtered
    })

def xu_ly_plasmid(thong_tin_mau, taxid_mapping, duong_dan_output, duong_dan_fasta, thu_muc_output, report_chi_tiet, human_taxids):
    """Xử lý plasmid: Giữ AB + nhiễm trong filtered, copy nhiễm (≥1%) vào contaminated"""
    so_lieu = tinh_toan_so_lieu([e for e in thong_tin_mau if e['rank_code'] == 'S'])
    ab_taxids = so_lieu['ab_taxids']
    
    contig_ids_ab, contig_ids_khac, contig_ids_unclassified = doc_contig_ids_tu_output(
        duong_dan_output, taxid_mapping, ab_taxids, human_taxids, thong_tin_mau)
    
    contig_ids_ab_blast, contig_ids_nhiem_blast = xu_ly_unclassified(contig_ids_unclassified)
    contig_ids_ab.update(contig_ids_ab_blast)
    for cid in contig_ids_nhiem_blast:
        contig_ids_khac[cid] = {'taxid': 'unclassified', 'main_taxid': 'unclassified', 'name': 'unclassified (blast)'}
    
    ty_le = tinh_toan_ty_le_tu_contigs(contig_ids_ab, contig_ids_khac, ab_taxids)
    
    ten_file_base = os.path.basename(duong_dan_fasta).replace('.fasta', '')
    duong_dan_filtered = os.path.join(thu_muc_output, 'filtered', f"{ten_file_base}.filtered.fasta")
    
    # LUÔN giữ lại AB + nhiễm trong file filtered
    contigs_giu_lai = contig_ids_ab.union(contig_ids_khac.keys())
    
    if os.path.exists(duong_dan_filtered):
        os.remove(duong_dan_filtered)
    
    if not trich_xuat_contigs(duong_dan_fasta, contigs_giu_lai, duong_dan_filtered):
        print("KHÔNG THỂ TẠO FILE FILTERED PLASMID")
        return
    
    # Chỉ lưu nhiễm vào contaminated nếu tỷ lệ ≥1%
    if ty_le['ty_le_nhiem'] >= 1:
        ten_file_nhiem = f"{ten_file_base}.contaminated.fasta"
        duong_dan_nhiem = os.path.join(thu_muc_output, 'contaminated_plasmid', ten_file_nhiem)
        trich_xuat_contigs(duong_dan_fasta, contig_ids_khac.keys(), duong_dan_nhiem)
    
    report_chi_tiet.append({
        'Biosample': os.path.basename(duong_dan_fasta).split('.')[0],
        'File_Type': 'plasmid',
        'Total_Contigs': ty_le['total_contigs'],
        'AB_Contigs': ty_le['ab_contigs'],
        'Contaminated_Contigs': ty_le['contaminated_contigs'],
        'Contamination_Rate': f"{ty_le['ty_le_nhiem']:.2f}%",
        'Action': 'Filtered non-Acinetobacter (plasmid)',
        'Output_File': duong_dan_filtered
    })

def luu_report(report_data, thu_muc_output, ten_file):
    """Lưu report"""
    duong_dan_report = os.path.join(thu_muc_output, 'reports', ten_file)
    
    try:
        with open(duong_dan_report, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=report_data[0].keys())
            writer.writeheader()
            writer.writerows(report_data)
        print(f"Đã lưu report: {duong_dan_report}")
    except Exception as e:
        print(f"Lỗi khi lưu report: {str(e)}")

def main():
    """Hàm chính"""
    print("\n=== BẮT ĐẦU XỬ LÝ ===")
    print(f"Thời gian: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    thu_muc_output = tao_thu_muc_filter_contigs()
    report_tong_hop = []
    report_chi_tiet = []
    
    thong_tin_report, taxid_mapping, human_taxids = phan_tich_file_report(DUONG_DAN_REPORT)
    
    if not thong_tin_report:
        print("Không có dữ liệu để xử lý")
        return
    
    # Gom nhóm theo biosample và loại file
    report_theo_mau = defaultdict(list)
    for entry in thong_tin_report:
        key = (entry['biosample_id'], entry['file_type'])
        report_theo_mau[key].append(entry)
    
    # Xử lý từng file Kraken2 output
    for ten_file in os.listdir(THU_MUC_KRAKEN2):
        if not ten_file.endswith('.output'):
            continue
            
        duong_dan_output = os.path.join(THU_MUC_KRAKEN2, ten_file)
        parts = ten_file.split('_')
        biosample_id = parts[0]
        loai_file = 'chromosome' if 'chromosome' in ten_file.lower() else 'plasmid'
        
        thong_tin_mau = report_theo_mau.get((biosample_id, loai_file), [])
        if not thong_tin_mau:
            continue
            
        ten_fasta = f"{biosample_id}.{loai_file}.annotated.fasta"
        duong_dan_fasta = os.path.join(THU_MUC_GOC, biosample_id, ten_fasta)
        
        if not os.path.exists(duong_dan_fasta):
            continue
            
        if loai_file == 'chromosome':
            xu_ly_chromosome(thong_tin_mau, taxid_mapping, duong_dan_output, duong_dan_fasta, thu_muc_output, report_chi_tiet, human_taxids)
        else:
            xu_ly_plasmid(thong_tin_mau, taxid_mapping, duong_dan_output, duong_dan_fasta, thu_muc_output, report_chi_tiet, human_taxids)
    
    # Lưu report
    if report_chi_tiet:
        luu_report(report_chi_tiet, thu_muc_output, 'detailed_report.csv')
        
        # Tạo summary report
        summary = []
        for entry in report_chi_tiet:
            summary.append({
                'Biosample': entry['Biosample'],
                'File_Type': entry['File_Type'],
                'Total_Contigs': entry['Total_Contigs'],
                'AB_Contigs': entry['AB_Contigs'],
                'Contamination_Rate': entry['Contamination_Rate'],
                'Action': entry['Action']
            })
        
        luu_report(summary, thu_muc_output, 'summary_report.csv')
    
    print("\n=== HOÀN TẤT ===")

if __name__ == "__main__":
    main()