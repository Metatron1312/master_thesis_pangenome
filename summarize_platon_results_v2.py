import os
import json
import csv
from collections import defaultdict

def process_json_files(root_dir, output_file):
    plasmid_data = defaultdict(list)
    all_contigs = []
    
    for biosample_id in os.listdir(root_dir):
        biosample_path = os.path.join(root_dir, biosample_id)
        
        if not os.path.isdir(biosample_path):
            continue
            
        for filename in os.listdir(biosample_path):
            if not filename.endswith('.json'):
                continue
                
            json_path = os.path.join(biosample_path, filename)
            
            try:
                with open(json_path, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                    
                for contig_id, contig_data in data.items():
                    # Xử lý plasmid hits
                    plasmid_id = "No plasmid hit"
                    if 'plasmid_hits' in contig_data and contig_data['plasmid_hits']:
                        plasmid_id = contig_data['plasmid_hits'][0]['plasmid']['id']
                    
                    # Xử lý circular
                    is_circular = False
                    if 'is_circular' in contig_data:
                        circ_value = contig_data['is_circular']
                        if isinstance(circ_value, bool):
                            is_circular = circ_value
                        elif isinstance(circ_value, str):
                            is_circular = circ_value.lower() == 'true'
                    
                    # Xử lý replication
                    has_replicon = 'replication_hits' in contig_data and bool(contig_data['replication_hits'])
                    
                    # Xử lý inc types - CHỈNH SỬA QUAN TRỌNG Ở ĐÂY
                    inc_types = []
                    if 'inc_types' in contig_data and contig_data['inc_types']:
                        for inc in contig_data['inc_types']:
                            if 'type' in inc:
                                inc_types.append(inc['type'])
                    inc_type_str = ', '.join(inc_types) if inc_types else 'None'
                    
                    contig_info = {
                        'Plasmid ID': plasmid_id,
                        'Biosample ID': biosample_id,
                        'Contig ID': contig_id,
                        'Length': contig_data.get('length', 'N/A'),
                        'Inc Type': inc_type_str,
                        'ORFs': len(contig_data.get('orfs', {})),
                        'Circular': 'Yes' if is_circular else 'No',
                        'Replicon': 'Yes' if has_replicon else 'No',
                        'OriT': len(contig_data.get('orit_hits', [])),
                        'Mobilization': len(contig_data.get('mobilization_hits', [])),
                        'Conjugation': len(contig_data.get('conjugation_hits', [])),
                        'AMR': len(contig_data.get('amr_hits', [])),
                        'rRNA': len(contig_data.get('rrnas', [])),
                        'Coverage': contig_data.get('coverage', 'N/A'),
                        'Protein Score': contig_data.get('protein_score', 'N/A')
                    }
                    
                    plasmid_data[plasmid_id].append(contig_info)
                    all_contigs.append(contig_info)
                    
            except Exception as e:
                print(f"Lỗi khi xử lý {json_path}: {str(e)}")
                continue
    
    # Ghi kết quả vào CSV
    with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = [
            'Plasmid ID', 'Biosample ID', 'Contig ID', 'Length', 'Inc Type',
            'ORFs', 'Circular', 'Replicon', 'OriT', 'Mobilization',
            'Conjugation', 'AMR', 'rRNA', 'Coverage', 'Protein Score'
        ]
        
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        current_plasmid = None
        for contig in sorted(all_contigs, key=lambda x: x['Plasmid ID']):
            if contig['Plasmid ID'] != current_plasmid:
                # Thêm dòng header nhóm
                group_row = {
                    'Plasmid ID': f"=== PLASMID GROUP: {contig['Plasmid ID']} ===",
                    'Biosample ID': f"Total contigs: {len(plasmid_data[contig['Plasmid ID']])}",
                    'Contig ID': '',
                    'Length': '',
                    'Inc Type': '',
                    'ORFs': '',
                    'Circular': '',
                    'Replicon': '',
                    'OriT': '',
                    'Mobilization': '',
                    'Conjugation': '',
                    'AMR': '',
                    'rRNA': '',
                    'Coverage': '',
                    'Protein Score': ''
                }
                writer.writerow(group_row)
                current_plasmid = contig['Plasmid ID']
            
            writer.writerow(contig)

if __name__ == '__main__':
    input_directory = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/Platon_results"
    output_csv = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/Platon_results/summarized_platon_results.csv"
    
    if not os.path.exists(input_directory):
        print(f"LỖI: Thư mục đầu vào không tồn tại: {input_directory}")
    else:
        process_json_files(input_directory, output_csv)
        print(f"Đã tạo báo cáo thành công: {output_csv}")