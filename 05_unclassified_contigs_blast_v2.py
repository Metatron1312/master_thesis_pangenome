import os
import glob
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

def parse_kraken_output(kraken_folder):
    """
    Lọc các contig có rank_code là U (unclassified) và trích xuất biosample_id
    Đồng thời xác định loại file (chromosome/plasmid) chứa contig
    """
    unclassified_contigs = {}
    
    for file in glob.glob(os.path.join(kraken_folder, "*.output")):
        filename = os.path.basename(file)
        if "_chromosome.output" in filename:
            seq_type = "chromosome"
            biosample_id = filename.split("_chromosome.output")[0]
        elif "_plasmid.output" in filename:
            seq_type = "plasmid"
            biosample_id = filename.split("_plasmid.output")[0]
        else:
            continue
        
        with open(file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3 and parts[0] == 'U':
                    contig_id = parts[1]
                    key = (biosample_id, seq_type)
                    if key not in unclassified_contigs:
                        unclassified_contigs[key] = []
                    unclassified_contigs[key].append(contig_id)
    
    return unclassified_contigs

def find_fasta_file(target_folder, biosample_id, seq_type):
    """
    Tìm file fasta tương ứng với biosample_id và loại sequence (chromosome/plasmid)
    """
    sample_dir = os.path.join(target_folder, biosample_id)
    if not os.path.exists(sample_dir):
        return None
    
    pattern = f"*{biosample_id}*.{seq_type}.annotated.fasta"
    fasta_files = glob.glob(os.path.join(sample_dir, pattern))
    return fasta_files[0] if fasta_files else None

def extract_contig_sequence(fasta_file, contig_id):
    """Trích xuất trình tự contig từ file fasta"""
    if not fasta_file or not os.path.exists(fasta_file):
        return None
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == contig_id:
            return record
    return None

def run_blast(sequence, query_id):
    """Thực hiện BLAST với trình tự contig"""
    print(f"Running BLAST for {query_id}...")
    try:
        result_handle = NCBIWWW.qblast(
            program="blastn",
            database="nt",
            sequence=sequence.seq,
            hitlist_size=10,  # Lấy nhiều kết quả để lọc sau
            entrez_query='none',
            format_type="XML"
        )
        blast_records = NCBIXML.parse(result_handle)
        return list(blast_records)
    except Exception as e:
        print(f"BLAST failed for {query_id}: {str(e)}")
        return None

def filter_blast_hits(hits, query_length):
    """Lọc kết quả BLAST theo tiêu chí chất lượng"""
    filtered = []
    
    for hit in hits:
        # Tính toán các thông số bổ sung
        coverage = hit['Alignment_Length'] / query_length
        hit['Query_Length'] = query_length
        hit['Coverage'] = round(coverage * 100, 2)
        
        # Tiêu chuẩn lọc chất lượng
        if (hit['E_value'] < 0.001 and              # E-value ngặt nghèo
            hit['Bit_Score'] > 50 and               # Bit Score đủ cao
            hit['Percentage_Identity'] > 70 and     # Độ tương đồng >70%
            coverage > 0.7):                        # Độ phủ >70%
            
            filtered.append(hit)
    
    # Sắp xếp theo Bit Score giảm dần
    filtered.sort(key=lambda x: x['Bit_Score'], reverse=True)
    
    # Chỉ giữ lại 3 kết quả tốt nhất sau khi lọc
    return filtered[:3]

def process_contig(biosample_id, seq_type, contig_id, target_folder):
    """Xử lý một contig cụ thể"""
    fasta_file = find_fasta_file(target_folder, biosample_id, seq_type)
    if not fasta_file:
        print(f"No {seq_type} fasta file found for {biosample_id}")
        return None
    
    contig_seq = extract_contig_sequence(fasta_file, contig_id)
    if not contig_seq:
        print(f"Contig {contig_id} not found in {seq_type} fasta for {biosample_id}")
        return None
    
    query_id = f"{biosample_id}_{seq_type}_{contig_id}"
    blast_results = run_blast(contig_seq, query_id)
    
    if not blast_results:
        return None
    
    # Xử lý và lọc kết quả BLAST
    hits = []
    for blast_record in blast_results:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                hits.append({
                    'Query_ID': query_id,
                    'Sample_ID': biosample_id,
                    'Sequence_Type': seq_type,
                    'Contig_ID': contig_id,
                    'Hit_ID': alignment.hit_id,
                    'Accession': alignment.accession,
                    'Description': alignment.hit_def,
                    'E_value': hsp.expect,
                    'Bit_Score': hsp.bits,
                    'Identity': hsp.identities,
                    'Alignment_Length': hsp.align_length,
                    'Query_Start': hsp.query_start,
                    'Query_End': hsp.query_end,
                    'Subject_Start': hsp.sbjct_start,
                    'Subject_End': hsp.sbjct_end,
                    'Percentage_Identity': round((hsp.identities / hsp.align_length) * 100, 2),
                    'Alignment': str(hsp.sbjct)
                })
    
    # Áp dụng bộ lọc chất lượng
    return filter_blast_hits(hits, len(contig_seq))

def save_to_csv(results, output_file):
    """Lưu kết quả vào file CSV"""
    if not results:
        print("No results to save.")
        return
    
    df = pd.DataFrame(results)
    
    # Sắp xếp các cột để dễ đọc
    columns_order = [
        'Query_ID', 'Sample_ID', 'Sequence_Type', 'Contig_ID', 'Query_Length',
        'Hit_ID', 'Accession', 'Description',
        'E_value', 'Bit_Score', 'Percentage_Identity', 'Identity', 'Alignment_Length', 'Coverage',
        'Query_Start', 'Query_End', 'Subject_Start', 'Subject_End',
        'Alignment'
    ]
    
    df = df[columns_order]
    df.to_csv(output_file, index=False)
    print(f"Results successfully saved to {output_file}")

def main(kraken_folder, target_folder, output_file):
    """Hàm chính thực hiện toàn bộ quy trình"""
    unclassified_contigs = parse_kraken_output(kraken_folder)
    
    if not unclassified_contigs:
        print("No unclassified contigs found in Kraken2 output files.")
        return
    
    all_results = []
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = []
        for (biosample_id, seq_type), contigs in unclassified_contigs.items():
            for contig_id in contigs:
                futures.append(
                    executor.submit(
                        process_contig,
                        biosample_id,
                        seq_type,
                        contig_id,
                        target_folder
                    )
                )
        
        for future in futures:
            result = future.result()
            if result:
                all_results.extend(result)
    
    save_to_csv(all_results, output_file)

if __name__ == "__main__":
    KRAKEN_FOLDER = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/Kraken2_results_v2"
    TARGET_FOLDER = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/Platon_results_annotated"
    OUTPUT_FILE = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/Kraken2_results_v2/blast_results_filtered.csv"
    
    main(KRAKEN_FOLDER, TARGET_FOLDER, OUTPUT_FILE)
