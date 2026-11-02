import os
import glob
import time
import logging
from datetime import datetime
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
import socket

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('blast_analysis.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Set socket timeout globally
socket.setdefaulttimeout(60)  # 60 seconds timeout for all socket operations

def parse_kraken_output(kraken_folder):
    """Parse Kraken2 output files to find unclassified contigs"""
    unclassified_contigs = {}
    processed_files = 0
    
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
        
        processed_files += 1
        contig_count = 0
        
        with open(file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3 and parts[0] == 'U':
                    contig_id = parts[1]
                    key = (biosample_id, seq_type)
                    if key not in unclassified_contigs:
                        unclassified_contigs[key] = []
                    unclassified_contigs[key].append(contig_id)
                    contig_count += 1
        
        logger.info(f"Processed {filename}: found {contig_count} unclassified contigs")
    
    logger.info(f"Total files processed: {processed_files}")
    logger.info(f"Total samples with unclassified contigs: {len(unclassified_contigs)}")
    return unclassified_contigs

def find_fasta_file(target_folder, biosample_id, seq_type):
    """Find fasta file for given biosample and sequence type"""
    sample_dir = os.path.join(target_folder, biosample_id)
    if not os.path.exists(sample_dir):
        logger.warning(f"Directory not found: {sample_dir}")
        return None
    
    patterns = [
        f"*_contigs.with_metadata.{seq_type}.fasta",
        f"*.{seq_type}.fasta",
        f"*{seq_type}*.fasta"
    ]
    
    for pattern in patterns:
        fasta_files = glob.glob(os.path.join(sample_dir, pattern))
        if fasta_files:
            logger.debug(f"Found {seq_type} fasta file: {fasta_files[0]}")
            return fasta_files[0]
    
    logger.warning(f"No {seq_type} fasta file found for {biosample_id} in {sample_dir}")
    return None

def extract_contig_sequence(fasta_file, contig_id):
    """Extract contig sequence from fasta file"""
    if not fasta_file or not os.path.exists(fasta_file):
        logger.error(f"Fasta file not found or inaccessible: {fasta_file}")
        return None
    
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id == contig_id:
                logger.debug(f"Found contig {contig_id} in {fasta_file}")
                return record
        logger.warning(f"Contig {contig_id} not found in {fasta_file}")
        return None
    except Exception as e:
        logger.error(f"Error parsing fasta file {fasta_file}: {str(e)}")
        return None

def run_blast(sequence, query_id, max_retries=3):
    """Run BLAST search with retry mechanism"""
    for attempt in range(max_retries):
        try:
            logger.info(f"Running BLAST for {query_id} (attempt {attempt+1})...")
            
            # Remove timeout parameter as it's not supported in some Biopython versions
            result_handle = NCBIWWW.qblast(
                program="blastn",
                database="nt",
                sequence=sequence.seq,
                hitlist_size=10,
                entrez_query='none',
                format_type="XML"
            )
            
            blast_records = list(NCBIXML.parse(result_handle))
            logger.info(f"BLAST successful for {query_id}")
            return blast_records
            
        except Exception as e:
            logger.warning(f"Attempt {attempt+1} failed for {query_id}: {str(e)}")
            if attempt == max_retries - 1:
                logger.error(f"BLAST failed for {query_id} after {max_retries} attempts")
                return None
            time.sleep(10 * (attempt + 1))  # Exponential backoff for retry

def filter_blast_hits(hits, query_length):
    """Filter BLAST results based on quality criteria"""
    filtered = []
    
    for hit in hits:
        try:
            coverage = hit['Alignment_Length'] / query_length
            hit['Query_Length'] = query_length
            hit['Coverage'] = round(coverage * 100, 2)
            
            if (hit['E_value'] < 0.001 and
                hit['Bit_Score'] > 50 and
                hit['Percentage_Identity'] > 70 and
                coverage > 0.7):
                
                filtered.append(hit)
        except KeyError as e:
            logger.error(f"Missing key in hit data: {str(e)}")
            continue
    
    filtered.sort(key=lambda x: x['Bit_Score'], reverse=True)
    return filtered[:3]

def process_contig(biosample_id, seq_type, contig_id, target_folder):
    """Process a single contig"""
    logger.info(f"Processing {biosample_id} {seq_type} {contig_id}")
    
    fasta_file = find_fasta_file(target_folder, biosample_id, seq_type)
    if not fasta_file:
        return None
    
    contig_seq = extract_contig_sequence(fasta_file, contig_id)
    if not contig_seq:
        return None
    
    query_id = f"{biosample_id}_{seq_type}_{contig_id}"
    blast_results = run_blast(contig_seq, query_id)
    
    if not blast_results:
        return None
    
    hits = []
    for blast_record in blast_results:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                try:
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
                except Exception as e:
                    logger.error(f"Error processing alignment for {query_id}: {str(e)}")
                    continue
    
    return filter_blast_hits(hits, len(contig_seq))

def save_to_csv(results, output_file):
    """Save results to CSV file"""
    if not results:
        logger.warning("No results to save.")
        return
    
    try:
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        df = pd.DataFrame(results)
        
        columns_order = [
            'Query_ID', 'Sample_ID', 'Sequence_Type', 'Contig_ID', 'Query_Length',
            'Hit_ID', 'Accession', 'Description',
            'E_value', 'Bit_Score', 'Percentage_Identity', 'Identity', 'Alignment_Length', 'Coverage',
            'Query_Start', 'Query_End', 'Subject_Start', 'Subject_End',
            'Alignment'
        ]
        
        columns_order = [col for col in columns_order if col in df.columns]
        df = df[columns_order]
        
        df.to_csv(output_file, index=False)
        logger.info(f"Successfully saved {len(df)} results to {output_file}")
    except Exception as e:
        logger.error(f"Failed to save results: {str(e)}")
        raise

def main(kraken_folder, target_folder, output_file, max_workers=4):
    """Main analysis function"""
    start_time = datetime.now()
    logger.info(f"Starting analysis at {start_time}")
    
    try:
        unclassified_contigs = parse_kraken_output(kraken_folder)
        
        if not unclassified_contigs:
            logger.warning("No unclassified contigs found in Kraken2 output files.")
            return
        
        total_contigs = sum(len(v) for v in unclassified_contigs.values())
        logger.info(f"Total unclassified contigs to process: {total_contigs}")
        
        all_results = []
        processed = 0
        failed = 0
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
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
            
            for future in as_completed(futures):
                try:
                    result = future.result()
                    if result:
                        all_results.extend(result)
                        processed += 1
                    else:
                        failed += 1
                except Exception as e:
                    logger.error(f"Error processing contig: {str(e)}")
                    failed += 1
                
                if (processed + failed) % 10 == 0:
                    logger.info(
                        f"Progress: {processed + failed}/{total_contigs} "
                        f"({(processed + failed)/total_contigs:.1%}), "
                        f"Success: {processed}, Failed: {failed}"
                    )
        
        save_to_csv(all_results, output_file)
        
    except Exception as e:
        logger.error(f"Fatal error in main: {str(e)}")
        raise
    finally:
        duration = datetime.now() - start_time
        logger.info(f"Analysis completed in {duration}")

if __name__ == "__main__":
    try:
        # Configure paths
        KRAKEN_FOLDER = "/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur/kraken2_results"
        TARGET_FOLDER = "/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur/platon_results"
        OUTPUT_FILE = "/Users/shinra/Bioinformatics/WGS_AMR/data_pasteur/filtered_contigs/blast_results_filtered.csv"
        
        # Run analysis
        main(
            kraken_folder=KRAKEN_FOLDER,
            target_folder=TARGET_FOLDER,
            output_file=OUTPUT_FILE,
            max_workers=4
        )
    except Exception as e:
        logger.critical(f"Script failed: {str(e)}")