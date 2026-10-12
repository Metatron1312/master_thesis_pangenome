#!/usr/bin/env python3
"""
AMR GENE READ EXTRACTION - GENERALIZED ADAPTIVE VERSION (LOWER THRESHOLDS)
Xử lý tất cả AMR genes với adaptive thresholds - ĐÃ HẠ THRESHOLD ĐỂ GIỮ NHIỀU GENE
"""
import os
import subprocess
import csv
import glob
import pysam
from collections import defaultdict
import json
import time
import logging
from pathlib import Path

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("amr_gene_extraction_general.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class AMRGeneExtractorGeneral:
    def __init__(self, config):
        self.config = config
        self.setup_directories()
        self.validation_report = {}
        
        # GENERALIZED ADAPTIVE THRESHOLDS
        self.plasmid_patterns = self.config.get('plasmid_patterns', ['plasmid', 'circular'])
        
    def setup_directories(self):
        """Create necessary directories"""
        os.makedirs(self.config['output_dir'], exist_ok=True)
        os.makedirs(f"{self.config['output_dir']}/validation_reports", exist_ok=True)
        os.makedirs(f"{self.config['output_dir']}/amr_genes_only", exist_ok=True)
    
    def get_contig_type(self, contig_name):
        """Xác định contig type để áp dụng adaptive thresholds"""
        contig_lower = contig_name.lower()
        
        # Kiểm tra nếu là plasmid
        for pattern in self.plasmid_patterns:
            if pattern in contig_lower:
                return "plasmid"
        
        # Mặc định là chromosome
        return "chromosome"
    
    def get_gene_length_category(self, gene_length):
        """Phân loại gene dựa trên length"""
        if gene_length < 500:
            return "short"
        elif gene_length < 1500:
            return "medium"
        else:
            return "long"
    
    def get_adaptive_thresholds(self, contig_type, gene_length, total_reads_in_region):
        """Lấy thresholds phù hợp dựa trên contig type, gene length và coverage"""
        
        # BASE THRESHOLDS - ĐÃ GIẢM CHO CHROMOSOME
        if contig_type == "plasmid":
            # PLASMID: thresholds thấp
            base_min_coverage = 0.30
            base_min_reads = 3
            fetch_buffer = 1000
        else:
            # CHROMOSOME: ĐÃ GIẢM TỪ 0.50 → 0.35 → 0.15
            base_min_coverage = 0.15      # <--- ĐÃ GIẢM: từ 0.35 xuống 0.15
            base_min_reads = 2             # <--- ĐÃ GIẢM: từ 5 xuống 2
            fetch_buffer = 500
        
        # Điều chỉnh dựa trên gene length
        length_category = self.get_gene_length_category(gene_length)
        if length_category == "short":
            # Genes ngắn: coverage threshold cao hơn (vì khó align đủ)
            min_coverage = base_min_coverage + 0.10
        elif length_category == "long":
            # Genes dài: coverage threshold thấp hơn
            min_coverage = max(0.10, base_min_coverage - 0.10)  # <--- ĐÃ GIẢM: từ 0.20 xuống 0.10
        else:
            min_coverage = base_min_coverage
        
        # Điều chỉnh dựa trên tổng số reads trong region
        if total_reads_in_region < 1000:
            # Low coverage region: giảm threshold
            min_reads = max(1, base_min_reads - 1)  # <--- ĐÃ GIẢM: từ max(2) xuống max(1)
        elif total_reads_in_region > 10000:
            # High coverage region: tăng threshold
            min_reads = base_min_reads + 5
        else:
            min_reads = base_min_reads
        
        return {
            'min_read_coverage': min(min_coverage, 0.80),  # Không quá 80%
            'min_reads_per_gene': min_reads,
            'fetch_buffer': fetch_buffer
        }
    
    def load_amr_genes(self):
        """Load tất cả AMR genes"""
        logger.info("Loading all AMR genes from database...")
        
        amr_genes = []
        
        with open(self.config['bakta_rgi_db'], 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if row['Sample_ID'] == self.config['sample_id']:
                    gene_id = row['GeneID']
                    amr_gene = row['AMR_Gene'].strip()
                    
                    # Chỉ lấy genes có AMR_Gene field hợp lệ
                    if (amr_gene and 
                        amr_gene.lower() != 'nan' and 
                        amr_gene.upper() != 'N/A' and
                        any(c.isdigit() for c in amr_gene)):
                        
                        contig = row['Original_Contig']
                        product = row.get('Product', '')
                        gene_name = row.get('Gene_Name', '')
                        gene_length = int(row['End']) - int(row['Start']) + 1
                        
                        # Xác định contig type
                        contig_type = self.get_contig_type(contig)
                        
                        gene_info = {
                            'gene_id': gene_id,
                            'contig': contig,
                            'start': int(row['Start']),
                            'end': int(row['End']),
                            'strand': row['Strand'],
                            'product': product,
                            'amr_gene': amr_gene,
                            'source_type': row['Source_Type'],
                            'gene_name': gene_name,
                            'drug_class': row['Drug_Class'],
                            'gene_length': gene_length,
                            'contig_type': contig_type
                        }
                        amr_genes.append(gene_info)
        
        logger.info(f"Loaded {len(amr_genes)} AMR genes for sample {self.config['sample_id']}")
        
        # Thống kê
        plasmid_count = sum(1 for g in amr_genes if g['contig_type'] == 'plasmid')
        logger.info(f"  - Plasmid genes: {plasmid_count}")
        logger.info(f"  - Chromosome genes: {len(amr_genes) - plasmid_count}")
        
        return amr_genes
    
    def map_reads_to_contigs(self):
        """Map reads to contigs - giữ nguyên như cũ"""
        logger.info("Step 1: Mapping reads to contigs...")
        
        bam_file = f"{self.config['output_dir']}/{self.config['sample_id']}_to_contigs.bam"
        
        if os.path.exists(bam_file) and os.path.exists(bam_file + ".bai"):
            logger.info("BAM file already exists, skipping mapping")
            return bam_file
        
        # Run minimap2
        sam_file = f"{self.config['output_dir']}/{self.config['sample_id']}.sam"
        minimap2_cmd = [
            "minimap2", "-ax", "sr", "-t", "8", "--secondary=no", "-Y", "--eqx",
            self.config['contigs_file'],
            f"{self.config['reads_dir']}/{self.config['sample_id']}_1_trimmed.fastq.gz",
            f"{self.config['reads_dir']}/{self.config['sample_id']}_2_trimmed.fastq.gz",
            "-o", sam_file
        ]
        
        logger.info("Running minimap2...")
        result = subprocess.run(minimap2_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"Minimap2 failed: {result.stderr}")
            raise RuntimeError("Minimap2 mapping failed")
        
        # Convert SAM to BAM and sort
        logger.info("Converting SAM to BAM...")
        subprocess.run(["samtools", "view", "-b", "-o", bam_file, sam_file], check=True)
        subprocess.run(["samtools", "sort", "-o", bam_file, bam_file], check=True)
        subprocess.run(["samtools", "index", bam_file], check=True)
        
        # Remove temporary SAM file
        os.remove(sam_file)
        
        logger.info("Read mapping completed successfully")
        return bam_file
    
    def calculate_overlap(self, read_start, read_end, gene_start, gene_end):
        """Calculate overlap between read and gene region"""
        return max(0, min(read_end, gene_end) - max(read_start, gene_start))
    
    def calculate_read_coverage_in_gene(self, read, gene_start, gene_end):
        """Calculate what percentage of READ LENGTH aligns within gene boundaries"""
        read_start = read.reference_start
        read_end = read.reference_end
        read_length = read.query_length
        
        if read_length == 0:
            return 0
        
        # Calculate overlap with gene
        overlap = self.calculate_overlap(read_start, read_end, gene_start, gene_end)
        
        # Return ratio of READ that overlaps with gene
        return overlap / read_length
    
    def write_reads_to_fastq(self, reads, output_file):
        """Write reads to FASTQ file"""
        with open(output_file, 'w') as f:
            for read in reads:
                if read.is_read1:
                    read_id = f"@{read.query_name}/1"
                else:
                    read_id = f"@{read.query_name}/2"
                
                f.write(f"{read_id}\n")
                f.write(f"{read.query_sequence}\n")
                f.write("+\n")
                f.write(f"{read.qual}\n")
    
    def estimate_total_reads_in_region(self, bam, contig, start, end, buffer=1000):
        """Ước lượng tổng số reads trong region để điều chỉnh thresholds"""
        try:
            fetch_start = max(0, start - buffer)
            fetch_end = end + buffer
            count = 0
            for _ in bam.fetch(contig, fetch_start, fetch_end):
                count += 1
                if count >= 10000:  # Chỉ đếm đến 10000 để tiết kiệm thời gian
                    break
            return count
        except:
            return 1000  # Default value
    
    def extract_amr_gene_reads(self, bam_file, amr_genes):
        """Step 2: Extract reads cho tất cả AMR genes với adaptive thresholds"""
        logger.info(f"Step 2: Extracting reads for {len(amr_genes)} AMR genes...")
        
        bam = pysam.AlignmentFile(bam_file, "rb")
        processed_genes = 0
        successful_extractions = 0
        
        # Statistics
        stats = {
            'plasmid': {'total': 0, 'success': 0, 'reads_extracted': 0},
            'chromosome': {'total': 0, 'success': 0, 'reads_extracted': 0}
        }
        
        for gene_info in amr_genes:
            gene_id = gene_info['gene_id']
            contig = gene_info['contig']
            start = gene_info['start']
            end = gene_info['end']
            gene_length = gene_info['gene_length']
            contig_type = gene_info['contig_type']
            
            # Ước lượng tổng reads trong region
            total_reads_estimate = self.estimate_total_reads_in_region(bam, contig, start, end)
            
            # Lấy adaptive thresholds
            thresholds = self.get_adaptive_thresholds(contig_type, gene_length, total_reads_estimate)
            min_read_coverage = thresholds['min_read_coverage']
            min_reads_per_gene = thresholds['min_reads_per_gene']
            fetch_buffer = thresholds['fetch_buffer']
            
            # Update statistics
            stats[contig_type]['total'] += 1
            
            # Create gene directory
            gene_dir = f"{self.config['output_dir']}/amr_genes_only/{gene_id}"
            os.makedirs(gene_dir, exist_ok=True)
            
            r1_reads = []
            r2_reads = []
            total_reads_in_region = 0
            passed_reads = 0
            
            try:
                # Fetch reads từ vùng gene
                fetch_start = max(0, start - fetch_buffer)
                fetch_end = end + fetch_buffer
                
                for read in bam.fetch(contig, fetch_start, fetch_end):
                    total_reads_in_region += 1
                    
                    # Skip unmapped hoặc low quality reads
                    if read.is_unmapped or read.mapping_quality < self.config['min_mapping_quality']:
                        continue
                    
                    # Tính coverage
                    read_coverage = self.calculate_read_coverage_in_gene(read, start, end)
                    
                    # Lấy read nếu ≥ threshold
                    if read_coverage >= min_read_coverage:
                        passed_reads += 1
                        if read.is_read1:
                            r1_reads.append(read)
                        elif read.is_read2:
                            r2_reads.append(read)
                
                # Write reads
                if r1_reads:
                    self.write_reads_to_fastq(r1_reads, f"{gene_dir}/{gene_id}_R1.fastq")
                if r2_reads:
                    self.write_reads_to_fastq(r2_reads, f"{gene_dir}/{gene_id}_R2.fastq")
                
                # Store validation info
                self.validation_report[gene_id] = {
                    'amr_gene': gene_info['amr_gene'],
                    'total_reads_in_region': total_reads_in_region,
                    'passed_reads': passed_reads,
                    'r1_reads': len(r1_reads),
                    'r2_reads': len(r2_reads),
                    'total_extracted_reads': len(r1_reads) + len(r2_reads),
                    'gene_length': gene_length,
                    'contig': contig,
                    'contig_type': contig_type,
                    'position': f"{start}-{end}",
                    'product': gene_info['product'],
                    'gene_name': gene_info['gene_name'],
                    'drug_class': gene_info['drug_class'],
                    'source_type': gene_info['source_type'],
                    'min_read_coverage_used': min_read_coverage,
                    'min_reads_per_gene_required': min_reads_per_gene,
                    'total_reads_estimate': total_reads_estimate,
                    'extraction_method': 'GENERALIZED_ADAPTIVE'
                }
                
                # Kiểm tra nếu extraction thành công
                total_extracted = len(r1_reads) + len(r2_reads)
                if total_extracted >= min_reads_per_gene:
                    successful_extractions += 1
                    stats[contig_type]['success'] += 1
                    stats[contig_type]['reads_extracted'] += total_extracted
                    status = "SUCCESS"
                else:
                    status = "LOW_COVERAGE"
                
                processed_genes += 1
                
                # Hiển thị progress mỗi 10 genes
                if processed_genes % 10 == 0:
                    logger.info(f"  [{processed_genes}/{len(amr_genes)}] {gene_id}: {status}")
                    logger.info(f"     Type: {contig_type}, Length: {gene_length}bp")
                    logger.info(f"     Passed: {passed_reads}/{total_reads_in_region} reads")
                    logger.info(f"     Threshold: {min_read_coverage*100}% coverage")
                    
            except Exception as e:
                logger.error(f"Error processing AMR gene {gene_id}: {e}")
                self.validation_report[gene_id] = {
                    'amr_gene': gene_info['amr_gene'],
                    'error': str(e),
                    'contig_type': contig_type,
                    'extraction_method': 'FAILED'
                }
        
        bam.close()
        
        # Summary statistics
        logger.info(f"\n📊 EXTRACTION SUMMARY:")
        for contig_type, stat in stats.items():
            if stat['total'] > 0:
                success_rate = (stat['success'] / stat['total']) * 100
                logger.info(f"  {contig_type.upper()}:")
                logger.info(f"    Total genes: {stat['total']}")
                logger.info(f"    Successfully extracted: {stat['success']} ({success_rate:.1f}%)")
                logger.info(f"    Reads extracted: {stat['reads_extracted']:,}")
        
        logger.info(f"\nCompleted read extraction for {processed_genes} AMR genes")
        logger.info(f"Overall success rate: {successful_extractions/processed_genes*100:.1f}%")
        
        return successful_extractions
    
    def validate_extracted_reads(self):
        """Validate extracted reads"""
        logger.info("Step 3: Validating extracted reads...")
        
        validation_results = {}
        
        for gene_id, report in self.validation_report.items():
            if 'error' in report:
                validation_results[gene_id] = {'status': 'ERROR', 'reason': report['error']}
                continue
            
            gene_dir = f"{self.config['output_dir']}/amr_genes_only/{gene_id}"
            r1_file = f"{gene_dir}/{gene_id}_R1.fastq"
            r2_file = f"{gene_dir}/{gene_id}_R2.fastq"
            
            # Check if files exist
            r1_exists = os.path.exists(r1_file) and os.path.getsize(r1_file) > 0
            r2_exists = os.path.exists(r2_file) and os.path.getsize(r2_file) > 0
            
            min_reads_required = report.get('min_reads_per_gene_required', 10)
            
            if not r1_exists and not r2_exists:
                validation_results[gene_id] = {'status': 'FAILED', 'reason': 'No reads extracted'}
            elif report['r1_reads'] < min_reads_required or report['r2_reads'] < min_reads_required:
                validation_results[gene_id] = {
                    'status': 'LOW_COVERAGE', 
                    'reason': f"Low read count (R1: {report['r1_reads']}, R2: {report['r2_reads']})"
                }
            else:
                validation_results[gene_id] = {'status': 'PASSED'}
            
            # Update report
            report['validation_status'] = validation_results[gene_id]['status']
            if 'reason' in validation_results[gene_id]:
                report['validation_reason'] = validation_results[gene_id]['reason']
        
        return validation_results
    
    def generate_summary_report(self):
        """Generate comprehensive summary report"""
        logger.info("Generating summary report...")
        
        # Calculate statistics
        total_genes = len(self.validation_report)
        passed_genes = sum(1 for r in self.validation_report.values() 
                          if r.get('validation_status') == 'PASSED')
        failed_genes = total_genes - passed_genes
        
        # Read distribution
        total_reads = sum(r.get('total_extracted_reads', 0) for r in self.validation_report.values())
        
        # Phân tích theo contig type
        type_stats = defaultdict(lambda: {'total': 0, 'passed': 0, 'reads': 0})
        
        for report in self.validation_report.values():
            if 'error' in report:
                continue
            contig_type = report.get('contig_type', 'unknown')
            status = report.get('validation_status', 'UNKNOWN')
            
            type_stats[contig_type]['total'] += 1
            if status == 'PASSED':
                type_stats[contig_type]['passed'] += 1
            type_stats[contig_type]['reads'] += report.get('total_extracted_reads', 0)
        
        summary = {
            'processing_date': time.strftime('%Y-%m-%d %H:%M:%S'),
            'sample_id': self.config['sample_id'],
            'total_amr_genes': total_genes,
            'genes_passed': passed_genes,
            'genes_failed': failed_genes,
            'success_rate': (passed_genes / total_genes * 100) if total_genes > 0 else 0,
            'total_reads_extracted': total_reads,
            'stats_by_contig_type': dict(type_stats),
            'configuration': self.config
        }
        
        # Save reports
        summary_file = f"{self.config['output_dir']}/validation_reports/summary_report.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        detailed_file = f"{self.config['output_dir']}/validation_reports/detailed_report.json"
        with open(detailed_file, 'w') as f:
            json.dump(self.validation_report, f, indent=2)
        
        # Gene list file
        gene_list_file = f"{self.config['output_dir']}/validation_reports/gene_list.tsv"
        with open(gene_list_file, 'w') as f:
            f.write("gene_id\tamr_gene\tgene_name\tcontig_type\t")
            f.write("min_coverage\trequired_reads\textracted_reads\tstatus\n")
            
            for gene_id, report in self.validation_report.items():
                if 'error' in report:
                    continue
                f.write(f"{gene_id}\t{report.get('amr_gene', '')}\t{report.get('gene_name', '')}\t")
                f.write(f"{report.get('contig_type', '')}\t{report.get('min_read_coverage_used', 0)}\t")
                f.write(f"{report.get('min_reads_per_gene_required', 0)}\t")
                f.write(f"{report.get('total_extracted_reads', 0)}\t{report.get('validation_status', 'UNKNOWN')}\n")
        
        return summary
    
    def print_final_summary(self, summary):
        """Print final summary"""
        logger.info("=" * 80)
        logger.info("AMR GENE READ EXTRACTION - GENERALIZED VERSION (LOW THRESHOLD)")
        logger.info("=" * 80)
        logger.info(f"Sample ID: {summary['sample_id']}")
        logger.info(f"Total AMR genes processed: {summary['total_amr_genes']}")
        logger.info(f"Genes passed validation: {summary['genes_passed']}")
        logger.info(f"Success rate: {summary['success_rate']:.1f}%")
        logger.info(f"Total reads extracted: {summary['total_reads_extracted']:,}")
        
        if 'stats_by_contig_type' in summary:
            logger.info(f"\n📊 STATISTICS BY CONTIG TYPE:")
            for contig_type, stats in summary['stats_by_contig_type'].items():
                if stats['total'] > 0:
                    pass_rate = (stats['passed'] / stats['total']) * 100
                    logger.info(f"  {contig_type.upper()}:")
                    logger.info(f"    Genes: {stats['total']} total, {stats['passed']} passed ({pass_rate:.1f}%)")
                    logger.info(f"    Reads extracted: {stats['reads']:,}")
        
        logger.info("=" * 80)
    
    def run_pipeline(self):
        """Run the complete pipeline"""
        logger.info("Starting AMR Gene Read Extraction Pipeline - Generalized Version (Low Threshold)")
        start_time = time.time()
        
        try:
            # Step 1: Load all AMR genes
            amr_genes = self.load_amr_genes()
            
            if not amr_genes:
                logger.error("No AMR genes found for this sample!")
                return
            
            # Step 2: Map reads to contigs
            bam_file = self.map_reads_to_contigs()
            
            # Step 3: Extract reads với adaptive thresholds
            successful_extractions = self.extract_amr_gene_reads(bam_file, amr_genes)
            
            # Step 4: Validate results
            self.validate_extracted_reads()
            
            # Step 5: Generate reports
            summary = self.generate_summary_report()
            
            # Step 6: Final summary
            processing_time = time.time() - start_time
            summary['processing_time_seconds'] = round(processing_time, 2)
            summary['successful_extractions'] = successful_extractions
            
            self.print_final_summary(summary)
            
            logger.info(f"Pipeline completed in {processing_time:.2f} seconds")
            
        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            raise

def main():
    """Main function - generalized version"""
    # Configuration
    config = {
        'sample_id': "SAMEA8218225",  # Thay đổi sample ID tại đây
        'reads_dir': "/Users/shinra/SAMEA8218225/fastp_trimmed",
        'contigs_file': "/Users/shinra/Bioinformatics/WGS_AMR/data_public/combined_contigs_and_refseq_v2/SAMEA8218225.filtered.complete.fasta",
        'bakta_rgi_db': "/Users/shinra/Bioinformatics/WGS_AMR/data_public/Graph_annotation/bakta_rgi_integration_107_samples_v5.tsv",
        'output_dir': "/Users/shinra/SAMEA8218225/amr_gene_reads_v2",
        
        # Quality thresholds
        'min_mapping_quality': 20,
        
        # Plasmid detection patterns
        'plasmid_patterns': ['plasmid', 'circular'],
        
        # Adaptive thresholds configuration
        'adaptive_thresholds': True
    }
    
    # Run pipeline
    extractor = AMRGeneExtractorGeneral(config)
    extractor.run_pipeline()

if __name__ == "__main__":
    main()