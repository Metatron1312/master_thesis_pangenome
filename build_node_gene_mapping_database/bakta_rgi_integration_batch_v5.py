#!/usr/bin/env python3
"""
BAKTA-RGI INTEGRATION - SINGLE DATABASE VERSION
Tạo 1 file database duy nhất cho tất cả samples
"""

import pandas as pd
import re
from pathlib import Path
import logging
from datetime import datetime
import sys

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

def clean_value(value):
    """Clean string value"""
    if pd.isna(value) or str(value).strip().lower() == 'nan':
        return ''
    return str(value).strip()

def has_overlap(start1, end1, start2, end2, threshold=20):
    """Check if two positions overlap"""
    try:
        overlap = min(end1, end2) - max(start1, start2)
        return overlap >= threshold
    except:
        return False

def process_single_sample(sample_id, bakta_dir, rgi_dir):
    """Process a single sample and return DataFrame"""
    
    # File paths
    rgi_file = rgi_dir / sample_id / "rgi_results.txt"
    bakta_tsv = bakta_dir / sample_id / f"{sample_id}.tsv"
    bakta_log = bakta_dir / sample_id / f"{sample_id}.log"
    
    # 1. Parse contig mapping từ Bakta log
    contig_mapping = {}
    if bakta_log.exists():
        with open(bakta_log, 'r') as f:
            for line in f:
                if 'orig-id=' in line:
                    match = re.search(r'id=([^,]+),.*orig-id=([^,]+)', line)
                    if match:
                        bakta_contig, orig_contig = match.groups()
                        contig_mapping[bakta_contig] = orig_contig
    
    # 2. Load RGI data
    try:
        df_rgi = pd.read_csv(rgi_file, sep='\t', dtype=str)
    except Exception as e:
        logger.warning(f"  ❌ Cannot read RGI file for {sample_id}: {e}")
        return pd.DataFrame()
    
    # Process RGI entries
    rgi_entries = []
    for _, row in df_rgi.iterrows():
        # Clean contig name
        contig = re.sub(r'_\d+$', '', row['Contig'])
        
        # ARO number
        aro = clean_value(row['ARO'])
        if not aro or not aro.isdigit():
            continue
        
        # Strand
        orientation = clean_value(row['Orientation'])
        strand = '-' if orientation in ['-1', '-'] else '+'
        
        # Best hit và drug class
        best_hit = clean_value(row['Best_Hit_ARO'])
        drug_class = clean_value(row.get('Drug Class', ''))
        
        try:
            rgi_entries.append({
                'contig': contig,
                'start': int(clean_value(row['Start']) or 0),
                'end': int(clean_value(row['Stop']) or 0),
                'strand': strand,
                'aro': aro,
                'best_hit': best_hit,
                'drug_class': drug_class
            })
        except ValueError:
            continue
    
    # 3. Load Bakta genes
    bakta_genes = []
    if bakta_tsv.exists():
        with open(bakta_tsv, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 7 and parts[1] in ['CDS', 'cds']:
                    try:
                        gene_name = clean_value(parts[6] if len(parts) > 6 else '')
                        product = clean_value(parts[7] if len(parts) > 7 else '')
                        
                        bakta_genes.append({
                            'contig': parts[0],
                            'start': int(parts[2]),
                            'end': int(parts[3]),
                            'strand': parts[4],
                            'gene_id': parts[5],
                            'gene_name': gene_name,
                            'product': product
                        })
                    except ValueError:
                        continue
    
    # 4. Match RGI với Bakta và tạo results
    results = []
    matched_rgi_indices = set()
    
    # Match từng Bakta gene
    for bakta in bakta_genes:
        original_contig = contig_mapping.get(bakta['contig'], bakta['contig'])
        
        # Tìm RGI match
        matched_rgi = None
        matched_idx = -1
        
        for idx, rgi in enumerate(rgi_entries):
            if idx in matched_rgi_indices:
                continue
            
            if (rgi['contig'] == original_contig and 
                has_overlap(rgi['start'], rgi['end'], bakta['start'], bakta['end'])):
                matched_rgi = rgi
                matched_idx = idx
                break
        
        if matched_rgi:
            # BAKTA+RGI entry
            matched_rgi_indices.add(matched_idx)
            
            # Gene name: ưu tiên từ Bakta, nếu không có thì dùng từ RGI
            gene_name = bakta['gene_name'] if bakta['gene_name'] else matched_rgi['best_hit']
            
            results.append({
                'Sample_ID': sample_id,
                'Original_Contig': original_contig,
                'Bakta_Contig': bakta['contig'],
                'Start': matched_rgi['start'],
                'End': matched_rgi['end'],
                'Strand': matched_rgi['strand'],
                'GeneID': bakta['gene_id'],
                'Product': bakta['product'] or matched_rgi['best_hit'],
                'Gene_Name': gene_name,
                'AMR_Gene': matched_rgi['aro'],
                'Drug_Class': matched_rgi['drug_class'],
                'Source_Type': "BAKTA+RGI"
            })
        else:
            # BAKTA-only entry
            results.append({
                'Sample_ID': sample_id,
                'Original_Contig': original_contig,
                'Bakta_Contig': bakta['contig'],
                'Start': bakta['start'],
                'End': bakta['end'],
                'Strand': bakta['strand'],
                'GeneID': bakta['gene_id'],
                'Product': bakta['product'],
                'Gene_Name': bakta['gene_name'],
                'AMR_Gene': "",
                'Drug_Class': "",
                'Source_Type': "BAKTA"
            })
    
    # 5. Thêm RGI-only entries
    for idx, rgi in enumerate(rgi_entries):
        if idx in matched_rgi_indices:
            continue
        
        # Tìm Bakta contig
        bakta_contig = "UNMAPPED"
        for bakta_c, orig_c in contig_mapping.items():
            if orig_c == rgi['contig']:
                bakta_contig = bakta_c
                break
        
        gene_id = f"RGI_{rgi['aro']}_{rgi['start']}_{rgi['end']}"
        
        results.append({
            'Sample_ID': sample_id,
            'Original_Contig': rgi['contig'],
            'Bakta_Contig': bakta_contig,
            'Start': rgi['start'],
            'End': rgi['end'],
            'Strand': rgi['strand'],
            'GeneID': gene_id,
            'Product': rgi['best_hit'],
            'Gene_Name': rgi['best_hit'],
            'AMR_Gene': rgi['aro'],
            'Drug_Class': rgi['drug_class'],
            'Source_Type': "RGI_ONLY"
        })
    
    # Tạo DataFrame
    if results:
        df = pd.DataFrame(results)
        
        # Clean final values
        for col in ['Gene_Name', 'AMR_Gene', 'Drug_Class']:
            df[col] = df[col].apply(clean_value)
        
        return df
    else:
        return pd.DataFrame()

def main():
    """Main function - tạo 1 database duy nhất"""
    
    logger.info("=" * 80)
    logger.info("BAKTA-RGI INTEGRATION - SINGLE DATABASE")
    logger.info("=" * 80)
    
    # Paths
    bakta_dir = Path("/Users/shinra/Bioinformatics/WGS_AMR/data_public/Bakta_results")
    rgi_dir = Path("/Users/shinra/Bioinformatics/WGS_AMR/data_public/Rgi_results")
    output_dir = Path("/Users/shinra/Bioinformatics/WGS_AMR/data_public/Graph_annotation")
    
    # Tạo output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Output file (1 file duy nhất)
    output_file = output_dir / "bakta_rgi_integration_107_samples_v5.tsv"
    
    # Tìm tất cả samples
    logger.info("🔍 Scanning for samples...")
    samples = []
    for sample_dir in bakta_dir.iterdir():
        if sample_dir.is_dir():
            sample_id = sample_dir.name
            required_files = [
                sample_dir / f"{sample_id}.tsv",
                sample_dir / f"{sample_id}.log",
                rgi_dir / sample_id / "rgi_results.txt"
            ]
            if all(f.exists() for f in required_files):
                samples.append(sample_id)
    
    logger.info(f"📊 Found {len(samples)} samples with complete data")
    
    if not samples:
        logger.error("❌ No samples found!")
        return
    
    # Xử lý từng sample và collect kết quả
    all_data = []
    stats = {
        'processed': 0,
        'failed': 0,
        'failed_list': [],
        'total_entries': 0,
        'total_with_aro': 0
    }
    
    logger.info("\n🚀 Starting processing...")
    
    for i, sample_id in enumerate(samples, 1):
        logger.info(f"[{i}/{len(samples)}] Processing {sample_id}")
        
        try:
            df_sample = process_single_sample(sample_id, bakta_dir, rgi_dir)
            
            if not df_sample.empty:
                all_data.append(df_sample)
                stats['processed'] += 1
                
                # Sample statistics
                sample_entries = len(df_sample)
                sample_with_aro = len(df_sample[df_sample['AMR_Gene'] != ''])
                
                stats['total_entries'] += sample_entries
                stats['total_with_aro'] += sample_with_aro
                
                logger.info(f"   ✅ {sample_entries} entries ({sample_with_aro} with ARO)")
            else:
                stats['failed'] += 1
                stats['failed_list'].append(sample_id)
                logger.warning(f"   ⚠️ No data generated")
                
        except Exception as e:
            stats['failed'] += 1
            stats['failed_list'].append(sample_id)
            logger.error(f"   ❌ Error: {e}")
    
    # Combine all data
    if all_data:
        logger.info("\n💾 Combining all data into single database...")
        df_all = pd.concat(all_data, ignore_index=True)
        
        # Save to single file
        df_all.to_csv(output_file, sep='\t', index=False)
        
        # Final statistics
        total_entries = len(df_all)
        total_with_aro = len(df_all[df_all['AMR_Gene'] != ''])
        
        source_counts = df_all['Source_Type'].value_counts()
        
        logger.info("\n" + "=" * 80)
        logger.info("🏁 PROCESSING COMPLETED!")
        logger.info("=" * 80)
        logger.info(f"📊 Samples processed: {stats['processed']}/{len(samples)}")
        logger.info(f"📁 Total entries: {total_entries:,}")
        logger.info(f"💊 Entries with ARO numbers: {total_with_aro:,}")
        
        logger.info("\n📈 SOURCE TYPE DISTRIBUTION:")
        for source, count in source_counts.items():
            percentage = (count / total_entries) * 100
            logger.info(f"   {source}: {count:,} ({percentage:.1f}%)")
        
        # Check specific sample
        check_sample = "SAMD00178086"
        if check_sample in samples:
            df_check = df_all[df_all['Sample_ID'] == check_sample]
            check_aro = len(df_check[df_check['AMR_Gene'] != ''])
            logger.info(f"\n🔍 CHECK {check_sample}:")
            logger.info(f"   Entries: {len(df_check):,}")
            logger.info(f"   With ARO: {check_aro:,}")
            
            if check_aro > 0:
                aro_entries = df_check[df_check['AMR_Gene'] != ''].head(3)
                for _, row in aro_entries.iterrows():
                    logger.info(f"     - {row['Original_Contig']}:{row['Start']}-{row['End']}")
                    logger.info(f"       Gene: {row['Gene_Name'][:50]}...")
                    logger.info(f"       ARO: {row['AMR_Gene']}")
        
        if stats['failed'] > 0:
            logger.info(f"\n❌ Failed samples: {stats['failed']}")
            for sample in stats['failed_list'][:10]:
                logger.info(f"   - {sample}")
        
        logger.info(f"\n💾 SINGLE DATABASE FILE: {output_file}")
        logger.info("=" * 80)
        
        # Tạo summary đơn giản
        summary_file = output_dir / "processing_summary.txt"
        with open(summary_file, 'w') as f:
            f.write("BAKTA-RGI INTEGRATION SUMMARY\n")
            f.write("=" * 50 + "\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Total samples: {len(samples)}\n")
            f.write(f"Processed successfully: {stats['processed']}\n")
            f.write(f"Failed: {stats['failed']}\n")
            f.write(f"Total entries: {total_entries}\n")
            f.write(f"Entries with ARO: {total_with_aro}\n")
            f.write(f"Database file: {output_file}\n")
        
        logger.info(f"📝 Summary saved to: {summary_file}")
        
    else:
        logger.error("❌ No data generated for any sample!")

if __name__ == "__main__":
    main()