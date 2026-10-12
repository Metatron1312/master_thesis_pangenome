#!/usr/bin/env python3
"""
BATCH MAP AMR CONTIGS - CLEAN VERSION
ONLY creates: /output/{gene_id}/{gene_id}_mapped.gam
"""

import os
import sys
import subprocess
import glob
import time
import tempfile
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed

def setup_config():
    """Clean configuration - no all_contigs references"""
    return {
        'GRAPH_GBZ': "/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/vg_graph/graph/106_strains_with_refseq_complete_vg_autoindex/complete_106_strains_with_refseq_pangenome.giraffe.gbz",
        'MIN_INDEX': "/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/vg_graph/graph/106_strains_with_refseq_complete_vg_autoindex/complete_106_strains_with_refseq_pangenome.shortread.withzip.min",
        'DIST_INDEX': "/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/vg_graph/graph/106_strains_with_refseq_complete_vg_autoindex/complete_106_strains_with_refseq_pangenome.dist",
        
        'SAMPLE_ID': "SAMEA8218225",
        'ASSEMBLIES_ROOT': "/Users/shinra/SAMEA8218225/amr_gene_assemblies",
        'OUTPUT_ROOT': "/Users/shinra/SAMEA8218225/arm_contigs_mapped_to_graph",
        
        'THREADS': 8,
        'MAX_CONCURRENT': 4,
    }

def clean_output_directory(output_root):
    """Clean output directory - remove ALL existing content"""
    print(f"🧹 Cleaning output directory: {output_root}")
    
    if os.path.exists(output_root):
        # Remove everything in the directory
        for item in os.listdir(output_root):
            item_path = os.path.join(output_root, item)
            try:
                if os.path.isdir(item_path):
                    shutil.rmtree(item_path)
                    print(f"   Removed directory: {item}")
                else:
                    os.remove(item_path)
                    print(f"   Removed file: {item}")
            except Exception as e:
                print(f"   ⚠️  Could not remove {item}: {e}")
    else:
        os.makedirs(output_root, exist_ok=True)
    
    print("✅ Output directory cleaned")

def find_gene_directories(assemblies_root):
    """Find gene directories with contig files"""
    print(f"🔍 Scanning: {assemblies_root}")
    
    gene_dirs = []
    
    for item in os.listdir(assemblies_root):
        item_path = os.path.join(assemblies_root, item)
        
        if os.path.isdir(item_path):
            contig_files = glob.glob(os.path.join(item_path, "*_contigs.fa"))
            if contig_files:
                gene_dirs.append({
                    'path': item_path,
                    'name': item,
                    'contig_files': contig_files
                })
    
    print(f"✅ Found {len(gene_dirs)} gene directories")
    return gene_dirs

def map_single_gene(gene_info, config, job_id):
    """Map contigs for a single gene"""
    gene_name = gene_info['name']
    
    print(f"[{job_id}] 🧬 {gene_name}")
    
    # Create ONLY gene directory
    gene_output_dir = os.path.join(config['OUTPUT_ROOT'], gene_name)
    os.makedirs(gene_output_dir, exist_ok=True)
    
    # Create temporary FASTQ
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=True) as temp_fastq:
        try:
            # Convert contigs to FASTQ
            contig_count = 0
            for contig_file in gene_info['contig_files']:
                with open(contig_file, 'r') as f:
                    lines = f.readlines()
                
                # Parse FASTA
                header = ""
                sequence = ""
                for line in lines:
                    line = line.strip()
                    if line.startswith('>'):
                        if header and sequence:
                            # Write previous contig
                            contig_id = header[1:].split()[0]
                            temp_fastq.write(f"@{gene_name}_{contig_id}\n")
                            temp_fastq.write(f"{sequence}\n")
                            temp_fastq.write("+\n")
                            temp_fastq.write(f"{'F' * len(sequence)}\n")
                            contig_count += 1
                        
                        header = line
                        sequence = ""
                    else:
                        sequence += line
                
                # Write last contig
                if header and sequence:
                    contig_id = header[1:].split()[0]
                    temp_fastq.write(f"@{gene_name}_{contig_id}\n")
                    temp_fastq.write(f"{sequence}\n")
                    temp_fastq.write("+\n")
                    temp_fastq.write(f"{'F' * len(sequence)}\n")
                    contig_count += 1
            
            temp_fastq.flush()
            
            if contig_count == 0:
                return {'gene': gene_name, 'success': False, 'error': 'No contigs'}
            
            print(f"[{job_id}]   ✓ {contig_count} contigs")
            
        except Exception as e:
            return {'gene': gene_name, 'success': False, 'error': f'FASTQ: {e}'}
        
        # Map with vg giraffe
        gam_file = os.path.join(gene_output_dir, f"{gene_name}_mapped.gam")
        start_time = time.time()
        
        cmd = [
            "vg", "giraffe",
            "-Z", config['GRAPH_GBZ'],
            "-m", config['MIN_INDEX'],
            "-d", config['DIST_INDEX'],
            "-f", temp_fastq.name,
            "-t", str(config['THREADS']),
            "--align-from-chains"
        ]
        
        try:
            with open(gam_file, 'w') as out_f:
                process = subprocess.Popen(
                    cmd,
                    stdout=out_f,
                    stderr=subprocess.PIPE,
                    text=True
                )
                
                # Wait for completion
                stdout, stderr = process.communicate()
                mapping_time = time.time() - start_time
                
                if process.returncode != 0:
                    return {
                        'gene': gene_name,
                        'success': False,
                        'error': f"vg failed: {stderr[:100]}",
                        'time': mapping_time
                    }
                
                # Check output
                if not os.path.exists(gam_file) or os.path.getsize(gam_file) == 0:
                    return {
                        'gene': gene_name,
                        'success': False,
                        'error': "Empty GAM",
                        'time': mapping_time
                    }
                
                size_kb = os.path.getsize(gam_file) / 1024
                print(f"[{job_id}]   ✅ {size_kb:.1f} KB, {mapping_time:.1f}s")
                
                return {
                    'gene': gene_name,
                    'success': True,
                    'contigs': contig_count,
                    'time': mapping_time,
                    'size_kb': size_kb,
                    'gam_file': gam_file
                }
                
        except Exception as e:
            return {'gene': gene_name, 'success': False, 'error': f'Mapping: {e}'}

def main():
    print("=" * 80)
    print("BATCH MAP AMR CONTIGS - CLEAN VERSION")
    print("=" * 80)
    
    # Setup
    config = setup_config()
    
    print(f"\n📁 Sample: {config['SAMPLE_ID']}")
    print(f"📁 Input: {config['ASSEMBLIES_ROOT']}")
    print(f"📁 Output: {config['OUTPUT_ROOT']}")
    print(f"⚙️  Max concurrent: {config['MAX_CONCURRENT']}")
    
    # Clean output directory FIRST
    clean_output_directory(config['OUTPUT_ROOT'])
    
    # Check graph files
    print(f"\n🔍 Checking graph files...")
    for key in ['GRAPH_GBZ', 'MIN_INDEX', 'DIST_INDEX']:
        if not os.path.exists(config[key]):
            print(f"❌ Missing: {os.path.basename(config[key])}")
            sys.exit(1)
        print(f"   ✓ {os.path.basename(config[key])}")
    
    # Find genes
    gene_dirs = find_gene_directories(config['ASSEMBLIES_ROOT'])
    if not gene_dirs:
        print("❌ No genes found")
        sys.exit(1)
    
    print(f"\n📋 Found {len(gene_dirs)} genes")
    for i, gene in enumerate(gene_dirs[:5], 1):
        print(f"   {i}. {gene['name']} ({len(gene['contig_files'])} contigs)")
    if len(gene_dirs) > 5:
        print(f"   ... and {len(gene_dirs) - 5} more")
    
    # Confirm
    print(f"\n⚠️  Will map {len(gene_dirs)} genes")
    if input("Continue? (y/N): ").strip().lower() not in ['y', 'yes']:
        print("Aborted")
        sys.exit(0)
    
    # Process
    print(f"\n🚀 Starting...")
    start_time = time.time()
    results = []
    
    with ThreadPoolExecutor(max_workers=config['MAX_CONCURRENT']) as executor:
        # Submit jobs
        futures = {
            executor.submit(map_single_gene, gene, config, i+1): gene
            for i, gene in enumerate(gene_dirs)
        }
        
        # Process results
        completed = 0
        for future in as_completed(futures):
            gene = futures[future]
            try:
                result = future.result(timeout=300)
                results.append(result)
                completed += 1
                
                # Progress
                elapsed = time.time() - start_time
                if completed % 5 == 0 or completed == len(gene_dirs):
                    percent = completed / len(gene_dirs) * 100
                    print(f"\n📈 {completed}/{len(gene_dirs)} ({percent:.1f}%)")
                
            except Exception as e:
                print(f"\n❌ {gene['name']}: {e}")
                results.append({
                    'gene': gene['name'],
                    'success': False,
                    'error': str(e)
                })
                completed += 1
    
    total_time = time.time() - start_time
    
    # Results
    successful = [r for r in results if r['success']]
    failed = [r for r in results if not r['success']]
    
    print("\n" + "=" * 80)
    print("COMPLETED ✅")
    print("=" * 80)
    
    print(f"\n📊 RESULTS:")
    print(f"  • Total: {len(gene_dirs)} genes")
    print(f"  • Success: {len(successful)}")
    print(f"  • Failed: {len(failed)}")
    print(f"  • Time: {total_time:.1f}s")
    print(f"  • Output: {config['OUTPUT_ROOT']}")
    
    # Show structure
    print(f"\n📁 FINAL STRUCTURE:")
    print(f"  {config['OUTPUT_ROOT']}/")
    
    # List created directories
    created_dirs = [d for d in os.listdir(config['OUTPUT_ROOT']) 
                   if os.path.isdir(os.path.join(config['OUTPUT_ROOT'], d))]
    
    for i, dir_name in enumerate(sorted(created_dirs)[:3], 1):
        gam_path = os.path.join(config['OUTPUT_ROOT'], dir_name, f"{dir_name}_mapped.gam")
        if os.path.exists(gam_path):
            size_kb = os.path.getsize(gam_path) / 1024
            print(f"    ├── {dir_name}/")
            print(f"    │   └── {dir_name}_mapped.gam ({size_kb:.1f} KB)")
    
    if len(created_dirs) > 3:
        print(f"    └── ... {len(created_dirs) - 3} more directories")
    
    # Verify NO all_contigs directory
    all_contigs_path = os.path.join(config['OUTPUT_ROOT'], "all_contigs")
    if os.path.exists(all_contigs_path):
        print(f"\n⚠️  WARNING: Found 'all_contigs' directory (will be removed)")
        shutil.rmtree(all_contigs_path)
    
    print(f"\n✅ Ready for node extraction!")
    print(f"   Each gene: {config['OUTPUT_ROOT']}/{{gene_id}}/{{gene_id}}_mapped.gam")
    print("=" * 80)

if __name__ == "__main__":
    main()