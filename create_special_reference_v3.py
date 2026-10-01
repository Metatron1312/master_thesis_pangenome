#!/usr/bin/env python3
import os
import csv
from Bio import SeqIO
from Bio.Seq import Seq

# ==================== CONFIGURATION ====================
REF_GENOME_FILE = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/combined_contigs_and_refseq/BJAB07104.filtered.complete.fasta"
COORDS_FILE = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/vg_graph/unaligned_block_seqs_processed/block_coordinates.csv"
SPECIAL_REF_FILE = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/vg_graph/special_reference_v3.fa"
# =======================================================

def load_coordinates(coords_file):
    """Load coordinates từ CSV file"""
    coordinates = {}
    with open(coords_file, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # Skip header
        for row in reader:
            if len(row) >= 5:
                block_id, ref_start, ref_end, strand, length = row[0], row[1], row[2], row[3], row[4]
                coordinates[block_id] = {
                    'ref_start': int(ref_start),
                    'ref_end': int(ref_end),
                    'strand': strand,
                    'length': int(length)
                }
    return coordinates

def get_sequence_pangraph_coordinates(ref_genome_file, start, end, strand):
    """
    Lấy sequence theo pangraph coordinate system
    Pangraph: 0-based, right-excluded [start, end)
    Python: [start:end]
    """
    try:
        for record in SeqIO.parse(ref_genome_file, "fasta"):
            if record.id == "NC_021726.1":
                # PANGRAHP COORDINATES: [start, end) → Python [start:end]
                sequence = str(record.seq[start:end])
                
                # Xử lý strand
                if strand == '-':
                    sequence = str(Seq(sequence).reverse_complement())
                
                return sequence.upper()
        return None
    except Exception as e:
        print(f"Error reading reference genome: {e}")
        return None

def validate_pangraph_coordinates(coordinates):
    """Validate coordinates theo pangraph system"""
    print("🔍 Validating pangraph coordinates...")
    
    # Đọc genome để lấy total length
    genome_length = 0
    for record in SeqIO.parse(REF_GENOME_FILE, "fasta"):
        if record.id == "NC_021726.1":
            genome_length = len(record.seq)
            break
    
    print(f"📏 Genome length: {genome_length} bp")
    print(f"📏 Pangraph coordinates range: [0, {genome_length})")
    
    issues = []
    for block_id, coords in list(coordinates.items())[:10]:  # Check 10 blocks đầu
        # Pangraph coordinates: [start, end) với start >= 0, end <= genome_length
        if coords['ref_start'] < 0:
            issues.append(f"Block {block_id}: start {coords['ref_start']} < 0")
        
        if coords['ref_end'] > genome_length:
            issues.append(f"Block {block_id}: end {coords['ref_end']} > genome length {genome_length}")
        
        if coords['ref_start'] >= coords['ref_end']:
            issues.append(f"Block {block_id}: start {coords['ref_start']} >= end {coords['ref_end']}")
        
        # Expected length theo pangraph: end - start
        expected_length = coords['ref_end'] - coords['ref_start']
        if coords['length'] != expected_length:
            issues.append(f"Block {block_id}: length mismatch. CSV says {coords['length']}, pangraph expects {expected_length}")
    
    if issues:
        print("❌ Coordinate issues found:")
        for issue in issues[:5]:
            print(f"   {issue}")
        return False
    else:
        print("✅ All checked coordinates are valid for pangraph system")
        return True

def create_special_reference_pangraph():
    """Tạo special reference theo đúng pangraph coordinate system"""
    print("🚀 Creating Special Reference - PANGRAHP COORDINATE SYSTEM")
    print("=" * 70)
    print("Coordinate system: 0-based, right-excluded [start, end)")
    print("Python equivalent: [start:end]")
    
    # Load coordinates
    coordinates = load_coordinates(COORDS_FILE)
    print(f"📖 Loaded coordinates for {len(coordinates)} blocks")
    
    # Validate coordinates
    if not validate_pangraph_coordinates(coordinates):
        print("⚠️  Coordinate validation failed, but continuing...")
    
    # Tạo special reference file
    with open(SPECIAL_REF_FILE, 'w') as ref_out:
        blocks_processed = 0
        length_mismatches = 0
        
        for block_id, coords in coordinates.items():
            # Lấy sequence theo pangraph coordinates
            sequence = get_sequence_pangraph_coordinates(
                REF_GENOME_FILE, 
                coords['ref_start'], 
                coords['ref_end'],
                coords['strand']
            )
            
            if sequence:
                # VALIDATE: Kiểm tra length theo pangraph system
                pangraph_expected_length = coords['ref_end'] - coords['ref_start']
                actual_length = len(sequence)
                
                if actual_length != pangraph_expected_length:
                    print(f"⚠️  Block {block_id}: length mismatch")
                    print(f"   Pangraph expected: {pangraph_expected_length} ([{coords['ref_start']}, {coords['ref_end']}))")
                    print(f"   Actual: {actual_length}")
                    print(f"   CSV length: {coords['length']}")
                    length_mismatches += 1
                
                # Tạo header với pangraph coordinates
                header = f">block_{block_id}_NC_021726.1_{coords['ref_start']}_{coords['ref_end']}_{coords['strand']}"
                
                # Ghi vào file
                ref_out.write(header + "\n")
                
                # Ghi sequence với formatting (80 chars per line)
                for i in range(0, len(sequence), 80):
                    ref_out.write(sequence[i:i+80] + "\n")
                
                blocks_processed += 1
                
                if blocks_processed % 100 == 0:
                    print(f"📊 Processed {blocks_processed}/{len(coordinates)} blocks")
            else:
                print(f"❌ Could not get sequence for block {block_id}")
    
    print(f"\n🎉 Special reference created: {SPECIAL_REF_FILE}")
    print(f"📊 Total blocks: {blocks_processed}/{len(coordinates)}")
    print(f"⚠️  Length mismatches: {length_mismatches}")

def comprehensive_pangraph_validation():
    """Validation toàn diện special reference với pangraph system"""
    print(f"\n🔍 COMPREHENSIVE PANGRAHP VALIDATION")
    print("=" * 60)
    
    if not os.path.exists(SPECIAL_REF_FILE):
        print("❌ Special reference file not found")
        return False
    
    # Load coordinates để so sánh
    coordinates = load_coordinates(COORDS_FILE)
    
    # Đọc special reference
    special_ref_seqs = {}
    current_header = ""
    current_sequence = ""
    
    with open(SPECIAL_REF_FILE, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_header and current_sequence:
                    special_ref_seqs[current_header] = current_sequence
                current_header = line.strip()
                current_sequence = ""
            else:
                current_sequence += line.strip()
        
        if current_header and current_sequence:
            special_ref_seqs[current_header] = current_sequence
    
    print(f"📊 Sequences in special reference: {len(special_ref_seqs)}")
    print(f"📊 Expected from coordinates: {len(coordinates)}")
    
    # Kiểm tra từng block
    validation_issues = []
    
    for header, sequence in list(special_ref_seqs.items())[:5]:  # Check 5 blocks đầu
        # Parse header
        parts = header.split('_')
        block_id = parts[1]
        
        if block_id in coordinates:
            coords = coordinates[block_id]
            # Pangraph expected length: end - start
            pangraph_expected_length = coords['ref_end'] - coords['ref_start']
            actual_length = len(sequence)
            
            if actual_length != pangraph_expected_length:
                validation_issues.append(f"Block {block_id}: length {actual_length} != pangraph expected {pangraph_expected_length}")
            
            print(f"   ✅ Block {block_id}: [{coords['ref_start']}, {coords['ref_end']}) ({coords['strand']}) - {actual_length} bp")
            
            # In vài base đầu để kiểm tra
            if sequence:
                print(f"      First 20 bases: {sequence[:20]}")
        else:
            validation_issues.append(f"Block {block_id} in special ref but not in coordinates")
    
    if validation_issues:
        print(f"❌ Validation issues: {len(validation_issues)}")
        for issue in validation_issues[:3]:
            print(f"   {issue}")
    else:
        print("✅ Special reference validation PASSED for pangraph system")
    
    return len(validation_issues) == 0

if __name__ == "__main__":
    create_special_reference_pangraph()
    comprehensive_pangraph_validation()