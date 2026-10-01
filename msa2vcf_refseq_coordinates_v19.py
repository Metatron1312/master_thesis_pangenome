#!/usr/bin/env python3
import os
import glob
import re
import json
from Bio import SeqIO
from collections import defaultdict
import sys

# ==================== CONFIGURATION ====================
INPUT_MSA_DIR = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/vg_graph/aligned_blocks_with_coords"
SPECIAL_REF_FILE = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/vg_graph/special_reference_v3.fa"
OUTPUT_VCF_DIR = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/vg_graph/msa2vcf_v19"
# =======================================================

def reverse_complement(seq):
    """Tạo reverse complement của sequence"""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N', '-': '-',
            'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    return ''.join(comp.get(base, base) for base in reversed(seq))

def load_special_reference_sequences():
    """Load sequences và headers đầy đủ từ special reference"""
    ref_sequences = {}
    ref_headers = {}  # Lưu trữ header đầy đủ cho mỗi block
    try:
        for record in SeqIO.parse(SPECIAL_REF_FILE, "fasta"):
            full_header = record.id
            # Extract block_id từ header (phần thứ 2 sau "block_")
            header_parts = full_header.split('_')
            if len(header_parts) >= 2:
                block_id = header_parts[1]
                ref_sequences[block_id] = str(record.seq).upper()
                ref_headers[block_id] = full_header  # Lưu header đầy đủ
        print(f"✅ Loaded {len(ref_sequences)} sequences from special reference")
        return ref_sequences, ref_headers
    except Exception as e:
        print(f"❌ Error loading special reference: {e}")
        return {}, {}

def extract_block_id_from_filename(filename):
    """Extract block ID từ filename"""
    patterns = [
        r'block_(\d+)_aligned\.fasta',
        r'block_(\d+)_single\.fasta'
    ]
    
    for pattern in patterns:
        match = re.search(pattern, filename)
        if match:
            return match.group(1)
    return None

def find_all_msa_files():
    """Tìm tất cả file MSA"""
    aligned_files = glob.glob(os.path.join(INPUT_MSA_DIR, "block_*_aligned.fasta"))
    single_files = glob.glob(os.path.join(INPUT_MSA_DIR, "block_*_single.fasta"))
    all_files = aligned_files + single_files
    print(f"📁 Found {len(all_files)} MSA files: {len(aligned_files)} aligned + {len(single_files)} single")
    return all_files

def parse_msa_sequences(msa_file):
    """Parse tất cả sequences từ MSA"""
    sequences = []
    try:
        for record in SeqIO.parse(msa_file, "fasta"):
            sequences.append({
                'header': record.description,
                'sequence': str(record.seq).upper()
            })
        return sequences
    except Exception as e:
        print(f"Error parsing MSA file: {e}")
        return []

def extract_pangraph_coordinates(header):
    """Extract pangraph coordinates và strand từ header"""
    try:
        json_start = header.find('{')
        json_end = header.find('}') + 1
        if json_start != -1 and json_end != -1:
            json_str = header[json_start:json_end]
            data = json.loads(json_str)
            
            pangraph_start = data.get("start", 0)
            pangraph_end = data.get("end", 0)
            strand = data.get("strand", "+")
            
            return pangraph_start, pangraph_end, strand
    except Exception as e:
        pass
    
    return None, None, "+"

def extract_strain_name_consistent(header):
    """Extract strain name từ header"""
    strain_name = ""
    
    if "path_name" in header and "{" in header:
        try:
            json_str = header[header.find("{"):header.find("}")+1]
            data = json.loads(json_str)
            strain_name = data.get("path_name", "unknown_strain")
        except:
            pass
    
    if not strain_name and "path_name" in header:
        match = re.search(r'"path_name":"([^"]+)"', header)
        if match:
            strain_name = match.group(1)
    
    if not strain_name:
        clean_header = header.split()[0] if header else "unknown_strain"
        strain_name = clean_header
    
    return strain_name

def build_alignment_to_ref_map(ref_sequence_aligned):
    """Xây dựng mapping từ alignment position sang reference position (1-based)"""
    ref_positions = []
    ref_pos = 0
    
    for i, base in enumerate(ref_sequence_aligned):
        if base != '-':
            ref_pos += 1
            ref_positions.append((i, ref_pos))
        else:
            ref_positions.append((i, ref_pos))
    
    return ref_positions

def has_reference_sequence(sequences):
    """Kiểm tra xem MSA có chứa reference sequence không"""
    for seq in sequences:
        if "NC_021726.1" in seq['header']:
            return True
    return False

def create_empty_vcf(vcf_file, contig_name, sample_names, ref_length):
    """Tạo VCF file trống (chỉ có header)"""
    with open(vcf_file, 'w') as vcf_out:
        vcf_header = f"""##fileformat=VCFv4.2
##contig=<ID={contig_name},length={ref_length}>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t""" + "\t".join(sample_names) + "\n"
        vcf_out.write(vcf_header)
    return 0

def detect_variants_from_msa_corrected(sequences, ref_sequence_special, strand):
    """Phát hiện variants từ MSA - CORRECTED VERSION"""
    variants = defaultdict(lambda: {'snp': defaultdict(list), 'ins': defaultdict(list), 'del': defaultdict(list)})
    
    # Tìm reference sequence trong MSA
    ref_seq_msa = None
    for seq in sequences:
        if "NC_021726.1" in seq['header']:
            ref_seq_msa = seq['sequence']
            break
    
    if not ref_seq_msa:
        return variants, 0
    
    # Tính block length từ reference sequence (không tính gaps)
    block_length = len(ref_seq_msa.replace('-', ''))
    
    # Xây dựng alignment mapping
    ref_positions = build_alignment_to_ref_map(ref_seq_msa)
    
    # Tạo mapping từ alignment position đến reference position trong block
    align_to_ref_pos = {}
    for align_pos, ref_pos in ref_positions:
        align_to_ref_pos[align_pos] = ref_pos
    
    for seq in sequences:
        # Bỏ qua reference sequence
        if "NC_021726.1" in seq['header']:
            continue
            
        sample_name = extract_strain_name_consistent(seq['header'])
        sample_seq = seq['sequence']
        
        i = 0
        while i < len(ref_seq_msa) and i < len(sample_seq):
            ref_base_msa = ref_seq_msa[i]
            sample_base = sample_seq[i]
            
            # CASE 1: SNP
            if (ref_base_msa != '-' and sample_base != '-' and 
                ref_base_msa != sample_base and sample_base not in ['N', '.']):
                
                ref_pos_in_block = align_to_ref_pos.get(i, 0)
                if ref_pos_in_block > 0:
                    vcf_position = ref_pos_in_block
                    variants[vcf_position]['snp'][sample_base].append(sample_name)
                
                i += 1
            
            # CASE 2: INSERTION - CORRECTED: Xử lý insertion ở đầu đặc biệt
            elif ref_base_msa == '-' and sample_base != '-':
                insertion_start = i
                insertion_bases = ""
                
                # Thu thập insertion bases
                j = i
                while (j < len(ref_seq_msa) and j < len(sample_seq) and 
                       ref_seq_msa[j] == '-' and sample_seq[j] != '-'):
                    insertion_bases += sample_seq[j]
                    j += 1
                
                if insertion_bases:
                    # TÌM VỊ TRÍ CHÍNH XÁC CHO INSERTION
                    vcf_position = 1  # Mặc định cho insertion ở đầu
                    
                    # Tìm base reference đầu tiên
                    first_ref_base = None
                    first_ref_pos = 0
                    for k in range(len(ref_seq_msa)):
                        if ref_seq_msa[k] != '-':
                            first_ref_base = ref_seq_msa[k]
                            first_ref_pos = align_to_ref_pos.get(k, 1)
                            break
                    
                    # Nếu insertion ở trước base đầu tiên, dùng position 1
                    if first_ref_base and insertion_start < k:
                        vcf_position = 1
                    else:
                        # Insertion ở giữa, tìm base reference trước đó
                        for pos in range(insertion_start-1, -1, -1):
                            if ref_seq_msa[pos] != '-':
                                vcf_position = align_to_ref_pos.get(pos, 1)
                                break
                    
                    # Kiểm tra trùng lặp
                    existing_ins = False
                    for existing_bases in variants[vcf_position]['ins']:
                        if existing_bases == insertion_bases:
                            variants[vcf_position]['ins'][existing_bases].append(sample_name)
                            existing_ins = True
                            break
                    
                    if not existing_ins:
                        variants[vcf_position]['ins'][insertion_bases].append(sample_name)
                
                i = j
                continue
            
            # CASE 3: DELETION - CORRECTED: Đảm bảo phát hiện chính xác
            elif ref_base_msa != '-' and sample_base == '-':
                deletion_start = i
                deletion_bases = ""
                
                # Tìm reference position bắt đầu deletion
                deletion_ref_pos_in_block = align_to_ref_pos.get(i, 0)
                
                if deletion_ref_pos_in_block == 0:
                    i += 1
                    continue
                
                # Thu thập deletion bases từ REFERENCE trong MSA
                j = i
                while (j < len(ref_seq_msa) and j < len(sample_seq) and 
                       ref_seq_msa[j] != '-' and sample_seq[j] == '-'):
                    deletion_bases += ref_seq_msa[j]
                    j += 1
                
                if deletion_bases and deletion_ref_pos_in_block > 0:
                    vcf_position = deletion_ref_pos_in_block
                    
                    # Kiểm tra deletion không vượt quá reference length
                    if vcf_position <= len(ref_sequence_special):
                        # Kiểm tra trùng lặp
                        existing_del = False
                        for existing_bases in variants[vcf_position]['del']:
                            if existing_bases == deletion_bases:
                                variants[vcf_position]['del'][existing_bases].append(sample_name)
                                existing_del = True
                                break
                        
                        if not existing_del:
                            variants[vcf_position]['del'][deletion_bases].append(sample_name)
                
                i = j
                continue
            
            else:
                i += 1
    
    return variants, block_length

def write_variants_to_vcf_corrected(vcf_file, variants, contig_name, sample_names, ref_sequence_special, strand, block_length):
    """Ghi variants vào file VCF - CORRECTED VERSION"""
    written_variants = 0
    
    with open(vcf_file, 'w') as vcf_out:
        vcf_header = f"""##fileformat=VCFv4.2
##contig=<ID={contig_name},length={len(ref_sequence_special)}>
##INFO=<ID=STRAND,Number=1,Type=String,Description="Original strand">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t""" + "\t".join(sample_names) + "\n"
        vcf_out.write(vcf_header)
        
        sorted_positions = sorted(variants.keys())
        
        for vcf_pos in sorted_positions:
            if vcf_pos < 1 or vcf_pos > len(ref_sequence_special):
                continue
                
            pos_variants = variants[vcf_pos]
            
            # Xử lý SNPs
            for alt_base, samples in pos_variants['snp'].items():
                if not samples:
                    continue
                
                ref_base = ref_sequence_special[vcf_pos - 1] if vcf_pos - 1 < len(ref_sequence_special) else 'N'
                
                if ref_base != 'N' and ref_base != alt_base:
                    gt_fields = []
                    sample_found_count = 0
                    
                    for sample in sample_names:
                        if sample in samples:
                            gt_fields.append("1/1")
                            sample_found_count += 1
                        else:
                            gt_fields.append("0/0")
                    
                    if sample_found_count > 0:
                        info_field = f"STRAND={strand}"
                        vcf_line = f"{contig_name}\t{vcf_pos}\t.\t{ref_base}\t{alt_base}\t.\tPASS\t{info_field}\tGT\t" + "\t".join(gt_fields)
                        vcf_out.write(vcf_line + "\n")
                        written_variants += 1
            
            # Xử lý INSERTIONS - CORRECTED: Luôn có ALT base
            for ins_bases, samples in pos_variants['ins'].items():
                if not samples:
                    continue
                
                # Xác định REF base cho insertion
                if vcf_pos == 1:
                    # Insertion ở đầu: dùng base đầu tiên làm anchor
                    ref_base = ref_sequence_special[0] if len(ref_sequence_special) > 0 else 'N'
                else:
                    # Insertion ở giữa: dùng base tại vị trí insertion làm anchor
                    ref_base = ref_sequence_special[vcf_pos - 1] if vcf_pos - 1 < len(ref_sequence_special) else 'N'
                
                if ref_base != 'N':
                    alt_base = ref_base + ins_bases
                    
                    gt_fields = []
                    sample_found_count = 0
                    
                    for sample in sample_names:
                        if sample in samples:
                            gt_fields.append("1/1")
                            sample_found_count += 1
                        else:
                            gt_fields.append("0/0")
                    
                    if sample_found_count > 0:
                        info_field = f"INDEL;STRAND={strand}"
                        vcf_line = f"{contig_name}\t{vcf_pos}\t.\t{ref_base}\t{alt_base}\t.\tPASS\t{info_field}\tGT\t" + "\t".join(gt_fields)
                        vcf_out.write(vcf_line + "\n")
                        written_variants += 1
            
            # Xử lý DELETIONS - CORRECTED: LUÔN có ALT base
            for del_bases, samples in pos_variants['del'].items():
                if not samples:
                    continue
                
                deletion_length = len(del_bases)
                
                # Kiểm tra deletion có nằm trong reference
                if vcf_pos + deletion_length - 1 <= len(ref_sequence_special):
                    ref_bases = ref_sequence_special[vcf_pos-1:vcf_pos-1+deletion_length]
                    
                    # ALT base LUÔN là base đầu tiên của REF
                    alt_base = ref_sequence_special[vcf_pos - 1] if vcf_pos - 1 < len(ref_sequence_special) else 'N'
                    
                    # ĐẢM BẢO: ref_bases khớp với del_bases VÀ có ALT base
                    if (ref_bases != 'N' * deletion_length and 
                        ref_bases == del_bases and 
                        alt_base != 'N'):
                        
                        gt_fields = []
                        sample_found_count = 0
                        
                        for sample in sample_names:
                            if sample in samples:
                                gt_fields.append("1/1")
                                sample_found_count += 1
                            else:
                                gt_fields.append("0/0")
                        
                        if sample_found_count > 0:
                            info_field = f"INDEL;STRAND={strand}"
                            vcf_line = f"{contig_name}\t{vcf_pos}\t.\t{ref_bases}\t{alt_base}\t.\tPASS\t{info_field}\tGT\t" + "\t".join(gt_fields)
                            vcf_out.write(vcf_line + "\n")
                            written_variants += 1
    
    return written_variants

def msa_to_vcf_corrected(msa_file, ref_sequences, ref_headers, vcf_file):
    """Chuyển đổi file MSA thành VCF - CORRECTED VERSION với chrom name đầy đủ"""
    try:
        sequences = parse_msa_sequences(msa_file)
        if not sequences:
            print(f"  ⚠️  No sequences in {os.path.basename(msa_file)}")
            return -1
        
        if not has_reference_sequence(sequences):
            print(f"  ⚠️  No reference sequence in {os.path.basename(msa_file)} - SKIP")
            return -2
        
        block_id = extract_block_id_from_filename(os.path.basename(msa_file))
        if not block_id:
            print(f"  ❌ Cannot extract block_id from filename")
            return -1
        
        if block_id not in ref_sequences:
            print(f"  ❌ Block {block_id} not in special reference")
            return -1
        
        ref_sequence_special = ref_sequences[block_id]
        
        # SỬA QUAN TRỌNG: Sử dụng header đầy đủ từ special reference làm contig name
        contig_name = ref_headers[block_id]
        print(f"  🎯 Using contig name: {contig_name}")
        
        # Tìm reference sequence trong MSA
        ref_seq_msa = None
        pangraph_start = 0
        pangraph_end = 0
        strand = "+"
        
        for seq in sequences:
            if "NC_021726.1" in seq['header']:
                ref_seq_msa = seq['sequence']
                pangraph_start, pangraph_end, strand = extract_pangraph_coordinates(seq['header'])
                break
        
        if not ref_seq_msa:
            print(f"  ❌ No reference sequence in MSA for block {block_id}")
            return -1
        
        print(f"  📍 Pangraph coordinates: [{pangraph_start}, {pangraph_end}) strand: {strand}")
        
        # Tạo sample names (bỏ qua reference)
        sample_names = []
        for seq in sequences:
            if "NC_021726.1" not in seq['header']:
                strain_name = extract_strain_name_consistent(seq['header'])
                if strain_name and strain_name not in sample_names:
                    sample_names.append(strain_name)
        
        # Phát hiện variants
        variants, block_length = detect_variants_from_msa_corrected(
            sequences, ref_sequence_special, strand
        )
        
        # Thống kê variants
        total_snps = sum(len(samples) for pos_data in variants.values() for samples in pos_data['snp'].values())
        total_ins = sum(len(samples) for pos_data in variants.values() for samples in pos_data['ins'].values())
        total_del = sum(len(samples) for pos_data in variants.values() for samples in pos_data['del'].values())
        
        print(f"  🔍 Variants detected: {total_snps} SNPs, {total_ins} INS, {total_del} DEL")
        
        # Ghi VCF
        if variants:
            written_count = write_variants_to_vcf_corrected(
                vcf_file, variants, contig_name, sample_names, 
                ref_sequence_special, strand, block_length
            )
            print(f"  ✅ {written_count} variants written for block {block_id}")
            return written_count
        else:
            written_count = create_empty_vcf(vcf_file, contig_name, sample_names, len(ref_sequence_special))
            print(f"  ✅ Empty VCF created for block {block_id} (no variants)")
            return 0
        
    except Exception as e:
        print(f"  ❌ Error processing {os.path.basename(msa_file)}: {str(e)}")
        import traceback
        traceback.print_exc()
        return -1

def main():
    """Hàm main"""
    input_dir = INPUT_MSA_DIR
    output_dir = OUTPUT_VCF_DIR
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    print("🚀 MSA TO VCF - FINAL CORRECTED VERSION WITH FULL CHROM NAMES")
    print("=" * 60)
    print("CRITICAL FIXES:")
    print("✅ ALL variants have ALT bases")
    print("✅ Insertions at start properly handled") 
    print("✅ Deletions always have correct ALT bases")
    print("✅ VCF 4.2 compliant")
    print("✅ USING FULL CHROM NAMES FROM SPECIAL REFERENCE")  # QUAN TRỌNG
    
    # Load special reference - NHẬN CẢ sequences VÀ headers
    print("📖 Loading special reference sequences and headers...")
    ref_sequences, ref_headers = load_special_reference_sequences()
    if not ref_sequences:
        return
    
    # Tìm tất cả file MSA
    msa_files = find_all_msa_files()
    
    successful_with_variants = 0
    successful_empty = 0
    no_reference = 0
    failed = 0
    
    for i, msa_file in enumerate(msa_files, 1):
        block_id = extract_block_id_from_filename(os.path.basename(msa_file))
        
        if not block_id:
            print(f"  ❌ Cannot extract block_id from {os.path.basename(msa_file)}")
            failed += 1
            continue
        
        vcf_file = os.path.join(output_dir, f"block_{block_id}.vcf")
        
        print(f"\n[{i}/{len(msa_files)}] Processing {os.path.basename(msa_file)}...")
        
        # TRUYỀN THÊM ref_headers vào hàm
        result = msa_to_vcf_corrected(msa_file, ref_sequences, ref_headers, vcf_file)
        
        if result > 0:
            successful_with_variants += 1
        elif result == 0:
            successful_empty += 1
        elif result == -2:
            no_reference += 1
        else:
            failed += 1
        
        if i % 10 == 0:
            print(f"📊 Progress: {i}/{len(msa_files)} files processed")
    
    total_vcf_files = successful_with_variants + successful_empty
    print("=" * 60)
    print(f"🎉 FINAL RESULTS:")
    print(f"   Total MSA files: {len(msa_files)}")
    print(f"   VCF files created: {total_vcf_files}")
    print(f"     - With variants: {successful_with_variants}")
    print(f"     - Empty: {successful_empty}")
    print(f"   No reference (skipped): {no_reference}")
    print(f"   Failed: {failed}")

if __name__ == "__main__":
    main()