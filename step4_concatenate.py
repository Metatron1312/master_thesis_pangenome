import os
import glob
from collections import defaultdict
from Bio import SeqIO

# ============================================================
# CẤU HÌNH
# ============================================================

ALIGNED_DIR = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/iqtree/core_genome/core_block_aligned"
CORE_BLOCKS = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/iqtree/core_genome/core_blocks.txt"
OUTPUT_DIR = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/phylogenetic/iqtree/core_genome"
OUTPUT_FASTA = os.path.join(OUTPUT_DIR, "core_genome_aln.fa")

# ============================================================
# BƯỚC 1: ĐỌC DANH SÁCH CORE BLOCK
# ============================================================

print("=" * 60)
print("BƯỚC 1: ĐỌC DANH SÁCH CORE BLOCK")
print("=" * 60)

with open(CORE_BLOCKS) as f:
    core_blocks = [line.strip() for line in f if line.strip()]

print(f"Số core block: {len(core_blocks)}")

# ============================================================
# BƯỚC 2: NỐI CÁC ALIGNMENT LẠI
# ============================================================

print("\n" + "=" * 60)
print("BƯỚC 2: NỐI CÁC ALIGNMENT LẠI")
print("=" * 60)

# Dictionary: sample_id -> concatenated sequence
concatenated = defaultdict(str)

# Biến đếm
total_columns = 0
blocks_processed = 0
blocks_failed = 0

for i, block_id in enumerate(core_blocks):
    aligned_file = os.path.join(ALIGNED_DIR, f"block_{block_id}.aln.fa")
    
    if not os.path.exists(aligned_file):
        print(f"  CẢNH BÁO: Không tìm thấy {aligned_file}")
        blocks_failed += 1
        continue
    
    # Đọc alignment
    records = list(SeqIO.parse(aligned_file, "fasta"))
    
    if len(records) == 0:
        print(f"  CẢNH BÁO: File rỗng {aligned_file}")
        blocks_failed += 1
        continue
    
    # Kiểm tra tất cả sequence có cùng độ dài không
    lengths = [len(record.seq) for record in records]
    if len(set(lengths)) != 1:
        print(f"  CẢNH BÁO: Block {block_id} có độ dài KHÔNG ĐỀU: {set(lengths)}")
        blocks_failed += 1
        continue
    
    block_length = lengths[0]
    total_columns += block_length
    blocks_processed += 1
    
    # Thêm sequence vào concatenated
    for record in records:
        sample_id = record.id
        concatenated[sample_id] += str(record.seq)
    
    if (i + 1) % 50 == 0:
        print(f"  Đã xử lý {i + 1}/{len(core_blocks)} block")

print(f"\nHoàn tất nối alignment!")
print(f"  Số block xử lý: {blocks_processed}")
print(f"  Số block thất bại: {blocks_failed}")
print(f"  Tổng số cột: {total_columns:,}")
print(f"  Số sample: {len(concatenated)}")

# ============================================================
# BƯỚC 3: KIỂM TRA ĐỘ DÀI
# ============================================================

print("\n" + "=" * 60)
print("BƯỚC 3: KIỂM TRA ĐỘ DÀI")
print("=" * 60)

lengths = [len(seq) for seq in concatenated.values()]
if len(set(lengths)) != 1:
    print(f"⚠️  CẢNH BÁO: Các sequence có độ dài KHÁC NHAU!")
    print(f"  Độ dài: {set(lengths)}")
else:
    print(f"✅ Tất cả sequence có cùng độ dài: {lengths[0]:,} cột")

# ============================================================
# BƯỚC 4: LƯU FILE ALIGNMENT
# ============================================================

print("\n" + "=" * 60)
print("BƯỚC 4: LƯU FILE ALIGNMENT")
print("=" * 60)

with open(OUTPUT_FASTA, "w") as f:
    for sample_id, seq in concatenated.items():
        f.write(f">{sample_id}\n")
        # Chia dòng 80 ký tự
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")

print(f"Đã lưu: {OUTPUT_FASTA}")
print(f"Kích thước: {os.path.getsize(OUTPUT_FASTA) / 1024 / 1024:.2f} MB")
print(f"Số mẫu: {len(concatenated)}")
print(f"Số cột: {lengths[0] if len(set(lengths)) == 1 else 'KHÔNG ĐỀU'}")

# ============================================================
# BƯỚC 5: KIỂM TRA FILE OUTPUT
# ============================================================

print("\n" + "=" * 60)
print("BƯỚC 5: KIỂM TRA FILE OUTPUT")
print("=" * 60)

# Đếm số mẫu
with open(OUTPUT_FASTA) as f:
    n_samples = sum(1 for line in f if line.startswith('>'))

print(f"Số mẫu trong file: {n_samples}")

# Đếm số cột
with open(OUTPUT_FASTA) as f:
    f.readline()  # Bỏ header
    first_seq = f.readline().strip()
    # Đọc thêm để đếm đủ
    while len(first_seq) < 1000:
        line = f.readline().strip()
        if not line or line.startswith('>'):
            break
        first_seq += line

print(f"100 ký tự đầu của sequence đầu tiên: {first_seq[:100]}")

# ============================================================
# TÓM TẮT
# ============================================================

print("\n" + "=" * 60)
print("TÓM TẮT")
print("=" * 60)
print(f"File output: {OUTPUT_FASTA}")
print(f"Số mẫu: {n_samples}")
print(f"Số cột: {lengths[0] if len(set(lengths)) == 1 else 'KHÔNG ĐỀU'}")
print(f"Kích thước: {os.path.getsize(OUTPUT_FASTA) / 1024 / 1024:.2f} MB")
