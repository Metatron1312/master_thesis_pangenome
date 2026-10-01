#!/usr/bin/env python3
"""
MAP READS TO GRAPH - STANDALONE SCRIPT
=======================================
Script để map paired-end reads vào pangenome graph sử dụng vg giraffe
Đã hardcode đường dẫn cho SAMD00178099
"""

import os
import sys
import subprocess
import json
from pathlib import Path
from datetime import datetime

# ========== CẤU HÌNH ĐƯỜNG DẪN ==========
# Đường dẫn đến graph files
GRAPH_CONFIG = {
    "graph_gbz": "/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/vg_graph/graph/106_strains_with_refseq_complete_vg_autoindex/complete_106_strains_with_refseq_pangenome.giraffe.gbz",
    "min_index": "/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/vg_graph/graph/106_strains_with_refseq_complete_vg_autoindex/complete_106_strains_with_refseq_pangenome.shortread.withzip.min",
    "dist_index": "/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/vg_graph/graph/106_strains_with_refseq_complete_vg_autoindex/complete_106_strains_with_refseq_pangenome.dist",
    "threads": 4
}

# ========== ĐƯỜNG DẪN READS (HARDCODED) ==========
READS_CONFIG = {
    "sample_id": "SAMD00178099",
    "r1_file": "/Users/shinra/gene_detection_v2/validate_database/SAMD00178099/fastp_trimmed/SAMD00178099_1_trimmed.fastq.gz",
    "r2_file": "/Users/shinra/gene_detection_v2/validate_database/SAMD00178099/fastp_trimmed/SAMD00178099_2_trimmed.fastq.gz",
    "output_dir": "/Users/shinra/gene_detection_v2/validate_database/SAMD00178099/reads_mapped_to_graph"
}

def check_graph_files():
    """Kiểm tra các file graph có tồn tại không"""
    print("\n🔍 Checking graph files...")
    missing = []
    for name, path in GRAPH_CONFIG.items():
        if name == 'threads':
            continue
        if not os.path.exists(path):
            missing.append(f"  ❌ {name}: {path}")
        else:
            print(f"  ✅ {name}: {os.path.basename(path)}")
    
    if missing:
        print("\n❌ MISSING GRAPH FILES:")
        print("\n".join(missing))
        return False
    
    print("\n✅ All graph files are ready")
    return True

def check_reads_files():
    """Kiểm tra các file reads có tồn tại không"""
    print("\n🔍 Checking reads files...")
    
    if not os.path.exists(READS_CONFIG["r1_file"]):
        print(f"  ❌ R1 file not found: {READS_CONFIG['r1_file']}")
        return False
    else:
        file_size = os.path.getsize(READS_CONFIG["r1_file"]) / (1024**2)  # MB
        print(f"  ✅ R1: {os.path.basename(READS_CONFIG['r1_file'])} ({file_size:.1f} MB)")
    
    if not os.path.exists(READS_CONFIG["r2_file"]):
        print(f"  ❌ R2 file not found: {READS_CONFIG['r2_file']}")
        return False
    else:
        file_size = os.path.getsize(READS_CONFIG["r2_file"]) / (1024**2)  # MB
        print(f"  ✅ R2: {os.path.basename(READS_CONFIG['r2_file'])} ({file_size:.1f} MB)")
    
    return True

def map_reads_to_graph(output_gam_name=None):
    """
    Map paired-end reads vào graph sử dụng vg giraffe
    
    Args:
        output_gam_name: Tên file GAM output (không bao gồm extension)
                         Nếu None, sẽ dùng sample_id làm tên
    
    Returns:
        dict: Mapping results with status and alignment count
    """
    
    # Kiểm tra files
    if not check_reads_files():
        return {'status': 'error', 'error': 'Reads files missing'}
    
    if not check_graph_files():
        return {'status': 'error', 'error': 'Graph files missing'}
    
    # Tạo output directory
    output_dir = Path(READS_CONFIG["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Tạo tên output file
    if output_gam_name:
        gam_file = output_dir / f"{output_gam_name}.gam"
    else:
        gam_file = output_dir / f"{READS_CONFIG['sample_id']}.gam"
    
    print("\n" + "="*80)
    print(f"MAPPING READS TO GRAPH")
    print("="*80)
    print(f"Sample: {READS_CONFIG['sample_id']}")
    print(f"R1: {READS_CONFIG['r1_file']}")
    print(f"R2: {READS_CONFIG['r2_file']}")
    print(f"Output: {gam_file}")
    print("="*80)
    
    # Kiểm tra nếu output đã tồn tại
    if os.path.exists(gam_file) and os.path.getsize(gam_file) > 0:
        # Đếm số alignments hiện có
        count_cmd = f"vg view -a '{gam_file}' -j 2>/dev/null | wc -l"
        count_result = subprocess.run(count_cmd, shell=True, capture_output=True, text=True)
        try:
            alignment_count = int(count_result.stdout.strip())
        except:
            alignment_count = 0
        
        print(f"\n⏭️  Output already exists: {alignment_count} alignments")
        print(f"   To re-run, delete: {gam_file}")
        
        return {
            'status': 'skipped',
            'alignments': alignment_count,
            'gam_file': str(gam_file)
        }
    
    # Chuẩn bị lệnh vg giraffe
    giraffe_cmd = [
        "vg", "giraffe",
        "-Z", GRAPH_CONFIG["graph_gbz"],
        "-m", GRAPH_CONFIG["min_index"],
        "-d", GRAPH_CONFIG["dist_index"],
        "-f", READS_CONFIG["r1_file"],
        "-f", READS_CONFIG["r2_file"],
        "-t", str(GRAPH_CONFIG["threads"]),
        "-o", "gam"
    ]
    
    print(f"\n🚀 Running vg giraffe...")
    print(f"   Command: {' '.join(giraffe_cmd[:5])} ... -t {GRAPH_CONFIG['threads']}")
    
    try:
        # Chạy vg giraffe
        start_time = datetime.now()
        with open(gam_file, 'w') as out_f:
            result = subprocess.run(giraffe_cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
        
        elapsed_time = (datetime.now() - start_time).total_seconds()
        
        if result.returncode == 0:
            # Đếm số alignments
            count_cmd = f"vg view -a '{gam_file}' -j 2>/dev/null | wc -l"
            count_result = subprocess.run(count_cmd, shell=True, capture_output=True, text=True)
            try:
                alignment_count = int(count_result.stdout.strip())
            except:
                alignment_count = 0
            
            # Lấy kích thước file
            file_size = os.path.getsize(gam_file) / (1024**2)  # MB
            
            print(f"\n✅ MAPPING SUCCESSFUL!")
            print(f"   • Alignments: {alignment_count:,}")
            print(f"   • Output size: {file_size:.2f} MB")
            print(f"   • Time: {elapsed_time:.1f} seconds")
            print(f"   • File: {gam_file}")
            
            return {
                'status': 'success',
                'alignments': alignment_count,
                'gam_file': str(gam_file),
                'file_size_mb': file_size,
                'time_seconds': elapsed_time,
                'returncode': result.returncode
            }
        else:
            # Xóa file output nếu lỗi
            if os.path.exists(gam_file):
                os.remove(gam_file)
            
            error_msg = result.stderr[:1000] if result.stderr else "Unknown error"
            print(f"\n❌ MAPPING FAILED!")
            print(f"   Error: {error_msg[:500]}...")
            
            return {
                'status': 'failed',
                'alignments': 0,
                'error': error_msg,
                'returncode': result.returncode
            }
            
    except Exception as e:
        print(f"\n❌ ERROR: {str(e)}")
        return {
            'status': 'error',
            'alignments': 0,
            'error': str(e)
        }

def extract_gam_info(gam_file):
    """
    Trích xuất thông tin từ GAM file
    
    Args:
        gam_file: Path to GAM file
    
    Returns:
        dict: Information about alignments
    """
    if not os.path.exists(gam_file):
        return {'error': 'GAM file not found'}
    
    print(f"\n📊 Analyzing GAM file: {gam_file}")
    
    # Đếm total alignments
    count_cmd = f"vg view -a '{gam_file}' -j 2>/dev/null | wc -l"
    count_result = subprocess.run(count_cmd, shell=True, capture_output=True, text=True)
    try:
        total_alignments = int(count_result.stdout.strip())
    except:
        total_alignments = 0
    
    # Lấy thông tin về mapping quality từ 100 alignments đầu tiên
    mq_scores = []
    view_cmd = f"vg view -a '{gam_file}' -j 2>/dev/null | head -100"
    view_result = subprocess.run(view_cmd, shell=True, capture_output=True, text=True)
    
    for line in view_result.stdout.strip().split('\n'):
        if not line.strip():
            continue
        try:
            data = json.loads(line)
            if 'mapping_quality' in data:
                mq_scores.append(data['mapping_quality'])
        except:
            pass
    
    if mq_scores:
        avg_mq = sum(mq_scores) / len(mq_scores)
        min_mq = min(mq_scores)
        max_mq = max(mq_scores)
        print(f"   • Mapping quality (sample of 100): avg={avg_mq:.1f}, min={min_mq}, max={max_mq}")
    else:
        avg_mq = min_mq = max_mq = 0
        print(f"   • No mapping quality data available")
    
    info = {
        'total_alignments': total_alignments,
        'file_size_mb': os.path.getsize(gam_file) / (1024**2),
        'avg_mapping_quality': avg_mq,
        'min_mapping_quality': min_mq,
        'max_mapping_quality': max_mq
    }
    
    print(f"   • Total alignments: {total_alignments:,}")
    print(f"   • File size: {info['file_size_mb']:.2f} MB")
    
    return info

def main():
    print("\n" + "="*80)
    print("READS TO GRAPH MAPPING PIPELINE")
    print("="*80)
    print(f"Sample: {READS_CONFIG['sample_id']}")
    print(f"Output directory: {READS_CONFIG['output_dir']}")
    print("="*80)
    
    # Thực hiện mapping
    result = map_reads_to_graph()
    
    # Lưu báo cáo
    report_file = Path(READS_CONFIG['output_dir']) / f"mapping_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    
    with open(report_file, 'w') as f:
        json.dump(result, f, indent=2, default=str)
    
    print(f"\n📄 Report saved: {report_file}")
    
    # Nếu thành công, trích xuất thêm thông tin
    if result['status'] == 'success':
        gam_info = extract_gam_info(result['gam_file'])
        
        # Cập nhật report với thông tin chi tiết
        result['detailed_info'] = gam_info
        with open(report_file, 'w') as f:
            json.dump(result, f, indent=2, default=str)
    
    # In kết luận
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    
    if result['status'] == 'success':
        print(f"✅ MAPPING COMPLETED SUCCESSFULLY")
        print(f"   • Output: {result['gam_file']}")
        print(f"   • Alignments: {result['alignments']:,}")
        if 'detailed_info' in result:
            print(f"   • File size: {result['detailed_info']['file_size_mb']:.2f} MB")
    elif result['status'] == 'skipped':
        print(f"⏭️  MAPPING SKIPPED (output already exists)")
        print(f"   • Output: {result['gam_file']}")
        print(f"   • Alignments: {result['alignments']:,}")
    else:
        print(f"❌ MAPPING FAILED")
        print(f"   • Error: {result.get('error', 'Unknown')[:200]}")
    
    print("="*80)
    
    return 0 if result['status'] in ['success', 'skipped'] else 1

if __name__ == "__main__":
    sys.exit(main())