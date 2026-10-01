#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
collect_mapping_stats.py

Thu thập stats mapping cho mẫu SAMN48847114: graph-based (vg) vs linear (bwa).
Xuất ra 1 file CSV. Nếu mẫu đã có trong CSV, cập nhật dòng cũ.

Cách dùng:
    python3 collect_mapping_stats.py
"""

import os
import sys
import shutil
import subprocess
import pandas as pd
from datetime import datetime

# ==================== CẤU HÌNH ====================
SAMPLE_ID = "SAMN48847114"

# GAM file (graph-based mapping)
GRAPH_GAM = "/Users/shinra/gene_detection_v2/validate_database/SAMN48847114/reads_mapped_to_graph/SAMN48847114.gam"

# Graph .gbz — PHẢI cùng graph đã dùng để tạo GAM
GRAPH_GBZ = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/vg_graph/graph/106_strains_with_refseq_complete_vg_autoindex/complete_106_strains_with_refseq_pangenome.giraffe.gbz"

# Reference FASTA (linear) — dùng cho breadth coverage
REFERENCE = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/combined_contigs_and_refseq/BJAB07104.filtered.complete.fasta"

# BAM đã map (output của run_bwa_mapping.sh)
LINEAR_BAM = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/Mapping_and_variant_assessment/linear_reference/read_mapping_comparison/SAMN48847114/SAMN48847114.sorted.bam"

# Output CSV
OUTPUT_CSV = "/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/mapping_and_variant_calling/Mapping_and_variant_assessment/linear_reference/read_mapping_comparison/mapping_comparison.csv"

THREADS = 8
# ==================================================


def run_command(cmd, check=True):
    """Chạy command shell, trả về stdout."""
    try:
        result = subprocess.run(cmd, shell=True, check=check,
                                capture_output=True, text=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"❌ Command thất bại: {cmd}")
        print(f"   stderr: {e.stderr[:500] if e.stderr else 'N/A'}")
        if check:
            raise
        return ""


# ==================== GRAPH (GAM) ====================
def parse_gam_stats(gam_file):
    """Mapping rate, properly paired rate, MAPQ mean/median từ vg stats -a."""
    print(f"🔍 vg stats: {os.path.basename(gam_file)}")
    output = run_command(f"vg stats -a {gam_file}")

    s = {}
    for line in output.split("\n"):
        line = line.strip()
        try:
            if line.startswith("Total alignments:"):
                s["total_alignments"] = int(line.split(":")[1])
            elif line.startswith("Total aligned:"):
                s["total_aligned"] = int(line.split(":")[1])
            elif line.startswith("Total paired:"):
                s["total_paired"] = int(line.split(":")[1])
            elif line.startswith("Total properly paired:"):
                s["total_properly_paired"] = int(line.split(":")[1])
            elif line.startswith("Mapping quality: mean"):
                s["mapq_mean"]   = float(line.split("mean")[1].split(",")[0])
                s["mapq_median"] = float(line.split("median")[1].split(",")[0])
        except (ValueError, IndexError):
            continue

    s["mapping_rate"] = round(
        s.get("total_aligned", 0) / s["total_alignments"] * 100, 2
    ) if s.get("total_alignments", 0) > 0 else 0.0

    s["properly_paired_rate"] = round(
        s.get("total_properly_paired", 0) / s["total_paired"] * 100, 2
    ) if s.get("total_paired", 0) > 0 else 0.0

    return s


def parse_vg_depth(graph_gbz, gam_file):
    """
    vg depth -g → <mean> <stddev>
    Đo số alignment đi qua mỗi node được sample ngẫu nhiên của graph.
    KHÔNG so sánh trực tiếp với linear mean_depth.
    """
    print(f"🔍 vg depth: {os.path.basename(gam_file)}")
    cmd = f"vg depth {graph_gbz} -g {gam_file} -t {THREADS}"
    out = run_command(cmd)

    parts = out.strip().split()
    if len(parts) >= 2:
        try:
            mean = float(parts[0])
            std = float(parts[1])
            cv = round(std / mean, 4) if mean > 0 else None
            return {
                "mean_depth": round(mean, 2),
                "depth_stddev": round(std, 2),
                "depth_cv": cv,
            }
        except ValueError:
            pass
    print(f"   ⚠️ Không parse được vg depth output: {out[:120]}")
    return {"mean_depth": None, "depth_stddev": None, "depth_cv": None}


# ==================== LINEAR (BAM) ====================
def parse_flagstat(bam):
    """Mapping rate + properly paired rate từ samtools flagstat."""
    print(f"🔍 flagstat: {os.path.basename(bam)}")
    out = run_command(f"samtools flagstat {bam}")

    s = {}
    for line in out.split("\n"):
        if "in total" in line:
            s["total_reads"] = int(line.split()[0])
        elif "mapped (" in line and "primary" not in line:
            s["mapping_rate"] = float(line.split("(")[1].split("%")[0])
        elif "properly paired" in line and "(" in line:
            s["properly_paired_rate"] = float(line.split("(")[1].split("%")[0])
    return s


def parse_mapq(bam):
    """
    MAPQ mean + median từ samtools stats histogram.
    Tính cả mean và median từ histogram — không phụ thuộc version samtools.
    """
    print(f"🔍 samtools stats: {os.path.basename(bam)}")
    out = run_command(f"samtools stats {bam}")

    mapqs = []
    for line in out.split("\n"):
        # Format: "MAPQ\t<value>\t<count>"
        if line.startswith("MAPQ\t"):
            parts = line.split("\t")
            if len(parts) >= 3:
                try:
                    mapqs.append((int(parts[1]), int(parts[2])))
                except ValueError:
                    continue

    if not mapqs:
        return {"mapq_mean": 0.0, "mapq_median": 0.0}

    total = sum(c for _, c in mapqs)
    mean = sum(m * c for m, c in mapqs) / total

    # Median
    half = total / 2
    cum = 0
    median = 0
    for m, c in mapqs:
        cum += c
        if cum >= half:
            median = m
            break

    return {"mapq_mean": round(mean, 2), "mapq_median": float(median)}


def compute_coverage(bam):
    """Breadth of coverage (%) + mean depth trên reference phẳng."""
    print(f"🔍 coverage: {os.path.basename(bam)}")
    cmd = (
        f"samtools depth -a {bam} | "
        f"awk '{{sum+=$3; n++; if($3>0) cov++}} "
        f'END {{if(n>0) printf "%.2f\\t%.2f", cov*100/n, sum/n; else print "0\\t0"}}\''
    )
    out = run_command(cmd)
    try:
        breadth, depth = out.split("\t")
        return {
            "breadth_coverage_pct": float(breadth),
            "mean_depth": float(depth),
        }
    except (ValueError, IndexError):
        return {"breadth_coverage_pct": 0.0, "mean_depth": 0.0}


# ==================== MAIN ====================
def main():
    print("=== THU THẬP MAPPING STATS: GRAPH vs LINEAR ===")
    print(f"Thời gian: {datetime.now():%Y-%m-%d %H:%M:%S}")
    print(f"Mẫu: {SAMPLE_ID}\n")

    for tool in ["vg", "samtools"]:
        if shutil.which(tool) is None:
            print(f"❌ Thiếu '{tool}' trong PATH")
            sys.exit(1)

    row = {"sample_id": SAMPLE_ID}

    # ---- GRAPH ----
    if os.path.exists(GRAPH_GAM):
        g = parse_gam_stats(GRAPH_GAM)
        for k in ["mapping_rate", "properly_paired_rate", "mapq_mean", "mapq_median"]:
            row[f"graph_{k}"] = g.get(k)

        if os.path.exists(GRAPH_GBZ):
            d = parse_vg_depth(GRAPH_GBZ, GRAPH_GAM)
            row["graph_mean_depth"]   = d.get("mean_depth")
            row["graph_depth_stddev"] = d.get("depth_stddev")
            row["graph_depth_cv"]     = d.get("depth_cv")
        else:
            print(f"⚠️ Không tìm thấy GBZ: {GRAPH_GBZ}")
            row["graph_mean_depth"] = None
            row["graph_depth_stddev"] = None
            row["graph_depth_cv"] = None
    else:
        print(f"❌ Không tìm thấy GAM: {GRAPH_GAM}")

    # ---- LINEAR ----
    if os.path.exists(LINEAR_BAM):
        fs = parse_flagstat(LINEAR_BAM)
        mq = parse_mapq(LINEAR_BAM)
        cov = compute_coverage(LINEAR_BAM)

        row["linear_mapping_rate"]         = fs.get("mapping_rate")
        row["linear_properly_paired_rate"] = fs.get("properly_paired_rate")
        row["linear_mapq_mean"]            = mq.get("mapq_mean")
        row["linear_mapq_median"]          = mq.get("mapq_median")
        row["linear_breadth_coverage_pct"] = cov.get("breadth_coverage_pct")
        row["linear_mean_depth"]           = cov.get("mean_depth")
    else:
        print(f"❌ Không tìm thấy BAM: {LINEAR_BAM}")
        print(f"   → Chạy run_bwa_mapping.sh trước.")

    # ---- GHI / CẬP NHẬT CSV ----
    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)

    order = [
        "sample_id",
        "graph_mapping_rate", "linear_mapping_rate",
        "graph_properly_paired_rate", "linear_properly_paired_rate",
        "graph_mapq_mean", "linear_mapq_mean",
        "graph_mapq_median", "linear_mapq_median",
        "graph_mean_depth", "linear_mean_depth",
        "graph_depth_stddev", "graph_depth_cv",
        "linear_breadth_coverage_pct",
    ]

    df_new = pd.DataFrame([row])

    if os.path.exists(OUTPUT_CSV):
        df_old = pd.read_csv(OUTPUT_CSV)

        if SAMPLE_ID in df_old["sample_id"].values:
            df_old = df_old[df_old["sample_id"] != SAMPLE_ID]
            print(f"🔁 Mẫu {SAMPLE_ID} đã có trong CSV — cập nhật dòng cũ.")

        df = pd.concat([df_old, df_new], ignore_index=True)
    else:
        print(f"📄 Tạo CSV mới: {OUTPUT_CSV}")
        df = df_new

    # Sắp xếp cột
    order = [c for c in order if c in df.columns]
    other = [c for c in df.columns if c not in order]
    df = df[order + other]

    df.to_csv(OUTPUT_CSV, index=False)
    print(f"\n🎉 Đã lưu: {OUTPUT_CSV}")
    print(f"📊 Tổng số mẫu trong CSV: {len(df)}\n")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()