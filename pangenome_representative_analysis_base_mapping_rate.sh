#!/bin/bash
# ============================================================
# map_full_106_with_identity_filter.sh
# Thêm bộ lọc identity (mặc định 70%)
# ============================================================

set -e

INPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/combined_contigs_and_refseq_v2"
LINEAR_REF="/Users/shinra/Bioinformatics/WGS_AMR/data_public/combined_contigs_and_refseq/BJAB07104.filtered.complete.fasta"
BLOCK_FASTA="/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/graph_and_stats/pangenome_106_strains_with_refseq_block_cons.fa"

OUTPUT_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/pangraph_106_samples_with_refseq/Representative_analysis/base_mapping_rate_v2"
THREADS=8

# --- BỘ LỌC ---
MIN_MAPQ=10            # MAPQ tối thiểu
MIN_IDENTITY=70        # Identity tối thiểu (%)
MIN_ALN_LEN=100        # Độ dài alignment tối thiểu (bp) — tùy chọn

EXCLUDE_SAMPLE="BJAB07104"

mkdir -p "$OUTPUT_DIR"
CSV="$OUTPUT_DIR/base_mapping_rate_identity70.csv"

echo "biosample_id,ref_type,total_bases,aligned_bases,base_mapping_rate(%),n_alignments,mean_identity(%),median_identity(%)" > "$CSV"

# Index
samtools faidx "$LINEAR_REF"
samtools faidx "$BLOCK_FASTA"

# Hàm merge với 3 bộ lọc
merge_bp_filtered() {
    local PAF=$1
    awk -v q=$MIN_MAPQ -v id=$MIN_IDENTITY -v L=$MIN_ALN_LEN '
        NF >= 12 && $12 >= q && ($10 / $11 * 100) >= id && $11 >= L {
            print $1"\t"$3"\t"$4
        }
    ' "$PAF" | sort -k1,1 -k2,2n | awk '
        { if ($1 == pc && $2 <= pe) { if ($3 > pe) pe = $3 }
          else { if (pc != "") sum += pe - ps; pc=$1; ps=$2; pe=$3 } }
        END { if (pc != "") sum += pe - ps; print sum + 0 }
    '
}

# Hàm tính mean/median identity
identity_stats() {
    local PAF=$1
    awk -v q=$MIN_MAPQ -v id=$MIN_IDENTITY -v L=$MIN_ALN_LEN '
        NF >= 12 && $12 >= q && ($10 / $11 * 100) >= id && $11 >= L {
            ids[NR] = $10 / $11 * 100
            sum += $10 / $11 * 100
            n++
        }
        END {
            if (n == 0) { print "0\t0"; exit }
            mean = sum / n
            asort(ids)
            if (n % 2 == 1) median = ids[(n+1)/2]
            else median = (ids[n/2] + ids[n/2+1]) / 2
            printf "%.2f\t%.2f", mean, median
        }
    ' "$PAF"
}

# Đếm mẫu
N_TOTAL=$(ls "$INPUT_DIR"/*.filtered.complete.fasta 2>/dev/null | \
    grep -v "/${EXCLUDE_SAMPLE}.filtered.complete.fasta" | wc -l | tr -d ' ')
echo "Tìm thấy $N_TOTAL mẫu (đã loại $EXCLUDE_SAMPLE)"
echo "Bộ lọc: MAPQ≥$MIN_MAPQ, identity≥$MIN_IDENTITY%, length≥$MIN_ALN_LEN bp"
echo ""

i=0
for GENOME in "$INPUT_DIR"/*.filtered.complete.fasta; do
    SAMPLE=$(basename "$GENOME" .filtered.complete.fasta)
    [[ "$SAMPLE" == "$EXCLUDE_SAMPLE" ]] && continue

    i=$((i+1))
    echo "[$i/$N_TOTAL] $SAMPLE"

    samtools faidx "$GENOME" 2>/dev/null || true
    TOTAL_BASES=$(awk '{s+=$2} END{print s}' "${GENOME}.fai")

    for REF_TYPE in linear pangenome; do
        [[ "$REF_TYPE" == "linear" ]] && REF="$LINEAR_REF" || REF="$BLOCK_FASTA"

        PAF="/tmp/${SAMPLE}.${REF_TYPE}.$$.paf"
        minimap2 -x asm5 --no-long-join -t $THREADS "$REF" "$GENOME" 2>/dev/null > "$PAF"

        ALIGNED=$(merge_bp_filtered "$PAF")
        RATE=$(echo "scale=4; $ALIGNED * 100 / $TOTAL_BASES" | bc)
        N_ALN=$(awk -v q=$MIN_MAPQ -v id=$MIN_IDENTITY -v L=$MIN_ALN_LEN \
                'NF>=12 && $12>=q && ($10/$11*100)>=id && $11>=L {n++} END{print n+0}' "$PAF")
        read -r MEAN_ID MED_ID <<< $(identity_stats "$PAF")

        echo "$SAMPLE,$REF_TYPE,$TOTAL_BASES,$ALIGNED,$RATE,$N_ALN,$MEAN_ID,$MED_ID" >> "$CSV"
        echo "   [$REF_TYPE] rate=${RATE}% | aln=${N_ALN} | id_mean=${MEAN_ID}% | id_med=${MED_ID}%"

        rm -f "$PAF"
    done
done

# Dọn .fai
for GENOME in "$INPUT_DIR"/*.filtered.complete.fasta; do
    rm -f "${GENOME}.fai"
done

echo ""
echo "✅ HOÀN THÀNH"
echo "CSV: $CSV"
echo "Số dòng: $(( $(wc -l < "$CSV") - 1 ))"