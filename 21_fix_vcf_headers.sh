#!/bin/bash
# fix_vcf_headers_corrected.sh
VCF_DIR="/Users/shinra/Bioinformatics/WGS_AMR/data_public/vg_graph_v2/msa2vcf_v6"
FIXED_DIR="${VCF_DIR}_fixed"

mkdir -p "$FIXED_DIR"

for vcf_file in "$VCF_DIR"/*.vcf; do
    if [ -f "$vcf_file" ]; then
        filename=$(basename "$vcf_file")
        block_id=$(echo "$filename" | sed 's/block_//' | sed 's/\.vcf//')
        output_file="$FIXED_DIR/$filename"
        
        echo "Fixing $filename with block ID: $block_id"
        
        # Extract sample names từ file gốc
        sample_names=$(grep "^#CHROM" "$vcf_file" | cut -f10- | tr '\t' ' ')
        
        # Nếu không có sample names, sử dụng mẫu mặc định
        if [ -z "$sample_names" ]; then
            sample_names="4627273680257107215 13258632512581240699 15002515446254948549 16581053303069541447 18406285861287570027"
        fi
        
        # Tạo VCF với header đầy đủ
        {
            echo "##fileformat=VCFv4.2"
            echo "##contig=<ID=${block_id}>"
            echo "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
            echo "##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL\">"
            echo "##INFO=<ID=START_INSERTION,Number=0,Type=Flag,Description=\"Insertion at the start of sequence\">"
            echo "##INFO=<ID=END_INSERTION,Number=0,Type=Flag,Description=\"Insertion at the end of sequence\">"
            
            # Dòng #CHROM với sample names
            echo -n "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT"
            for sample in $sample_names; do
                echo -n "	$sample"
            done
            echo ""  # New line
            
            # Variant data
            grep -v "^#" "$vcf_file"
        } > "$output_file"
        
        # Nén và index
        bgzip -f "$output_file"
        tabix -p vcf "${output_file}.gz"
        
        echo "Fixed and indexed: ${output_file}.gz"
    fi
done
