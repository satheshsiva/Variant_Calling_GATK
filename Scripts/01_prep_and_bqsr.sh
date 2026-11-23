#!/usr/bin/env bash
set -euo pipefail

# Directories & memory cap (safe for 16 GB RAM laptops)
REF="ref/chr20.fa"
RAW_BAM="data/na12878_chr20.bam"
OUT_DIR="data/processed"
JAVA_OPTS="-Xmx6g"

mkdir -p "$OUT_DIR"

echo "[1/5] Creating FASTA sequence dictionary..."
if [ ! -f "ref/chr20.dict" ]; then
    gatk --java-options "$JAVA_OPTS" CreateSequenceDictionary \
        -R "$REF" \
        -O "ref/chr20.dict"
fi

echo "[2/5] Fetching GRCh37 contig 20 known sites for BQSR..."
KNOWN_SNPS="ref/dbsnp_138.b37.chr20.vcf.gz"
KNOWN_INDELS="ref/Mills_1000G.b37.chr20.vcf.gz"

if [ ! -f "$KNOWN_SNPS" ]; then
    tabix -h https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/dbsnp_138.b37.vcf.gz 20 \
        | bgzip > "$KNOWN_SNPS" && tabix -p vcf "$KNOWN_SNPS"
fi

if [ ! -f "$KNOWN_INDELS" ]; then
    tabix -h https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Mills_and_1000G_gold_standard.indels.b37.vcf.gz 20 \
        | bgzip > "$KNOWN_INDELS" && tabix -p vcf "$KNOWN_INDELS"
fi

echo "[3/5] Marking PCR/optical duplicates..."
gatk --java-options "$JAVA_OPTS" MarkDuplicates \
    -I "$RAW_BAM" \
    -O "$OUT_DIR/na12878_chr20.dedup.bam" \
    -M "$OUT_DIR/na12878_chr20.dup_metrics.txt" \
    --CREATE_INDEX true

echo "[4/5] Building BQSR recalibration table..."
gatk --java-options "$JAVA_OPTS" BaseRecalibrator \
    -I "$OUT_DIR/na12878_chr20.dedup.bam" \
    -R "$REF" \
    --known-sites "$KNOWN_SNPS" \
    --known-sites "$KNOWN_INDELS" \
    -L 20 \
    -O "$OUT_DIR/recal_data.table"

echo "[5/5] Applying BQSR to generate analysis-ready BAM..."
gatk --java-options "$JAVA_OPTS" ApplyBQSR \
    -I "$OUT_DIR/na12878_chr20.dedup.bam" \
    -R "$REF" \
    --bqsr-recal-file "$OUT_DIR/recal_data.table" \
    -L 20 \
    -O "$OUT_DIR/na12878_chr20.recal.bam"

echo "Step 1 complete! Analysis-ready BAM: $OUT_DIR/na12878_chr20.recal.bam"
