#!/usr/bin/env bash
set -euo pipefail

REF="ref/chr20.fa"
RAW_VCF="Results/na12878_chr20.raw.vcf.gz"
TMP_DIR="Results/filter_tmp"
OUT_VCF="Results/na12878_chr20.filtered.pass.vcf.gz"
JAVA_OPTS="-Xmx6g"

mkdir -p "$TMP_DIR"

echo "[1/5] Extracting raw SNPs..."
gatk --java-options "$JAVA_OPTS" SelectVariants \
    -R "$REF" \
    -V "$RAW_VCF" \
    --select-type-to-include SNP \
    -O "$TMP_DIR/raw_snps.vcf.gz"

echo "[2/5] Applying GATK Hard Filters to SNPs..."
gatk --java-options "$JAVA_OPTS" VariantFiltration \
    -R "$REF" \
    -V "$TMP_DIR/raw_snps.vcf.gz" \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O "$TMP_DIR/filtered_snps.vcf.gz"

echo "[3/5] Extracting raw Indels (and mixed variants)..."
gatk --java-options "$JAVA_OPTS" SelectVariants \
    -R "$REF" \
    -V "$RAW_VCF" \
    --select-type-to-include INDEL \
    --select-type-to-include MIXED \
    -O "$TMP_DIR/raw_indels.vcf.gz"

echo "[4/5] Applying GATK Hard Filters to Indels..."
gatk --java-options "$JAVA_OPTS" VariantFiltration \
    -R "$REF" \
    -V "$TMP_DIR/raw_indels.vcf.gz" \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O "$TMP_DIR/filtered_indels.vcf.gz"

echo "[5/5] Merging filtered SNPs + Indels and retaining PASS variants..."
gatk --java-options "$JAVA_OPTS" MergeVcfs \
    -I "$TMP_DIR/filtered_snps.vcf.gz" \
    -I "$TMP_DIR/filtered_indels.vcf.gz" \
    -O "$TMP_DIR/merged_filtered.vcf.gz"

gatk --java-options "$JAVA_OPTS" SelectVariants \
    -R "$REF" \
    -V "$TMP_DIR/merged_filtered.vcf.gz" \
    --exclude-filtered true \
    -O "$OUT_VCF"

echo "Step 3 complete! Clean PASS VCF: $OUT_VCF"
