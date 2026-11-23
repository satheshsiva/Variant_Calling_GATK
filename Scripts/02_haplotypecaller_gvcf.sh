#!/usr/bin/env bash
set -euo pipefail

REF="ref/chr20.fa"
RECAL_BAM="data/processed/na12878_chr20.recal.bam"
DBSNP="ref/dbsnp_138.b37.chr20.vcf.gz"
OUT_DIR="Results"
JAVA_OPTS="-Xmx6g"

mkdir -p "$OUT_DIR"

echo "[1/2] Running HaplotypeCaller in GVCF mode on contig 20..."
gatk --java-options "$JAVA_OPTS" HaplotypeCaller \
    -R "$REF" \
    -I "$RECAL_BAM" \
    -L 20 \
    -ERC GVCF \
    --native-pair-hmm-threads 4 \
    -O "$OUT_DIR/na12878_chr20.g.vcf.gz"

echo "[2/2] Running GenotypeGVCFs to produce raw genotyped VCF..."
gatk --java-options "$JAVA_OPTS" GenotypeGVCFs \
    -R "$REF" \
    -V "$OUT_DIR/na12878_chr20.g.vcf.gz" \
    -D "$DBSNP" \
    -L 20 \
    -O "$OUT_DIR/na12878_chr20.raw.vcf.gz"

echo "Step 2 complete! Raw genotyped VCF: $OUT_DIR/na12878_chr20.raw.vcf.gz"
