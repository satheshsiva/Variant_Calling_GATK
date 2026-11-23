#!/usr/bin/env bash
set -euo pipefail

REF="ref/chr20.fa"
PASS_VCF="Results/na12878_chr20.filtered.pass.vcf.gz"
NORM_VCF="Results/na12878_chr20.norm.vcf.gz"
CLINVAR="ref/clinvar.vcf.gz"
VEP_VCF="Results/na12878_chr20.vep_clinvar.vcf.gz"
SUMMARY_TSV="Results/candidate_clinical_variants.tsv"

echo "[1/3] Splitting multiallelic sites and left-normalizing indels..."
bcftools norm -m -any -f "$REF" "$PASS_VCF" -Oz -o "$NORM_VCF"
tabix -f -p vcf "$NORM_VCF"

echo "[2/3] Annotating normalized variants with Ensembl VEP (GRCh37 online DB) + local ClinVar..."
vep -i "$NORM_VCF" \
    -o "$VEP_VCF" \
    --vcf \
    --compress_output bgzip \
    --database \
    --assembly GRCh37 \
    --port 3337 \
    --symbol \
    --canonical \
    --hgvs \
    --check_existing \
    --sift b \
    --polyphen b \
    --custom "$CLINVAR",ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
    --force_overwrite

tabix -f -p vcf "$VEP_VCF"

echo "[3/3] Extracting candidate clinical/high-impact variants summary TSV..."
echo -e "CHROM\tPOS\tREF\tALT\tGENE\tCONSEQUENCE\tIMPACT\tHGVSp\tgnomADe_AF\tClinVar_CLNSIG\tClinVar_CLNDN" > "$SUMMARY_TSV"

bcftools +split-vep "$VEP_VCF" \
    -f "%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Consequence\t%IMPACT\t%HGVSp\t%gnomADe_AF\t%ClinVar_CLNSIG\t%ClinVar_CLNDN\n" \
    -d \
    -i 'IMPACT="HIGH" || ClinVar_CLNSIG~"Pathogenic" || ClinVar_CLNSIG~"Likely_pathogenic"' \
    >> "$SUMMARY_TSV"

echo "Step 4 complete! Annotated VCF: $VEP_VCF"
echo "Clinical summary table: $SUMMARY_TSV"
