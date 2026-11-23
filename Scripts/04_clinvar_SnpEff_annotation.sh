#!/usr/bin/env bash
set -euo pipefail

REF="ref/chr20.fa"
PASS_VCF="Results/na12878_chr20.filtered.pass.vcf.gz"
NORM_VCF="Results/na12878_chr20.norm.vcf.gz"
CLINVAR="ref/clinvar.vcf.gz"

CLINVAR_VCF="Results/na12878_chr20.clinvar.vcf.gz"
SNPEFF_VCF="Results/na12878_chr20.clinvar.snpeff.vcf.gz"
SUMMARY_TSV="Results/candidate_clinical_variants.tsv"

echo "[1/3] Splitting multiallelic sites and left-normalizing indels..."
bcftools norm -m -any -f "$REF" "$PASS_VCF" -Oz -o "$NORM_VCF"
bcftools index -f -t "$NORM_VCF"

echo "[2/3] Running Two-Stage Offline Annotation (Local ClinVar + SnpEff)..."

echo "  -> Stage 2A: Attaching local ClinVar annotations via bcftools annotate..."
bcftools annotate \
    -a "$CLINVAR" \
    -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/GENEINFO,INFO/RS,ID \
    "$NORM_VCF" \
    -Oz -o "$CLINVAR_VCF"
bcftools index -f -t "$CLINVAR_VCF"

echo "  -> Stage 2B: Adding SnpEff transcript annotations (GRCh37.75)..."
snpEff -Xmx4g -v GRCh37.75 "$CLINVAR_VCF" | bgzip > "$SNPEFF_VCF"
bcftools index -f -t "$SNPEFF_VCF"

echo "[3/3] Extracting candidate clinical/high-impact variants summary TSV..."
echo -e "CHROM\tPOS\tID\tREF\tALT\tGENE\tCONSEQUENCE\tIMPACT\tHGVSp\tClinVar_CLNSIG\tClinVar_CLNDN" > "$SUMMARY_TSV"

bcftools query \
    -i 'INFO/ANN ~ "HIGH" || INFO/CLNSIG ~ "Pathogenic" || INFO/CLNSIG ~ "Likely_pathogenic"' \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNDN\t%INFO/ANN\n' \
    "$SNPEFF_VCF" \
| awk -F'\t' 'BEGIN {OFS="\t"} {
    n = split($8, anns, ",");
    for (i = 1; i <= n; i++) {
        split(anns[i], f, "|");
        # f[2]=Annotation, f[3]=Annotation_Impact, f[4]=Gene_Name, f[11]=HGVS.p
        if (f[3] == "HIGH" || $6 ~ /Pathogenic|Likely_pathogenic/) {
            hgvsp = (f[11] != "" ? f[11] : ".");
            clnsig = ($6 != "" ? $6 : ".");
            clndn = ($7 != "" ? $7 : ".");
            print $1, $2, $3, $4, $5, f[4], f[2], f[3], hgvsp, clnsig, clndn;
        }
    }
}' | sort -u >> "$SUMMARY_TSV"

echo "Step 4 complete! Final Annotated VCF: $SNPEFF_VCF"
echo "Clinical summary table: $SUMMARY_TSV"
