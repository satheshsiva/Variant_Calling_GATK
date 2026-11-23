# Step 5: Biological Interpretation & Alignment Visualization (`samtools tview`)

**Input Files:**
* `Results/candidate_clinical_variants.tsv` (50 extracted `HIGH`-impact / ClinVar candidate rows)
* `Results/na12878_chr20.clinvar.snpeff.vcf.gz` (Annotated VCF with `GT`, `AD`, and `DP` format fields)
* `data/processed/na12878_chr20.recal.bam` (BQSR-recalibrated alignment BAM)
* `ref/chr20.fa` (GRCh37 contig `20` reference FASTA)

---

## 1. Overview of Extracted Clinical & High-Impact Variants

Querying `Results/candidate_clinical_variants.tsv` across Chromosome 20 (`63.0 Mb`) in `NA12878` yields **50 candidate transcript rows** carrying `Annotation_Impact="HIGH"` (`frameshift_variant`, `stop_gained`, `splice_donor_variant`, `splice_acceptor_variant`) and **0 ClinVar `Pathogenic` or `Likely_pathogenic` variants**:

```bash
wc -l Results/candidate_clinical_variants.tsv
grep -iE "Pathogenic|Likely_pathogenic" Results/candidate_clinical_variants.tsv
```
*Biological Takeaway:* NA12878 is a healthy reference individual, and Chromosome 20 represents ~2% of the human genome. Finding zero pathogenic variants confirms high specificity against false clinical positives while highlighting why computational HIGH-impact predictions must always be evaluated alongside zygosity (GT), population frequency, and ClinVar curation.

## 2. Case Studies in Variant Interpretation
Inspect sample-level genotypes (GT), allelic depths (AD), and total read depth (DP) for standout candidate loci:

```Bash
bcftools query -r 20:126155,20:1895949-1895950,20:25281333,20:25283938 \
    -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD\t%DP]\n' \
    Results/na12878_chr20.clinvar.snpeff.vcf.gz | column -t
```

| Locus       | Gene    | Variant (HGVSp)     | GT   | AD (REF,ALT) | DP | ClinVar Status          | Interpretation                                                                                                |
|-------------|---------|---------------------|------|--------------|----|-------------------------|---------------------------------------------------------------------------------------------------------------|
| 20:25281333 | ABHD12  | G>A (splice donor)  | 1/1  | 0,26         | 26 | Benign (PHARC syndrome) | Homozygous alternate in a healthy donor confirms benign status despite HIGH splice prediction.                |
| 20:25283938 | ABHD12  | A>AG (p.Ile194fs)   | 1/1  | 0,24         | 24 | Benign (PHARC syndrome) | Clean homozygous 1-bp insertion (24/24 reads) tolerated across population cohorts.                            |
| 20:1895949  | SIRPA   | A>AGT (p.Asp95fs)   | 1\|1 | 0,27         | 27 | Unlisted (.)            | Phased (1\|1) with adjacent deletion at 20:1895950; +2 bp and -2 bp cancel out to preserve the reading frame. |
| 20:1895950  | SIRPA   | CCT>C (p.Leu96fs)   | 1\|1 | 0,27         | 27 | Unlisted (.)            | Combined with 20:1895949, forms a net 0-bp complex substitution rather than a true frameshift.                |
| 20:126155   | DEFB126 | GCAAA>G (p.Gln55fs) | 0/1  | 22,7         | 29 | Benign                  | Common benign beta-defensin frameshift polymorphism showing mild heterozygous allelic imbalance (24% alt).    |
| 20:25755907 | FAM182B | C>A (p.Glu14*)      | 0/1  | —            | —  | Unlisted (.)            | Appears across three rows (p.Glu14*, p.Glu17*, p.Glu25*) due to overlapping transcript isoforms in SnpEff.    |


## 3. Terminal Alignment Verification with samtools tview
Before reporting any candidate variant, inspect the raw read pileup in data/processed/na12878_chr20.recal.bam to rule out mapping artifacts, strand bias, or split-indel phasing illusions:
samtools tview Interactive Shortcuts
. — Toggle between dot/comma match display and explicit base letters (A/C/G/T)  

g — Jump directly to a specific coordinate (enter 20:25283938 and press Enter)  

r — Toggle read names vs. base alignment view  

q — Quit `samtools tview`  

### A. Inspecting Phased Adjacent Indels (SIRPA, 20:1895949–1895950)
```Bash
samtools tview -p 20:1895930 data/processed/na12878_chr20.recal.bam ref/chr20.fa
```
What to Look For: Every read spanning 1895949–1895950 carries both the GT insertion and the CT deletion on the same chromosome copy, visually confirming why HaplotypeCaller assigned phased genotypes (1|1) and why the reading frame remains intact downstream.

### B. Inspecting a Clean Homozygous Frameshift Indel (ABHD12, 20:25283938)
```Bash
samtools tview -p 20:25283920 data/processed/na12878_chr20.recal.bam ref/chr20.fa
```
What to Look For: Uniform +G insertion across all 24 forward and reverse reads with high base qualities and zero flanking mismatch noise.

### C. Inspecting Heterozygous Allelic Balance (DEFB126, 20:126155)
```Bash
samtools tview -p 20:126140 data/processed/na12878_chr20.recal.bam ref/chr20.fa
```
What to Look For: A mix of reference-matching reads and 4-bp deletion reads, illustrating how local haplotype assembly recovers heterozygous indels even when soft-clipping reduces raw pileup alternate depth (AD = 22,7).

__END__