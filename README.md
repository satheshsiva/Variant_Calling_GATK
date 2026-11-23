# GATK Germline Short Variant Discovery & Clinical Annotation Pipeline (Chr20)

A modular, reproducible **GATK4-based germline variant discovery workflow inspired by GATK Best Practices.** workflow for human germline short variant calling (SNPs and Indels), hard filtering, and functional/clinical annotation on **NA12878 Chromosome 20 (GRCh37 / b37)**.

Designed as a lightweight, laptop-friendly companion to standard pileup variant calling tutorials, this repository separates executable Bash modules (`scripts/`) from step-by-step biological and algorithmic documentation.

---

## Repository Structure

```text
GATK_pipe/
├── README.md
├── scripts/
│   ├── 01_prep_and_bqsr.sh             # Reference indexing, MarkDuplicates, and BQSR
│   ├── 02_haplotypecaller_gvcf.sh      # Local haplotype assembly (GVCF mode) & joint genotyping
│   ├── 03_hard_filtering.sh            # SNP/Indel separation & GATK Hard Filtering
│   ├── 04_vep_clinvar_annotation.sh    # Multiallelic normalization, Ensembl VEP & ClinVar annotation
│   └── 04_clinvar_SnpEff_annotation.sh # Multiallelic normalization, ClinVar & SnpEff annotation
|
├── 0_Precheck.md
├── 1_Setup_and_Preprocessing.md
├── 2_Variant_Calling_GVCF.md
├── 3_GATK_Hard_Filtering.md
├── 4_Annotation_and_Comparison.md
├── 5_Clinical_Interpretation_and_Visualization.md
└── Results/
    └── candidate_clinical_variants_SnpEff.tsv # Tracked summary of high-impact / ClinVar candidate variants
```
## Software Requirements
All command-line tools are managed via Conda/Mamba (conda-forge and bioconda channels):  
    **GATK4** (gatk4-lite v4.7.0.0, Picard v3.5.0, HTSJDK v5.0.0)  
    **SAMtools** / BCFtools / HTSlib (samtools, bcftools, bgzip, tabix)  
    **Ensembl VEP** (ensembl-vep v116.x)
    **SnpEff** (SnpEff v5.4c)

```Bash
mamba create -n gatk_env -c conda-forge -c bioconda gatk4-lite samtools bcftools htslib snpeff ensembl-vep -y
conda activate gatk_env

snpEff download GRCh37.75
```
## Quick Start
Place input reference and alignment files in ref/ and data/ (excluded from Git via .gitignore):
    ref/chr20.fa and ref/chr20.fa.fai (GRCh37 contig 20)  
    ref/clinvar.vcf.gz and ref/clinvar.vcf.gz.tbi (GRCh37 ClinVar release)  
    data/na12878_chr20.bam and data/na12878_chr20.bam.bai
    
Execute the modular pipeline sequentially:
```Bash
bash scripts/01_prep_and_bqsr.sh
bash scripts/02_haplotypecaller_gvcf.sh
bash scripts/03_hard_filtering.sh
bash scripts/04_vep_clinvar_annotation.sh
bash scripts/04_clinvar_SnpEff_annotation.sh
```

Benchmark Summary (NA12878 Contig 20, GRCh37)

|               Pipeline Stage               | Total Records |  SNPs  | Indels | Ts/Tv Ratio |        dbSNP Concordance        |
|:------------------------------------------:|:-------------:|:------:|:------:|:-----------:|:-------------------------------:|
| Raw Genotyped VCF (GenotypeGVCFs)          | 110,997       | 91,874 | 19,227 | 2.06        | —                               |
| Hard-Filtered PASS VCF (VariantFiltration) | 105,417       | 86,504 | 19,017 | 2.10        | 98.09% (SNPs) / 92.41% (Indels) |

### Detailed Guides
Step 1: Setup, Duplicate Marking & BQSR  
Step 2: GVCF Variant Discovery & Joint Genotyping  
Step 3: GATK Hard Filtering & QC Metrics  
Step 4: Normalization, Ensembl VEP, ClinVar & SnpEff Annotation  
Step 5: Biological Interpretation & Alignment Visualization  

## Limitations

- Analysis is restricted to chromosome 20.
- Single-sample germline variant calling is performed.
- Hard-filtering thresholds are demonstration parameters
  and require dataset-specific validation.
- ClinVar classifications are database annotations,
  not independent clinical diagnoses.
- Functional predictions from SnpEff and VEP require
  transcript-aware interpretation.
- Results should not be used for clinical decision-making.

```bash
gatk --version
#The Genome Analysis Toolkit (GATK) v4.7.0.0
#HTSJDK Version: 5.0.0
#Picard Version: 3.5.0

samtools --version
#samtools 1.21
#Using htslib 1.21

bcftools --version
#bcftools 1.14
#Using htslib 1.21

npEff -version
#SnpEff  5.4c    2026-02-23

vep --help | head
#----------------------------------#
# ENSEMBL VARIANT EFFECT PREDICTOR #
#----------------------------------#
#Versions:
#  ensembl              : 116.0d85231
#  ensembl-compara      : 116.4038708
#  ensembl-funcgen      : 116.90049ea
#  ensembl-io           : 116.6afb5dc
#  ensembl-variation    : 116.db7394a
```

__END__