## Workspace Hygiene and Environment Setup

1. Verify or Install gatk4 (via Conda/Mamba)  
Goal: Check if gatk is already available in your active Conda environment; if not, install gatk4-lite (lightweight GATK4 without heavy Spark dependencies) alongside samtools and bcftools.

```bash
# Check if gatk is already installed
which gatk && gatk --version

# If NOT found, create/update a clean environment:
# mamba create -n gatk_env -c conda-forge -c bioconda gatk4-lite samtools bcftools htslib snpeff ensembl-vep -y
# conda activate gatk_env
```
*Expected Output:* gatk --version prints The Genome Analysis Toolkit (GATK) v4.x.x.x (e.g., v4.6.2.0 or v4.7.0.0).

2. Download necessary inputs for the analysis  

```bash
# Stream ONLY chromosome 20 directly from the remote GIAB BAM
samtools view -b \
  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam \
  20 > na12878_chr20.bam

# Remove the cached remote .bai file that samtools downloads into your working directory
rm -f RMNISTHS_30xdownsample.bam.bai

# Index your new local chr20 BAM
samtools index na12878_chr20.bam

# 2. Download matching GRCh37 Chromosome 20 reference FASTA
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
samtools faidx hs37d5.fa.gz 20 > chr20.fa
rm -f hs37d5.fa.gz hs37d5.fa.gz.fai hs37d5.fa.gz.gzi
samtools faidx chr20.fa

# Code to directly download chr20 from web
# wget [https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.20.fa.gz](https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.20.fa.gz) && gunzip -c Homo_sapiens.GRCh37.dna.chromosome.20.fa.gz > chr20.fa && samtools faidx chr20.fa


# Download GRCh37 ClinVar VCF (~150 MB) and its tabix index
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi

# organize the files.
mkdir -p ref data
mv na12878_chr20.bam  na12878_chr20.bam.bai data/
mv chr20.fa  chr20.fa.fai  clinvar.vcf.gz  clinvar.vcf.gz.tbi ref/
```
  
3. Inspect BAM Header & Contig Naming (chr20 vs 20)  
*Goal:* Confirm whether your chr20.fa and na12878_chr20.bam use chr20 or 20 chromosome naming, and check if Read Groups (@RG) exist (GATK requires @RG tags in the BAM header).

```Bash
head -n 1 ref/chr20.fa
samtools view -H data/na12878_chr20.bam | grep -E "^@RG|SN:20" && tabix -l ref/clinvar.vcf.gz | grep -w "20\|chr20"
```
*Expected Output:* @SQ SN:20 LN:63025520  
One or more @RG lines (e.g., @RG ID:... SM:NA12878)  
20 from clinvar.vcf.gz (confirming GRCh37 contig alignment across FASTA, BAM, and ClinVar).

```
@SQ     SN:20   LN:63025520
@RG     ID:U0a_A2_L1    PL:illumina     PU:1    LB:1    SM:NA12878
@RG     ID:U0a_A2_L2    PL:illumina     PU:1    LB:1    SM:NA12878
@RG     ID:U0a_A_L1     PL:illumina     PU:1    LB:1    SM:NA12878 ...
```

__END__
