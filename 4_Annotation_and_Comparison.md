# Step 4: Normalization, Functional & Clinical Annotation

**Scripts:**
* `scripts/04_clinvar_SnpEff_annotation.sh` *(Fast offline track: `bcftools annotate` + `SnpEff`)*
* `scripts/04_vep_clinvar_annotation.sh` *(Ensembl VEP track: offline cache or `--database`)*

---

## 1. Biological & Technical Rationale

### Multiallelic Splitting & Left-Normalization (`bcftools norm`)
Before attaching clinical or transcript annotations, the filtered `PASS` VCF (`Results/na12878_chr20.filtered.pass.vcf.gz`) must be standardized:
* **Left-Alignment (`-f ref/chr20.fa`):** Indels inside repetitive sequences can be represented at multiple equivalent coordinates. Left-normalization shifts every indel to its leftmost parsimonious position so allele-exact queries against ClinVar never miss variants due to coordinate shifting.
* **Multiallelic Splitting (`-m -any`):** Splits sites carrying multiple alternate alleles (`1,529` sites on contig `20`) into independent biallelic records so functional consequences and clinical classifications are never conflated across alleles.

```bash
bcftools norm -m -any -f ref/chr20.fa \
    Results/na12878_chr20.filtered.pass.vcf.gz \
    -Oz -o Results/na12878_chr20.norm.vcf.gz
tabix -f -p vcf Results/na12878_chr20.norm.vcf.gz
```
## 2. Two Annotation Tracks
Depending on available local storage, network latency, and downstream requirements, choose either Track A or Track B:  

|      Feature     |      Track A: ClinVar + SnpEff (04_clinvar_SnpEff_annotation.sh)     |           Track B: Ensembl VEP + ClinVar (04_vep_clinvar_annotation.sh)           |
|:----------------:|:--------------------------------------------------------------------:|:---------------------------------------------------------------------------------:|
| Execution Mode   | 100% local / offline (after initial GRCh37.75 DB fetch)              | Offline (--cache) or remote MySQL (--database --port 3337)                        |
| Runtime on Chr20 | ~30 seconds                                                          | ~2–3 minutes (--cache --fork 4) / ~1–2 hours (--database)                         |
| Annotation Tags  | INFO/CLNSIG, INFO/CLNDN, INFO/RS, INFO/ANN                           | INFO/CSQ (SYMBOL, Consequence, IMPACT, HGVSp, SIFT, PolyPhen, ClinVar)            |
| Best For         | Rapid local turnaround, offline environments, low network dependency | Canonical transcript filtering, SIFT/PolyPhen deleteriousness scores, VEP plugins |
---

### Track A: Local ClinVar + SnpEff (scripts/04_clinvar_SnpEff_annotation.sh)
1. Attach Clinical Annotations (`bcftools annotate`)
Matches normalized alleles against a local GRCh37 ClinVar VCF (ref/clinvar.vcf.gz) and transfers clinical significance (CLNSIG), review status (CLNREVSTAT), disease names (CLNDN), gene info (GENEINFO), and dbSNP rsIDs (RS, ID):

```Bash
bcftools annotate \
    -a ref/clinvar.vcf.gz \
    -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/GENEINFO,INFO/RS,ID \
    Results/na12878_chr20.norm.vcf.gz \
    -Oz -o Results/na12878_chr20.clinvar.vcf.gz
bcftools index -f -t Results/na12878_chr20.clinvar.vcf.gz
```

2. Predict Transcript Consequences (`SnpEff`)
Annotates protein-coding consequences, Sequence Ontology terms, HIGH/MODERATE/LOW/MODIFIER impacts, and HGVS notation into INFO/ANN using the GRCh37.75 database:

```Bash
snpEff -Xmx4g -v GRCh37.75 Results/na12878_chr20.clinvar.vcf.gz \
    | bgzip > Results/na12878_chr20.clinvar.snpeff.vcf.gz
bcftools index -f -t Results/na12878_chr20.clinvar.snpeff.vcf.gz
```

3. Extract Candidate Clinical Variants
Splits `INFO/ANN` records and filters for HIGH-impact variants or ClinVar Pathogenic / Likely_pathogenic classifications:  

```Bash
bcftools +split-vep Results/na12878_chr20.clinvar.snpeff.vcf.gz \
    -a ANN \
    -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%Gene_Name\t%Annotation\t%Annotation_Impact\t%HGVS.p\t%INFO/CLNSIG\t%INFO/CLNDN\n" \
    -d \
    -i 'Annotation_Impact="HIGH" || INFO/CLNSIG~"Pathogenic" || INFO/CLNSIG~"Likely_pathogenic"' \
    >> Results/candidate_clinical_variants.tsv
```

### Track B: Ensembl VEP + Custom ClinVar (scripts/04_vep_clinvar_annotation.sh)
#### Option B1: Offline Local VEP Cache (Recommended for VEP)
```Bash
vep --input_file Results/na12878_chr20.norm.vcf.gz \
    --output_file Results/na12878_chr20.vep_clinvar.vcf.gz \
    --format vcf \
    --vcf \
    --compress_output bgzip \
    --cache \
    --offline \
    --dir_cache "${VEP_CACHE:-$HOME/.vep}" \
    --assembly GRCh37 \
    --fasta ref/chr20.fa \
    --symbol \
    --canonical \
    --hgvs \
    --af_gnomade \
    --sift b \
    --polyphen b \
    --custom ref/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
    --fork 4 \
    --force_overwrite
tabix -f -p vcf Results/na12878_chr20.vep_clinvar.vcf.gz
```

#### Option B2: Online Ensembl GRCh37 MySQL Database (Zero Local Cache Footprint)
(Note: Omit --fork and --af_gnomade when querying the public EBI MySQL server.)

```Bash
vep -i Results/na12878_chr20.norm.vcf.gz \
    -o Results/na12878_chr20.vep_clinvar.vcf.gz \
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
    --custom ref/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
    --force_overwrite
tabix -f -p vcf Results/na12878_chr20.vep_clinvar.vcf.gz
```

### Extract Candidate Clinical Variants (VEP CSQ Output)
```Bash
bcftools +split-vep Results/na12878_chr20.vep_clinvar.vcf.gz \
    -f "%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Consequence\t%IMPACT\t%HGVSp\t%ClinVar_CLNSIG\t%ClinVar_CLNDN\n" \
    -d \
    -i 'IMPACT="HIGH" || ClinVar_CLNSIG~"Pathogenic" || ClinVar_CLNSIG~"Likely_pathogenic"' \
    >> Results/candidate_clinical_variants.tsv
```

### Execution & Post-Filter QC
```bash
chmod +x 04_vep_clinvar_annotation.sh
bash 04_vep_clinvar_annotation.sh


chmod +x 04_clinvar_SnpEff_annotation.sh
bash 04_clinvar_SnpEff_annotation.sh
```

### Why `bcftools annotate` (ClinVar) + `SnpEff` Is the Default Track
While `scripts/04_vep_clinvar_annotation.sh` is retained for workflows requiring Ensembl VEP (`SIFT`, `PolyPhen`, canonical filtering), `scripts/04_clinvar_SnpEff_annotation.sh` is used as the default annotation pipeline for four practical reasons:
* **Deterministic Local Execution:** Runs 100% offline in ~30–45 seconds on Chromosome 20, avoiding multi-hour latency and remote MySQL socket timeouts common with online `vep --database` queries.
* **Lightweight Disk Footprint:** Uses the compact `SnpEff` (`GRCh37.75`) database rather than multi-gigabyte offline VEP cache archives, making it ideal for student laptops and WSL environments.
* **Decoupled Clinical Updates:** Separating `bcftools annotate` (`clinvar.vcf.gz`) from transcript prediction (`SnpEff`) lets clinical classifications (`CLNSIG`, `CLNDN`, `RS`) be refreshed in seconds without re-running functional annotation.
* **Clean Modularity:** Isolating VEP and `SnpEff` into separate scripts prevents redundant transcript tags in a single VCF and lets users choose the tool best suited to their compute setup.

__END__