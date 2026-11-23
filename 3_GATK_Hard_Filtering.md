# Step 3: GATK Hard Filtering & Post-Filter QC

**Script:** `scripts/03_hard_filtering.sh`

## 1. Why Hard Filtering Instead of VQSR?
Variant Quality Score Recalibration (`VQSR`) fits Gaussian mixture models across multiple variant annotations and requires **whole-genome variant counts** to converge reliably. On single-chromosome datasets (`chr20`), GATK recommends **Hard Filtering** applied separately to SNPs and Indels.

---

## 2. GATK Hard-Filtering Annotations

### SNP Thresholds
| Annotation | Filter Rule | Biological / Technical Meaning |
| :--- | :--- | :--- |
| **`QD`** (QualByDepth) | `QD < 2.0` | Variant confidence (`QUAL`) normalized by depth; catches artifacts driven solely by ultra-high read depth. |
| **`QUAL`** | `QUAL < 30.0` | Minimum Phred-scaled variant site quality threshold. |
| **`FS`** (FisherStrand) | `FS > 60.0` | Fisher's exact test p-value detecting strand bias (alternate allele skewed toward forward or reverse reads). |
| **`SOR`** (StrandOddsRatio) | `SOR > 3.0` | Symmetric odds-ratio test for strand bias; complements `FS` at high coverage. |
| **`MQ`** (RMSMappingQuality) | `MQ < 40.0` | Root-mean-square mapping quality across spanning reads; removes calls in repetitive/ambiguous regions. |
| **`MQRankSum`** | `MQRankSum < -12.5` | Wilcoxon rank-sum test comparing mapping quality of `ALT` reads vs. `REF` reads. |
| **`ReadPosRankSum`** | `ReadPosRankSum < -8.0` | Detects whether `ALT` alleles cluster near noisy read ends rather than distributing randomly. |

### Indel Thresholds
| Annotation | Filter Rule | Biological / Technical Meaning |
| :--- | :--- | :--- |
| **`QD`** | `QD < 2.0` | Depth-normalized variant confidence. |
| **`QUAL`** | `QUAL < 30.0` | Minimum site quality score. |
| **`FS`** | `FS > 200.0` | Strand-bias threshold tailored for indel error profiles. |
| **`ReadPosRankSum`** | `ReadPosRankSum < -20.0` | End-of-read clustering filter for insertions and deletions. |

---

## 3. Execution & Post-Filter QC

```bash
chmod +x scripts/03_hard_filtering.sh
bash scripts/03_hard_filtering.sh
```

```
bcftools stats Results/na12878_chr20.filtered.pass.vcf.gz | grep -E "^SN|^TSTV"
```
```
SN      0       number of samples:      1
SN      0       number of records:      105417
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 86504
SN      0       number of MNPs: 0
SN      0       number of indels:       19017
SN      0       number of others:       0
SN      0       number of multiallelic sites:   1529
SN      0       number of multiallelic SNP sites:       25
TSTV    0       58608   27921   2.10    58504   27798   2.10
```
**Benchmark Output** (`na12878_chr20.filtered.pass.vcf.gz`)
Clean PASS Records: 105,417 (86,504 SNPs, 19,017 Indels; ~5.0% low-confidence calls removed)

**Ts/Tv Ratio**: Improved from 2.06 to 2.10

### dbSNP Concordance (`CollectVariantCallingMetrics`)
Ensure the single-chromosome dbSNP VCF header dictionary matches ref/chr20.dict before running Picard metrics:

```Bash
gatk UpdateVcfSequenceDictionary \
    -I ref/dbsnp_138.b37.chr20.vcf.gz \
    --SEQUENCE_DICTIONARY ref/chr20.dict \
    -O ref/dbsnp_138.b37.chr20.fixed.vcf.gz
mv ref/dbsnp_138.b37.chr20.fixed.vcf.gz ref/dbsnp_138.b37.chr20.vcf.gz
tabix -f -p vcf ref/dbsnp_138.b37.chr20.vcf.gz

gatk CollectVariantCallingMetrics \
    -I Results/na12878_chr20.filtered.pass.vcf.gz \
    --DBSNP ref/dbsnp_138.b37.chr20.vcf.gz \
    -O Results/na12878_chr20.qc_metrics
```

```
TOTAL_SNPS  NUM_IN_DB_SNP  NOVEL_SNPS  FILTERED_SNPS  PCT_DBSNP  DBSNP_TITV  NOVEL_TITV  TOTAL_INDELS  NOVEL_INDELS  FILTERED_INDELS  PCT_DBSNP_INDELS  NUM_IN_DB_SNP_INDELS  DBSNP_INS_DEL_RATIO  NOVEL_INS_DEL_RATIO  TOTAL_MULTIALLELIC_SNPS  NUM_IN_DB_SNP_MULTIALLELIC  TOTAL_COMPLEX_INDELS  NUM_IN_DB_SNP_COMPLEX_INDELS  SNP_REFERENCE_BIAS  NUM_SINGLETONS
86206       84563          1643        0              0.980941   2.124326    1.409091    17682         1342          0                0.924104          16340                 0.957824             1.096875             194                      160                         1163                  1118                          0.542202            72041
```
dbSNP SNP Concordance (PCT_DBSNP): 98.09% (84,563 / 86,206 biallelic SNPs; DBSNP_TITV = 2.12)

dbSNP Indel Concordance (PCT_DBSNP_INDELS): 92.41% (16,340 / 17,682)

__END__