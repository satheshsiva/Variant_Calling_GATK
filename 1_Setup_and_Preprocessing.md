# Step 1: Workspace Setup, Duplicate Marking & BQSR

**Script:** `scripts/01_prep_and_bqsr.sh`

## 1. Biological & Algorithmic Rationale

### Reference Sequence Dictionary (`CreateSequenceDictionary`)
While `samtools faidx` creates a `.fai` byte-offset index for rapid sequence retrieval, GATK/HTSJDK strictly validates that the BAM header (`@SQ`) and reference FASTA share identical contig names (`SN:20`), lengths (`LN:63025520`), and sort order before processing reads.

### Duplicate Marking (`MarkDuplicates`)
* **Artifact:** PCR amplification during library preparation can duplicate random polymerase errors across multiple identical fragments. Similarly, patterned flow cells can misread a single cluster as two adjacent clusters (*optical duplicates*).
* **Mechanism:** `MarkDuplicates` identifies read pairs sharing identical 5' mapping coordinates and sets the SAM bitwise flag `0x400` (`1024`) on all except the highest base-quality fragment. Reads are flagged rather than deleted so downstream tools ignore duplicates during variant calling while preserving coverage statistics.

### Base Quality Score Recalibration (`BaseRecalibrator` + `ApplyBQSR`)
* **Artifact:** Sequencing instruments emit empirical Phred quality scores ($Q = -10 \log_{10} P_{\text{error}}$), but systematic machine errors vary by sequencing cycle and dinucleotide context.
* **Mechanism:**
  1. `BaseRecalibrator` masks known polymorphic sites (`--known-sites` from dbSNP and Mills/1000G indels) so real biological variants are not penalized as sequencing errors.
  2. Across remaining sites, it models empirical mismatch rates stratified by Read Group (`@RG`), reported quality score, cycle, and dinucleotide context.
  3. `ApplyBQSR` rewrites Phred scores in the BAM file to reflect calibrated error probabilities.

### Execution & Duplication Metrics
```Bash
chmod +x scripts/01_prep_and_bqsr.sh
bash scripts/01_prep_and_bqsr.sh
```

### Inspect the duplication materics 

```bash
grep -A 2 "^LIBRARY" data/processed/na12878_chr20.dup_metrics.txt | column -t
```
```
LIBRARY  UNPAIRED_READS_EXAMINED  READ_PAIRS_EXAMINED  SECONDARY_OR_SUPPLEMENTARY_RDS  UNMAPPED_READS  UNPAIRED_READ_DUPLICATES  READ_PAIR_DUPLICATES  READ_PAIR_OPTICAL_DUPLICATES  PERCENT_DUPLICATION  ESTIMATED_LIBRARY_SIZE
1        50458                    6705731              102401                          50466           6801                      13116                 104                           0.002454             1725608882
```

Read Pairs Examined: 6,705,731

Read Pair Duplicates: 13,116 (104 optical duplicates)

Duplication Rate: 0.002454 (~0.25%), confirming high library complexity (ESTIMATED_LIBRARY_SIZE = 1,725,608,882).

__END__