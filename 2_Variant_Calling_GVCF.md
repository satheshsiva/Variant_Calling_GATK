# Step 2: Germline Variant Discovery in GVCF Mode

**Script:** `scripts/02_haplotypecaller_gvcf.sh`

## 1. Biological & Algorithmic Rationale

### Local De Novo Assembly (`HaplotypeCaller -ERC GVCF`)
Unlike pileup-based callers that evaluate base counts column-by-column, `HaplotypeCaller` performs local assembly:
1. **Active Region Detection:** Scans contig `20` for windows exhibiting mismatch or indel evidence above background sequencing noise.
2. **De Bruijn Graph Assembly:** Reconstructs candidate haplotypes directly from reads spanning each active region.
3. **Pair-HMM Likelihood Calculation:** Aligns every read against every candidate haplotype using a Pair Hidden Markov Model to compute accurate per-read genotype likelihoods—substantially improving indel and complex variant sensitivity.
4. **GVCF Reference Compression (`-ERC GVCF`):** Emits both candidate variant sites (with a symbolic `<NON_REF>` allele) and homozygous-reference confidence intervals (`END=...`), enabling scalable multi-sample joint genotyping without re-reading BAM files.

### Joint Genotyping (`GenotypeGVCFs`)
`GenotypeGVCFs` evaluates genotype likelihoods across samples, strips `<NON_REF>` placeholders and homozygous reference blocks, applies the call-confidence threshold (`QUAL >= 10`), and annotates known dbSNP rsIDs (`-D`).

---

## 2. Execution & Raw Call QC

```bash
chmod +x scripts/02_haplotypecaller_gvcf.sh
bash scripts/02_haplotypecaller_gvcf.sh
```

## QC metrics (SN/TSTV)

```Bash
bcftools stats Results/na12878_chr20.raw.vcf.gz | grep -E "^SN|^TSTV"
```
***What to Look For in the Output***  
*SN (Summary Numbers):* Number of records / SNPs / Indels: On GRCh37 chr20 (~63 Mb), a human genome typically yields roughly 80,000–110,000 raw variants (~75k–90k SNPs and ~10k–18k Indels).
    
*TSTV (Transition / Transversion Ratio):* Transitions ($A \leftrightarrow G$, $C \leftrightarrow T$) occur roughly twice as frequently as transversions in human genomes.For raw whole-chromosome human calls, expect a Ts/Tv ratio around ~2.0–2.1 (raw unfiltered sets often sit slightly lower, around 1.95–2.05, and rise closer to 2.05–2.15 after hard filtering removes false positives).

```
SN      0       number of samples:      1
SN      0       number of records:      110997
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 91874
SN      0       number of MNPs: 0
SN      0       number of indels:       19227
SN      0       number of others:       0
SN      0       number of multiallelic sites:   1543
SN      0       number of multiallelic SNP sites:       29
TSTV    0       61887   30016   2.06    61776   29886   2.07
```

Benchmark Output (na12878_chr20.raw.vcf.gz)
Total Records: 110,997 (91,874 SNPs, 19,227 Indels, 1,543 multiallelic sites)

Ts/Tv Ratio: 2.06 (61,887 transitions / 30,016 transversions)

(Note: na12878_chr20.g.vcf.gz contains 10,836,560 total records, including 10,717,437 homozygous reference blocks, and is retained only for multi-sample joint genotyping—never passed directly to filtering or VEP.)

__END__