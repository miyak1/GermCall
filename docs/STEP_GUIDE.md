# GermCall step-by-step guide

For every stage: **what it does**, **why**, **what comes out** and **how to read it**. Stage names are the names you use with `--steps`.

The idea in one paragraph: sequencing gives millions of short reads per sample. We (1) check their quality, (2) cut off the bad parts, (3) place every read on the reference genome, (4) at every position, compare what the reads of a sample show with the reference to find that sample's variants, and (5) call all samples **together**, so that each sample gets a genotype at every variant site found in any sample.

---

## `prepare`
Checks that the programs, the sample table and all FASTQ files exist, and indexes the reference genome (BWA index, `.fai` and GATK `.dict`). These indices are needed by every later stage. Nothing to interpret, but read the printed line **Contigs: N | longest: ...**: it should match your genome (a chromosome-level assembly shows your chromosomes as the longest contigs).

## `qc`: FastQC on the raw reads
**What:** FastQC reads each FASTQ file and reports its quality. **Why:** to find problems (bad cycles, adapters, contamination) before they spread into the variant calls.
**Where:** `00_fastqc_raw/<sample>/*_fastqc.html` (open in a browser) and, for all samples at once, in the MultiQC report.
**How to read the main panels:**
- *Per base sequence quality*: quality (Phred score) along the read. Above 28 is good; a drop at the read end is normal and is removed by trimming. A sudden drop in the middle of all reads points to a sequencing problem.
- *Adapter content*: the share of reads containing sequencing adapters, which would later cause false mapping. Rising lines at the read end are common and are fixed by trimming.
- *Per sequence GC content*: should look like one smooth bell curve at your species' GC. A second peak suggests contamination (another organism).
- *Sequence duplication levels* and *overrepresented sequences*: a few sequences very frequent can mean adapters, rRNA/organelle DNA, or a low-complexity library.
- *Per base N content*: should be near 0 everywhere.
Warnings in FastQC are guidelines, not verdicts. Compare samples with each other: one sample that differs from all others is the one to inspect.

## `trim_align`: trimming, mapping, duplicates, statistics
**Trim Galore** removes adapters and low-quality ends (`trim_quality`, default 20) and drops reads that become too short (`trim_min_length`, default 50). The report says how many reads had adapters and how many bases were removed; losing more than about 10-20% of reads deserves a look at the raw FastQC.
**FastQC again** on the trimmed reads (`fastqc_after_trim`): adapter content should now be gone and the quality curve higher.
**BWA-MEM** aligns each read pair to the reference. The *read group* (`SM:<sample>`) written here is what names the sample in the VCF, so the sample name from your table is carried through all stages.
**MarkDuplicates** flags reads that come from the same original DNA fragment (PCR or optical duplicates): they would count the same molecule twice and make sites look better supported than they are.
**Statistics** (`02_bam/bam_summary.tsv`):
| Column | How to read it |
|---|---|
| Total_reads | reads after trimming that went into the alignment |
| Mapped_pct | 95% or more is typical for a good reference of the same species. **Below about 80%**: wrong reference, contamination, poor DNA, or a sample mix-up |
| Duplicate_pct | up to 10-20% is usual. Much higher means a low-complexity library (too little DNA, too many PCR cycles); duplicates are ignored, so effective depth is lower |
| Mean_depth | average number of reads covering a base. For genotyping: **10x or more** is comfortable; 5-10x works with more missing data; **below 5x** gives many uncertain or missing genotypes. This number also tells you whether pooled "missing" genotypes later will be common |

Compare samples: one with half the depth or a much lower mapping rate is an outlier to investigate or drop before calling.

## `gvcf`: HaplotypeCaller in GVCF mode
**What:** for each sample, GATK HaplotypeCaller examines regions where reads disagree with the reference, re-assembles the local haplotypes and decides the most likely genotype. With `-ERC GVCF` it writes **a record for every position**: the sites with a variant, and for the rest a block saying "this block of positions looks homozygous reference with this confidence".
**Why this form?** It stores, per sample, how much evidence there is at *every* position, including positions with no variant. That lets the joint step later tell "sample X is really reference here" apart from "sample X had no coverage here", which a plain variant-only VCF cannot do.
**Output:** `03_gvcf/<sample>.g.vcf.gz` (large; not for direct use).

## `joint`: joint genotyping
**What:** the gVCFs of all samples are combined (`GenomicsDBImport`, a database that scales to many samples; or `CombineGVCFs`, simpler but slower) and **GenotypeGVCFs** calls genotypes for all samples at once, one job per chromosome, after which the chromosomes are joined into one VCF.
**Why jointly?** Evidence is pooled across samples: a variant seen weakly in several samples is found; a site is called in all samples, so each sample gets a reference, heterozygous or alternative genotype, or a proper "missing" if it had no data. This avoids the main problem of merging single-sample VCFs.
**Output:** `04_joint/<prefix>_joint_raw.vcf.gz` (SNPs **and** indels, unfiltered) with `04_joint/<prefix>_joint_raw.stats.txt`.
**Interpret:** the printed line `Samples: N (expected N)` must match your sample table. The raw VCF contains many false positives; that is why the next stage filters.

## `filter`: hard filters on SNPs and indels
**What:** splits the raw VCF into SNPs and indels (they behave differently, so they get different thresholds) and removes the sites that fail GATK's generic quality criteria. A site is removed if it fails **any** one criterion. A site that lacks an annotation (for example no MQRankSum when nobody is heterozygous) is **not** removed by that criterion, which is GATK's default.

**SNPs**
| Filter | Removes SNPs with | Idea |
|---|---|---|
| QD < 2 | low quality **by depth** | weak call relative to the evidence |
| QUAL < 30 | low overall call confidence | |
| SOR > 3, FS > 60 | strand bias | variant seen mostly on one DNA strand: likely an artefact |
| MQ < 40 | low mapping quality | reads that map ambiguously (repeats) |
| MQRankSum < -12.5, ReadPosRankSum < -8 | alt-supporting reads are worse mapped, or sit at read ends | artefacts |

**Indels** (laxer on strand bias and read position, because indels are harder to call)
| Filter | Removes indels with |
|---|---|
| QD < 2, QUAL < 30 | same ideas as for SNPs |
| FS > 200 | strong strand bias |
| ReadPosRankSum < -20 | alt alleles near the read ends |

**Output:** `05_filtered/<prefix>_snps_pass.vcf.gz` (the file to give to GermGWAS), `<prefix>_indels_pass.vcf.gz`, the `*_filtered.vcf.gz` files (every site with its FILTER label, useful to see *why* a site was removed) and the statistics files.
**Interpret:** the screen shows `SNPs: N before filters -> M passing` and the same for indels. The `summary` stage lists how many sites each filter removed. The **transition/transversion ratio (Ts/Tv)** of the passing SNPs is a quality check: in plants it is typically around 2; a value near 0.5-1 means many random false calls. All thresholds are generic recommendations; if your distributions look very different, change `calling.filter_*` and `calling.filter_indel_*` in the config and rerun this stage.
Sites that carry both a SNP allele and an indel allele (GATK type MIXED) are in neither file, so `SNPs + indels` is a little smaller than the number of raw records; the `summary` report states how many.
These are *site-level* filters. Genotype-level filters (depth and quality of each genotype, missing data, minor allele frequency) are applied later by GermGWAS' QC step.

## `summary`: one report for the whole run
**What:** collects what the earlier stages measured and writes it in one place, `summary/run_summary.txt`, plus tables you can open in Excel:
| Section | Content |
|---|---|
| 1 Reference | file, total length, number of contigs, N50, GC %, unknown bases, chromosomes vs unplaced |
| 2 Raw reads | read pairs per sample (mean / min / max, total), raw bases, read length, GC % |
| 3 Trimming | bases kept, reads with adapters, pairs dropped for being too short |
| 4 Mapping | mapped %, properly paired %, duplicate %, mean depth, genome covered at least 1x, number of samples below 5x / 10x / 20x |
| 5 Variant calling | records, SNPs, indels, multi-allelic sites, Ts/Tv, SNPs per chromosome |
| 6 Hard filtering | every threshold (SNPs and indels), sites before / removed / passing, sites failing each filter |
| 7 Per-sample genotypes | missing-genotype % (raw and after SNP filtering), heterozygous %, worst samples |
| 8 Things to check | automatic warnings (low depth, low mapping rate, high missing data, odd heterozygosity, low Ts/Tv) |

**Tables** (all in `summary/`): `sample_summary.tsv` (one row per sample), `variant_summary.tsv`, `filter_criteria.tsv` (thresholds and the number of sites each removed), `variants_per_contig.tsv`, `reference_summary.tsv`.
**Run it** after `joint` (raw variants) and ideally after `filter` (so section 6 is filled): `--steps filter,summary`. It can be rerun at any time and only reads existing results.
**Interpret:** compare samples, not only averages. A sample with much lower depth, a high missing-genotype % or an unusually high heterozygous % is the one to look at before GWAS (mixed or contaminated material shows as high heterozygosity). Low-depth samples (below about 5x) give many uncertain genotypes; check section 7 and decide whether to keep them.

## `report`: MultiQC
One HTML page (`qc_reports/multiqc_report.html`) that puts FastQC, Trim Galore, mapping, duplicate and variant statistics of **all samples** side by side. Use it to spot the odd sample: in each plot, most samples form a group, and an outlier line is the one to investigate. Open it in a browser; no interpretation is needed beyond comparing samples.

---

## Handing over to GermGWAS
Give `05_filtered/<prefix>_snps_pass.vcf.gz` to GermGWAS as `input.merged_vcf`. GermGWAS then applies its own sample- and SNP-level filters (minor allele frequency, missingness, ...), so these two filter sets complement each other: GermCall removes technical artefacts, GermGWAS removes SNPs unsuitable for population analysis and association.
