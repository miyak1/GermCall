#!/usr/bin/env bash
# =============================================================================
#  GermCall - from raw reads to a joint-called VCF
#
#  Stages (run all, or pick with --steps):
#    prepare     check tools / samples, index the reference (bwa, samtools faidx, GATK dict)
#    qc          FastQC on the raw reads (read quality before trimming)
#    trim_align  Trim Galore -> (FastQC on trimmed reads) -> BWA-MEM -> sort -> MarkDuplicates
#                -> mapping statistics (reads, mapped %, duplicates %, mean depth)
#    gvcf        GATK HaplotypeCaller -ERC GVCF, one gVCF per sample
#    joint       combine the gVCFs (GenomicsDBImport or CombineGVCFs) + GenotypeGVCFs,
#                one job per chromosome, then concatenate -> <prefix>_joint_raw.vcf.gz
#    filter      split into SNPs and indels, apply GATK's generic hard filters to each
#                -> <prefix>_snps_pass.vcf.gz and <prefix>_indels_pass.vcf.gz
#    summary     one text report + tables: reference, raw reads, trimming, mapping, variant
#                counts (SNPs / indels), filter criteria and how many variants each filter removed
#    report      MultiQC: one HTML report with all read, mapping and variant statistics
#
#  Default run: prepare,qc,trim_align,gvcf,joint.   --all: also filter, summary and report.
#
#  Usage
#    bash germcall.sh --config config/my_config.yaml
#    bash germcall.sh --config config/my_config.yaml --steps qc,trim_align
#    bash germcall.sh --config config/my_config.yaml --all
#
#  Every finished piece is remembered (<workdir>/status/), so re-running resumes.
#  Delete a status file to redo one piece. Logs: <workdir>/logs/
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/config.sh"
init_config "$@"

STEPS_ARG="prepare,qc,trim_align,gvcf,joint"
for ((i=0; i<${#EXTRA_ARGS[@]}; i++)); do
    case "${EXTRA_ARGS[$i]}" in
        --steps) STEPS_ARG="${EXTRA_ARGS[$((i+1))]}" ;;
        --all)   STEPS_ARG="prepare,qc,trim_align,gvcf,joint,filter,summary,report" ;;
    esac
done
want() { [[ ",${STEPS_ARG}," == *",$1,"* ]]; }

# ---- configuration --------------------------------------------------------
REF="$(path_of "$(cfg calling.reference '')")"
SAMPLES_FILE="$(path_of "$(cfg calling.samples_file '')")"
PAIRED="$(cfg calling.paired true)"
PAR_S="$(cfg calling.parallel_samples 2)"; THR_S="$(cfg calling.threads_per_sample 4)"
PAR_I="$(cfg calling.parallel_intervals 4)"; MEM="$(cfg calling.java_mem_gb 8)"
TRIM_Q="$(cfg calling.trim_quality 20)"; TRIM_L="$(cfg calling.trim_min_length 50)"
TRIM_X="$(cfg calling.trim_extra '')"; FQC_TRIM="$(cfg calling.fastqc_after_trim true)"
COMBINE="$(cfg calling.combine_method genomicsdb)"
N_CHR="$(cfg calling.n_chromosomes '')"; INC_UNPL="$(cfg calling.include_unplaced true)"; PLOIDY="$(cfg calling.ploidy 2)"
KEEP_TRIM="$(cfg calling.keep_trimmed_fastq false)"; KEEP_BAM="$(cfg calling.keep_unmarked_bam false)"
KEEP_DB="$(cfg calling.keep_genomicsdb false)"
TRIMGALORE="$(tool trim_galore trim_galore)"; FASTQC="$(tool fastqc fastqc)"; MULTIQC="$(tool multiqc multiqc)"; BWA="$(tool bwa bwa)"
SAMTOOLS="$(tool samtools samtools)"; GATK="$(tool gatk gatk)"; BCFTOOLS="$(tool bcftools bcftools)"
BGZIP="$(tool bgzip bgzip)"; TABIX="$(tool tabix tabix)"
mapfile -t CHR_LIST < <(cfg_list calling.chromosomes | sed '/^$/d')

D_FQC="${WORKDIR}/00_fastqc_raw"; D_REP="${WORKDIR}/qc_reports"; D_TRIM="${WORKDIR}/01_trimmed"; D_BAM="${WORKDIR}/02_bam"; D_GVCF="${WORKDIR}/03_gvcf"
D_JOINT="${WORKDIR}/04_joint"; D_FILT="${WORKDIR}/05_filtered"; D_SUM="${WORKDIR}/summary"; D_LOG="${WORKDIR}/logs"; D_ST="${WORKDIR}/status"
mkdir -p "${D_FQC}" "${D_REP}" "${D_TRIM}" "${D_BAM}" "${D_GVCF}" "${D_JOINT}" "${D_FILT}" "${D_SUM}" "${D_LOG}" "${D_ST}"
GATK_OPTS="-Xmx${MEM}G -Djava.io.tmpdir=${WORKDIR}/tmp"; mkdir -p "${WORKDIR}/tmp"
run_gatk() { "${GATK}" --java-options "${GATK_OPTS}" "$@"; }

safe() { echo "$1" | sed 's/[^A-Za-z0-9._-]/_/g'; }
done_flag() { echo "${D_ST}/$1.done"; }

# run "$fn item" for every item, at most $max at a time; returns non-zero if any job failed
run_parallel() {
    local max="$1" fn="$2"; shift 2; local failed=0 item
    for item in "$@"; do
        "${fn}" "${item}" &
        while (( $(jobs -rp | wc -l) >= max )); do wait -n || failed=1; done
    done
    while (( $(jobs -rp | wc -l) > 0 )); do wait -n || failed=1; done
    return ${failed}
}

# ---- read the sample table ------------------------------------------------
declare -A R1 R2; SAMPLES=()
read_samples() {
    R1=(); R2=(); SAMPLES=()
    [[ -f "${SAMPLES_FILE}" ]] || { echo "ERROR: calling.samples_file not found: '${SAMPLES_FILE}'" >&2; exit 1; }
    local s a b
    while IFS=$'\t ' read -r s a b _; do
        [[ -z "${s}" || "${s}" == \#* || "${s}" == "sample" ]] && continue
        [[ "${s}" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "ERROR: sample name '${s}' may only contain letters, digits . _ -" >&2; exit 1; }
        [[ -n "${R1[$s]-}" ]] && { echo "ERROR: sample '${s}' is listed twice" >&2; exit 1; }
        [[ -f "${a}" ]] || { echo "ERROR: R1 file not found for ${s}: ${a}" >&2; exit 1; }
        if [[ "${PAIRED}" == "true" ]]; then
            [[ -f "${b}" ]] || { echo "ERROR: R2 file not found for ${s}: '${b}' (calling.paired is true)" >&2; exit 1; }
            R2[$s]="${b}"
        fi
        R1[$s]="${a}"; SAMPLES+=("${s}")
    done < "${SAMPLES_FILE}"
    [[ ${#SAMPLES[@]} -gt 0 ]] || { echo "ERROR: no samples in ${SAMPLES_FILE}" >&2; exit 1; }
}

REF_DICT="${REF%.*}.dict"

# ---- stage: prepare -------------------------------------------------------
stage_prepare() {
    banner "[prepare] checking tools, samples and reference"
    for t in "${BWA}" "${SAMTOOLS}" "${GATK}" "${BCFTOOLS}" "${BGZIP}" "${TABIX}"; do require_tool "${t}"; done
    require_tool "${TRIMGALORE}" "conda install -c bioconda trim-galore cutadapt"
    if want qc || [[ "${FQC_TRIM}" == "true" ]]; then require_tool "${FASTQC}" "conda install -c bioconda fastqc   (or set calling.fastqc_after_trim: false and skip stage qc)"; fi
    want report && require_tool "${MULTIQC}" "conda install -c bioconda multiqc"
    require_file "${REF}" "Set calling.reference"
    read_samples
    echo "  Samples   : ${#SAMPLES[@]} ($([[ "${PAIRED}" == "true" ]] && echo paired-end || echo single-end))"
    echo "  Reference : ${REF}"
    [[ -f "${REF}.bwt" || -f "${REF}.bwt.2bit.64" ]] || { echo "  indexing reference with bwa ..."; "${BWA}" index "${REF}" > "${D_LOG}/bwa_index.log" 2>&1; }
    [[ -f "${REF}.fai" ]] || "${SAMTOOLS}" faidx "${REF}"
    [[ -f "${REF_DICT}" ]] || run_gatk CreateSequenceDictionary -R "${REF}" -O "${REF_DICT}" > "${D_LOG}/dict.log" 2>&1
    echo "  Reference indices ready (bwa, .fai, .dict)"
    echo "  Contigs: $(wc -l < "${REF}.fai")  | longest: $(sort -k2,2nr "${REF}.fai" | head -3 | cut -f1,2 | tr '\t' ':' | tr '\n' ' ')"
}

# ---- stage: read QC (raw reads) -------------------------------------------------
fastqc_raw_one() {
    local s="$1" flag; flag="$(done_flag "fastqc_raw_${s}")"
    [[ -f "${flag}" ]] && { echo "  [${s}] raw FastQC already done - skipping"; return 0; }
    {
        mkdir -p "${D_FQC}/${s}"
        "${FASTQC}" -t "${THR_S}" -o "${D_FQC}/${s}" "${R1[$s]}" ${R2[$s]:+"${R2[$s]}"}
        touch "${flag}"
    } > "${D_LOG}/fastqc_raw_${s}.log" 2>&1 || { echo "  [${s}] FAILED - see ${D_LOG}/fastqc_raw_${s}.log" >&2; return 1; }
    echo "  [${s}] raw FastQC done"
}
stage_qc() {
    banner "[qc] FastQC on raw reads | ${#SAMPLES[@]} samples, ${PAR_S} at a time"
    run_parallel "${PAR_S}" fastqc_raw_one "${SAMPLES[@]}" || { echo "ERROR: some FastQC jobs failed (see logs)" >&2; exit 1; }
    echo "  Reports: ${D_FQC}/<sample>/*_fastqc.html   (the 'report' stage summarises all samples in one MultiQC page)"
}

# ---- stage: trim + align + markdup ------------------------------------------
trim_align_one() {
    local s="$1" flag; flag="$(done_flag "trim_align_${s}")"
    [[ -f "${flag}" ]] && { echo "  [${s}] trim/align already done - skipping"; return 0; }
    {
        echo ">>> ${s}: trimming"
        local tdir="${D_TRIM}/${s}"; mkdir -p "${tdir}"
        local tg=(--quality "${TRIM_Q}" --length "${TRIM_L}" --cores "${THR_S}" --basename "${s}" -o "${tdir}")
        local in1 in2=""
        if [[ "${PAIRED}" == "true" ]]; then
            # shellcheck disable=SC2086
            "${TRIMGALORE}" --paired "${tg[@]}" ${TRIM_X} "${R1[$s]}" "${R2[$s]}"
            in1="${tdir}/${s}_val_1.fq.gz"; in2="${tdir}/${s}_val_2.fq.gz"
        else
            # shellcheck disable=SC2086
            "${TRIMGALORE}" "${tg[@]}" ${TRIM_X} "${R1[$s]}"
            in1="${tdir}/${s}_trimmed.fq.gz"
        fi
        [[ -s "${in1}" ]] || { echo "ERROR: trimmed reads not found: ${in1}" >&2; exit 1; }
        if [[ "${FQC_TRIM}" == "true" ]]; then "${FASTQC}" -t "${THR_S}" -o "${tdir}" "${in1}" ${in2:+"${in2}"}; fi

        echo ">>> ${s}: BWA-MEM alignment + sort"
        local rg="@RG\\tID:${s}\\tSM:${s}\\tPL:ILLUMINA\\tLB:${s}"
        "${BWA}" mem -M -t "${THR_S}" -R "${rg}" "${REF}" "${in1}" ${in2:+"${in2}"} \
            | "${SAMTOOLS}" sort -@ "${THR_S}" -T "${WORKDIR}/tmp/${s}.sort" -o "${D_BAM}/${s}.sorted.bam" -
        "${SAMTOOLS}" index -@ "${THR_S}" "${D_BAM}/${s}.sorted.bam"

        echo ">>> ${s}: MarkDuplicates"
        run_gatk MarkDuplicates -I "${D_BAM}/${s}.sorted.bam" -O "${D_BAM}/${s}.dedup.bam" \
            -M "${D_BAM}/${s}.dup_metrics.txt"
        "${SAMTOOLS}" index -@ "${THR_S}" "${D_BAM}/${s}.dedup.bam"

        echo ">>> ${s}: statistics"
        mkdir -p "${D_BAM}/stats"
        "${SAMTOOLS}" flagstat -@ "${THR_S}" "${D_BAM}/${s}.dedup.bam" > "${D_BAM}/stats/${s}.flagstat.txt"
        "${SAMTOOLS}" coverage "${D_BAM}/${s}.dedup.bam" > "${D_BAM}/stats/${s}.coverage.txt"

        [[ "${KEEP_BAM}" == "true" ]]  || rm -f "${D_BAM}/${s}.sorted.bam" "${D_BAM}/${s}.sorted.bam.bai"
        [[ "${KEEP_TRIM}" == "true" ]] || rm -rf "${tdir}"/*.fq.gz
        touch "${flag}"
        echo ">>> ${s}: done"
    } > "${D_LOG}/trim_align_${s}.log" 2>&1 || { echo "  [${s}] FAILED - see ${D_LOG}/trim_align_${s}.log" >&2; return 1; }
    echo "  [${s}] trim/align done"
}

bam_summary() {
    local out="${D_BAM}/bam_summary.tsv"
    printf 'Sample\tTotal_reads\tMapped_pct\tDuplicate_pct\tMean_depth\n' > "${out}"
    local s f c
    for s in "${SAMPLES[@]}"; do
        f="${D_BAM}/stats/${s}.flagstat.txt"; c="${D_BAM}/stats/${s}.coverage.txt"
        [[ -f "${f}" && -f "${c}" ]] || continue
        awk -v s="${s}" 'FILENAME==ARGV[1] {
               if ($0 ~ / in total/) tot=$1; if ($0 ~ / duplicates/) dup=$1;
               if ($0 ~ / mapped \(/) { t2=$0; sub(/^.*\(/,"",t2); sub(/%.*$/,"",t2); mp=t2 } next }
             FNR>1 { len=$3-$2+1; w+=$7*len; L+=len }
             END { printf "%s\t%d\t%s\t%.2f\t%.2f\n", s, tot, (mp==""?"NA":mp), (tot>0?100*dup/tot:0), (L>0?w/L:0) }' "${f}" "${c}" >> "${out}"
    done
    echo "  Mapping summary: ${out}"
    awk -F'\t' '{printf "  %-12s %12s %11s %14s %11s\n",$1,$2,$3,$4,$5}' "${out}" | head -16
    local low; low="$(awk -F'\t' 'NR>1 && $5<5 {n++} END{print n+0}' "${out}")"
    [[ "${low}" -gt 0 ]] && echo "  NOTE: ${low} sample(s) have mean depth < 5x; low-depth samples give many uncertain genotypes."
    return 0
}

stage_trim_align() {
    banner "[trim_align] Trim Galore -> BWA-MEM -> sort -> MarkDuplicates | ${#SAMPLES[@]} samples, ${PAR_S} at a time"
    run_parallel "${PAR_S}" trim_align_one "${SAMPLES[@]}" || { echo "ERROR: some samples failed (see logs)" >&2; exit 1; }
    bam_summary
}

# ---- stage: gVCF ----------------------------------------------------------------
gvcf_one() {
    local s="$1" flag; flag="$(done_flag "gvcf_${s}")"
    [[ -f "${flag}" ]] && { echo "  [${s}] gVCF already done - skipping"; return 0; }
    [[ -f "${D_BAM}/${s}.dedup.bam" ]] || { echo "ERROR: ${D_BAM}/${s}.dedup.bam missing - run stage trim_align" >&2; return 1; }
    {
        run_gatk HaplotypeCaller -R "${REF}" -I "${D_BAM}/${s}.dedup.bam" \
            -O "${D_GVCF}/${s}.g.vcf.gz" -ERC GVCF -ploidy "${PLOIDY}" --native-pair-hmm-threads "${THR_S}"
        [[ -f "${D_GVCF}/${s}.g.vcf.gz.tbi" ]] || "${TABIX}" -p vcf "${D_GVCF}/${s}.g.vcf.gz"
        touch "${flag}"
    } > "${D_LOG}/gvcf_${s}.log" 2>&1 || { echo "  [${s}] FAILED - see ${D_LOG}/gvcf_${s}.log" >&2; return 1; }
    echo "  [${s}] gVCF done"
}
stage_gvcf() {
    banner "[gvcf] HaplotypeCaller -ERC GVCF | ${#SAMPLES[@]} samples, ${PAR_S} at a time"
    run_parallel "${PAR_S}" gvcf_one "${SAMPLES[@]}" || { echo "ERROR: some gVCFs failed (see logs)" >&2; exit 1; }
}

# ---- stage: joint genotyping -------------------------------------------------------
SAMPLE_MAP="${D_JOINT}/gvcf_sample_map.txt"
joint_one() {
    local c="$1" cs; cs="$(safe "${c}")"; local flag; flag="$(done_flag "joint_${cs}")"
    [[ -f "${flag}" ]] && { echo "  [${c}] already genotyped - skipping"; return 0; }
    {
        mkdir -p "${D_JOINT}/by_contig" "${WORKDIR}/tmp/${cs}"
        local vin
        if [[ "${COMBINE}" == "genomicsdb" ]]; then
            rm -rf "${D_JOINT}/db/${cs}"; mkdir -p "${D_JOINT}/db"
            run_gatk GenomicsDBImport --genomicsdb-workspace-path "${D_JOINT}/db/${cs}" \
                --sample-name-map "${SAMPLE_MAP}" -L "${c}" --reader-threads "${THR_S}" --tmp-dir "${WORKDIR}/tmp/${cs}"
            vin="gendb://${D_JOINT}/db/${cs}"
        else
            mkdir -p "${D_JOINT}/combined"
            local vargs=() s
            for s in "${SAMPLES[@]}"; do vargs+=(-V "${D_GVCF}/${s}.g.vcf.gz"); done
            run_gatk CombineGVCFs -R "${REF}" "${vargs[@]}" -L "${c}" -O "${D_JOINT}/combined/${cs}.g.vcf.gz"
            vin="${D_JOINT}/combined/${cs}.g.vcf.gz"
        fi
        run_gatk GenotypeGVCFs -R "${REF}" -V "${vin}" -L "${c}" -ploidy "${PLOIDY}" -O "${D_JOINT}/by_contig/${cs}.vcf.gz" --tmp-dir "${WORKDIR}/tmp/${cs}"
        [[ "${KEEP_DB}" == "true" ]] || rm -rf "${D_JOINT}/db/${cs}" "${D_JOINT}/combined/${cs}.g.vcf.gz"*
        rm -rf "${WORKDIR}/tmp/${cs}"
        touch "${flag}"
    } > "${D_LOG}/joint_${cs}.log" 2>&1 || { echo "  [${c}] FAILED - see ${D_LOG}/joint_${cs}.log" >&2; return 1; }
    echo "  [${c}] genotyped"
}
stage_joint() {
    banner "[joint] combine gVCFs (${COMBINE}) + GenotypeGVCFs, one job per contig, ${PAR_I} at a time"
    [[ -f "${REF}.fai" ]] || { echo "ERROR: ${REF}.fai missing - run stage prepare" >&2; exit 1; }
    local s; : > "${SAMPLE_MAP}"
    for s in "${SAMPLES[@]}"; do
        [[ -f "${D_GVCF}/${s}.g.vcf.gz" ]] || { echo "ERROR: gVCF missing for ${s} - run stage gvcf" >&2; exit 1; }
        printf '%s\t%s\n' "${s}" "${D_GVCF}/${s}.g.vcf.gz" >> "${SAMPLE_MAP}"
    done
    # contigs: explicit list, or the N longest as "chromosomes"; the rest only if include_unplaced
    local chroms=() others=() all c
    declare -A INFAI ISCHR
    mapfile -t all < <(cut -f1 "${REF}.fai")
    for c in "${all[@]}"; do INFAI[$c]=1; done
    if [[ ${#CHR_LIST[@]} -gt 0 ]]; then
        chroms=("${CHR_LIST[@]}")
        for c in "${chroms[@]}"; do [[ -n "${INFAI[$c]-}" ]] || { echo "ERROR: contig '${c}' (calling.chromosomes) is not in ${REF}.fai" >&2; exit 1; }; done
    else
        [[ "${N_CHR}" =~ ^[0-9]+$ && "${N_CHR}" -gt 0 ]] || {
            echo "ERROR: set calling.n_chromosomes (the number of chromosomes of your species, e.g. 10) or list them in calling.chromosomes" >&2; exit 1; }
        mapfile -t chroms < <(sort -k2,2nr "${REF}.fai" | head -n "${N_CHR}" | cut -f1)
    fi
    for c in "${chroms[@]}"; do ISCHR[$c]=1; done
    for c in "${all[@]}"; do [[ -n "${ISCHR[$c]-}" ]] || others+=("${c}"); done
    local jobs=("${chroms[@]}"); [[ "${INC_UNPL}" == "true" ]] && jobs+=("${others[@]}")
    declare -A ISJOB; for c in "${jobs[@]}"; do ISJOB[$c]=1; done
    echo "  Chromosomes (${#chroms[@]}): ${chroms[*]}"
    echo "  Other contigs: ${#others[@]} ($([[ "${INC_UNPL}" == "true" ]] && echo included || echo skipped))"
    run_parallel "${PAR_I}" joint_one "${jobs[@]}" || { echo "ERROR: some contigs failed (see logs)" >&2; exit 1; }

    # concatenate in reference order
    local list="${D_JOINT}/concat_list.txt"; : > "${list}"
    for c in "${all[@]}"; do
        [[ -n "${ISJOB[$c]-}" ]] && echo "${D_JOINT}/by_contig/$(safe "${c}").vcf.gz" >> "${list}"
    done
    local raw="${D_JOINT}/${PREFIX}_joint_raw.vcf.gz"
    "${BCFTOOLS}" concat -f "${list}" -O z -o "${raw}"
    "${BCFTOOLS}" index -f "${raw}"
    "${BCFTOOLS}" stats -s - "${raw}" > "${D_JOINT}/${PREFIX}_joint_raw.stats.txt"
    local ns; ns="$("${BCFTOOLS}" query -l "${raw}" | wc -l)"
    echo ""; banner "JOINT CALLING COMPLETE"
    echo "  Joint VCF : ${raw}"
    echo "  Samples   : ${ns} (expected ${#SAMPLES[@]})"
    local st="${D_JOINT}/${PREFIX}_joint_raw.stats.txt" nrec nsnp nind
    nrec="$(awk -F'\t' '$1=="SN" && $3=="number of records:" {print $4}' "${st}")"
    nsnp="$(awk -F'\t' '$1=="SN" && $3=="number of SNPs:" {print $4}' "${st}")"
    nind="$(awk -F'\t' '$1=="SN" && $3=="number of indels:" {print $4}' "${st}")"
    echo "  Records   : ${nrec}  (SNPs: ${nsnp}, indels: ${nind})"
    [[ "${ns}" -eq "${#SAMPLES[@]}" ]] || echo "  WARNING: sample count differs from the sample table!"
    echo "  NEXT      : hard filtering + summary report:  --steps filter,summary"
}

# ---- stage: hard filter ----------------------------------------------------------------
# GATK "generic hard filter" recommendations, applied separately to SNPs and indels (they
# behave differently). A site that lacks an annotation (e.g. no MQRankSum when nobody is
# heterozygous) is NOT removed by that filter - GATK's default.
stage_filter() {
    banner "[filter] SNPs and indels: GATK hard filters"
    local raw="${D_JOINT}/${PREFIX}_joint_raw.vcf.gz"; require_file "${raw}" "Run stage joint first."
    local snps="${D_FILT}/${PREFIX}_snps_raw.vcf.gz"   sflt="${D_FILT}/${PREFIX}_snps_filtered.vcf.gz"   spass="${D_FILT}/${PREFIX}_snps_pass.vcf.gz"
    local inds="${D_FILT}/${PREFIX}_indels_raw.vcf.gz" iflt="${D_FILT}/${PREFIX}_indels_filtered.vcf.gz" ipass="${D_FILT}/${PREFIX}_indels_pass.vcf.gz"
    run_gatk SelectVariants -R "${REF}" -V "${raw}" --select-type-to-include SNP -O "${snps}"
    run_gatk VariantFiltration -R "${REF}" -V "${snps}" -O "${sflt}" \
        --filter-name "QD" --filter-expression "QD < $(cfg calling.filter_QD 2.0)" \
        --filter-name "QUAL" --filter-expression "QUAL < $(cfg calling.filter_QUAL 30.0)" \
        --filter-name "SOR" --filter-expression "SOR > $(cfg calling.filter_SOR 3.0)" \
        --filter-name "FS" --filter-expression "FS > $(cfg calling.filter_FS 60.0)" \
        --filter-name "MQ" --filter-expression "MQ < $(cfg calling.filter_MQ 40.0)" \
        --filter-name "MQRankSum" --filter-expression "MQRankSum < $(cfg calling.filter_MQRankSum -12.5)" \
        --filter-name "ReadPosRankSum" --filter-expression "ReadPosRankSum < $(cfg calling.filter_ReadPosRankSum -8.0)"
    run_gatk SelectVariants -R "${REF}" -V "${sflt}" --exclude-filtered -O "${spass}"
    "${BCFTOOLS}" stats -s - "${spass}" > "${D_FILT}/${PREFIX}_snps_pass.stats.txt"

    run_gatk SelectVariants -R "${REF}" -V "${raw}" --select-type-to-include INDEL -O "${inds}"
    run_gatk VariantFiltration -R "${REF}" -V "${inds}" -O "${iflt}" \
        --filter-name "QD" --filter-expression "QD < $(cfg calling.filter_indel_QD 2.0)" \
        --filter-name "QUAL" --filter-expression "QUAL < $(cfg calling.filter_indel_QUAL 30.0)" \
        --filter-name "FS" --filter-expression "FS > $(cfg calling.filter_indel_FS 200.0)" \
        --filter-name "ReadPosRankSum" --filter-expression "ReadPosRankSum < $(cfg calling.filter_indel_ReadPosRankSum -20.0)"
    run_gatk SelectVariants -R "${REF}" -V "${iflt}" --exclude-filtered -O "${ipass}"
    "${BCFTOOLS}" stats -s - "${ipass}" > "${D_FILT}/${PREFIX}_indels_pass.stats.txt"

    echo ""; banner "FILTERING COMPLETE"
    echo "  SNPs   : $("${BCFTOOLS}" view -H "${snps}" | wc -l) before filters -> $("${BCFTOOLS}" view -H "${spass}" | wc -l) passing"
    echo "  Indels : $("${BCFTOOLS}" view -H "${inds}" | wc -l) before filters -> $("${BCFTOOLS}" view -H "${ipass}" | wc -l) passing"
    echo "  SNP VCF   : ${spass}"
    echo "  Indel VCF : ${ipass}"
    echo "  Use it in GermGWAS  : input.merged_vcf: ${spass}   (skip the 'merge' step)"
    echo "  NEXT      : summary report with the filter criteria:  --steps summary"
}

# ---- stage: summary report ---------------------------------------------------------------
stage_summary() {
    banner "[summary] reference, reads, trimming, mapping, variants, filters"
    require_tool "${BCFTOOLS}"; require_tool "${SAMTOOLS}"
    [[ ${#SAMPLES[@]} -gt 0 ]] || read_samples
    printf '%s\n' "${SAMPLES[@]}" > "${D_ST}/sample_names.txt"
    python3 "${GC_ROOT}/lib/make_summary.py" --config "${CONFIG}" --workdir "${WORKDIR}" --prefix "${PREFIX}" \
        --ref "${REF}" --samples "${D_ST}/sample_names.txt" --bcftools "${BCFTOOLS}" --samtools "${SAMTOOLS}" \
        --bwa "${BWA}" --gatk "${GATK}" --trimgalore "${TRIMGALORE}" --out "${D_SUM}"
    echo "  Report : ${D_SUM}/run_summary.txt"
    echo "  Tables : ${D_SUM}/*.tsv   (open in Excel / LibreOffice)"
}

# ---- stage: MultiQC report -------------------------------------------------------------
stage_report() {
    banner "[report] MultiQC: reads, trimming, mapping, duplicates, variants"
    mkdir -p "${D_REP}"
    local inputs=() d
    for d in "${D_FQC}" "${D_TRIM}" "${D_BAM}" "${D_JOINT}" "${D_FILT}"; do [[ -d "${d}" ]] && inputs+=("${d}"); done
    "${MULTIQC}" -f -o "${D_REP}" -n multiqc_report "${inputs[@]}" > "${D_LOG}/multiqc.log" 2>&1 \
        || { echo "ERROR: MultiQC failed - see ${D_LOG}/multiqc.log" >&2; exit 1; }
    echo "  Report: ${D_REP}/multiqc_report.html   (open in a web browser)"
}

# ---- run ----------------------------------------------------------------------------------------
banner "GermCall | ${PREFIX} | steps: ${STEPS_ARG}"
[[ ${#SAMPLES[@]} -gt 0 ]] || read_samples
want prepare    && stage_prepare
want qc         && stage_qc
want trim_align && stage_trim_align
want gvcf       && stage_gvcf
want joint      && stage_joint
want filter     && stage_filter
want summary    && stage_summary
want report     && stage_report
echo ""; echo "All requested stages finished."
