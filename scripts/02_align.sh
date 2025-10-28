#!/usr/bin/env bash
set -euo pipefail
source ../pipeline.config.sh
sample="$1"

outdir="${OUTDIR}/${sample}/align"
trimdir="${OUTDIR}/${sample}/fastq_trimmed"
log="${LOGDIR}/${sample}_02_align.log"
mkdir -p "${outdir}" "$(dirname "${log}")"

r1="${trimdir}/${sample}_R1.trim.fastq.gz"
r2="${trimdir}/${sample}_R2.trim.fastq.gz"

echo "[$(date)] STAR align for ${sample}" > "${log}"

# 2-pass STAR (write junctions to SJ.out.tab for 2nd pass)
STAR --runThreadN ${STAR_THREADS} \
  --genomeDir ${STAR_INDEX} \
  --readFilesIn "${r1}" "${r2}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${outdir}/${sample}." \
  --outSAMtype BAM Unsorted \
  --chimSegmentMin 12 \
  --chimJunctionOverhangMin 12 \
  --chimOutType Junctions \
  --outSAMattrRGline ID:${sample} SM:${sample} >> "${log}" 2>&1

# sort BAM
samtools sort -@ ${THREADS} -o "${outdir}/${sample}.sorted.bam" "${outdir}/${sample}.Aligned.out.bam"
samtools index "${outdir}/${sample}.sorted.bam"

echo "[$(date)] STAR align completed for ${sample}" >> "${log}"
