#!/usr/bin/env bash
set -euo pipefail
source ../pipeline.config.sh
sample="$1"
r1="$2"
r2="$3"

outdir="${OUTDIR}/${sample}/fastq_trimmed"
log="${LOGDIR}/${sample}_01_preprocess.log"
mkdir -p "${outdir}" "$(dirname "${log}")"

echo "[$(date)] Preprocessing ${sample}" > "${log}"

# fastqc raw
${FASTQC} -t ${THREADS} -o "${outdir}" "${r1}" "${r2}" >> "${log}" 2>&1 || true

# trimming (paired)
trim_r1="${outdir}/${sample}_R1.trim.fastq.gz"
trim_r2="${outdir}/${sample}_R2.trim.fastq.gz"

${CUTADAPT} ${TRIM_ARGS} -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o "${trim_r1}" -p "${trim_r2}" "${r1}" "${r2}" >> "${log}" 2>&1

# fastqc trimmed
${FASTQC} -t ${THREADS} -o "${outdir}" "${trim_r1}" "${trim_r2}" >> "${log}" 2>&1 || true

echo "[$(date)] Preprocessing done for ${sample}" >> "${log}"
