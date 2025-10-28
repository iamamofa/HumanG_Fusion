#!/usr/bin/env bash
set -euo pipefail
source ../pipeline.config.sh
sample="$1"

aligndir="${OUTDIR}/${sample}/align"
fusiondir="${OUTDIR}/${sample}/fusion"
log="${LOGDIR}/${sample}_03_fusion.log"
mkdir -p "${fusiondir}" "$(dirname "${log}")"

echo "[$(date)] Fusion detection for ${sample}" > "${log}"

# STAR-Fusion
# STAR-Fusion expects raw FASTQ or STAR chimeric outputs depending on version; here we call with --chimeric_junction
# If using the STAR-Fusion wrapper that calls STAR itself:
STAR-Fusion --left_fq "${OUTDIR}/${sample}/fastq_trimmed/${sample}_R1.trim.fastq.gz" \
            --right_fq "${OUTDIR}/${sample}/fastq_trimmed/${sample}_R2.trim.fastq.gz" \
            --CPU ${THREADS} \
            --genome_lib_dir ${STAR_FUSION_CTAT_LIB} \
            --output_dir "${fusiondir}/starfusion" >> "${log}" 2>&1 || echo "STAR-Fusion failed" >> "${log}"

# ARRIBA (requires chimeric SAM/BAM produced by STAR and a blacklist)
# Need: chimeric.out.junction and aligned.bam from STAR (adjust paths)
CHIMERIC="${aligndir}/${sample}.Chimeric.out.junction"
ALIGNED_BAM="${aligndir}/${sample}.sorted.bam"
BLACKLIST="/refs/arriba/blacklist_human_GRCh38.tsv"   # adjust
ANNOTATION="${GTF}"                                   # arriba uses gene annotation tsv/bed — may need conversion
if [[ -f "${CHIMERIC}" && -f "${ALIGNED_BAM}" ]]; then
  ${ARRIBA} -x "${CHIMERIC}" -b "${ALIGNED_BAM}" -o "${fusiondir}/arriba.tsv" -O "${fusiondir}/arriba.discarded.tsv" -a "${GENOME_FASTA}" -g "${ANNOTATION}" -t "${BLACKLIST}" > "${fusiondir}/arriba.log" 2>&1 || echo "arriba failed" >> "${log}"
else
  echo "arriba: missing STAR outputs for ${sample}" >> "${log}"
fi

# FusionCatcher (optional - slower but complementary)
if command -v ${FUSIONCATCHER} >/dev/null 2>&1; then
  mkdir -p "${fusiondir}/fusioncatcher"
  ${FUSIONCATCHER} -d "${fusiondir}/fusioncatcher" -p ${THREADS} -i "${OUTDIR}/${sample}/fastq_trimmed" >> "${log}" 2>&1 || echo "FusionCatcher failed" >> "${log}"
fi

echo "[$(date)] Fusion callers finished for ${sample}" >> "${log}"
