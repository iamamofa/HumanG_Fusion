#!/usr/bin/env bash
set -euo pipefail
source pipeline.config.sh

mkdir -p "${OUTDIR}" "${LOGDIR}" "${WORKDIR}" "${OUTDIR}/reports"

# helper: find pairs if samples file not provided
find_pairs() {
  # looks for *_R1*.fastq* and pairs with _R2*
  for r1 in ${SAMPLES_DIR}/*_R1*.fastq*; do
    sample=$(basename "${r1}" | sed -E 's/_R1.*//')
    r2="${SAMPLES_DIR}/${sample}_R2.fastq.gz"
    if [[ -f "${r2}" ]]; then
      echo -e "${sample}\t${r1}\t${r2}"
    else
      echo "Warning: R2 not found for ${sample}" >&2
    fi
  done
}

if [[ -n "${SAMPLES_FILE}" && -f "${SAMPLES_FILE}" ]]; then
  MAPFILE=("${SAMPLES_FILE}")
else
  # create temporary samples list
  TMP_SAMPLES="${WORKDIR}/samples.tsv"
  find_pairs > "${TMP_SAMPLES}"
  MAPFILE=("${TMP_SAMPLES}")
fi

# iterate
while IFS=$'\t' read -r sample r1 r2; do
  echo "=== Processing ${sample} ==="
  mkdir -p "${OUTDIR}/${sample}" "${OUTDIR}/${sample}/fastq_trimmed" "${OUTDIR}/${sample}/align" "${OUTDIR}/${sample}/fusion"

  # step 1: preprocess
  bash scripts/01_preprocess.sh "${sample}" "${r1}" "${r2}" || { echo "Preprocess failed for ${sample}"; exit 1; }

  # step 2: align
  bash scripts/02_align.sh "${sample}" || { echo "Align failed for ${sample}"; exit 1; }

  # step 3: fusion detection
  bash scripts/03_fusion_detection.sh "${sample}" || { echo "Fusion detection failed for ${sample}"; exit 1; }

  # step 4: postprocess & merge
  bash scripts/04_postprocess.sh "${sample}" || { echo "Postprocess failed for ${sample}"; exit 1; }

  # step 5: report
  bash scripts/05_report.sh "${sample}" || { echo "Report failed for ${sample}"; exit 1; }

  echo "=== Done ${sample} ==="
done < <(cat "${MAPFILE}")

echo "All samples processed. Reports in ${OUTDIR}/reports"
