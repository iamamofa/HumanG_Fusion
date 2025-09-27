#!/usr/bin/env bash
set -euo pipefail
source ../pipeline.config.sh
sample="$1"

fusiondir="${OUTDIR}/${sample}/fusion"
out="${fusiondir}/${sample}.merged_fusions.tsv"
log="${LOGDIR}/${sample}_04_postprocess.log"
mkdir -p "$(dirname "${log}")"

echo "[$(date)] Postprocessing ${sample}" > "${log}"

python3 ../tools/merge_fusions.py \
  --starfusion "${fusiondir}/starfusion/star-fusion.fusion_predictions.tsv" \
  --arriba "${fusiondir}/arriba.tsv" \
  --fusioncatcher "${fusiondir}/fusioncatcher/*/*/final-list_candidate-fusion-genes.txt" \
  --gtf "${GTF}" \
  --out "${out}" >> "${log}" 2>&1 || echo "Merge script had issues" >> "${log}"

echo "[$(date)] Postprocess done. Output: ${out}" >> "${log}"
