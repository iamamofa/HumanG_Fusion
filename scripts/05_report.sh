#!/usr/bin/env bash
set -euo pipefail
sample="$1"
source ../pipeline.config.sh

fusionfile="${OUTDIR}/${sample}/fusion/${sample}.merged_fusions.tsv"
reportdir="${OUTDIR}/reports"
mkdir -p "${reportdir}"

# create a minimal report
report="${reportdir}/${sample}_fusion_report.txt"
echo "Fusion Report for ${sample}" > "${report}"
echo "Generated: $(date --utc)" >> "${report}"
echo "" >> "${report}"
if [[ -f "${fusionfile}" ]]; then
  echo "Merged fusion calls:" >> "${report"
  cat "${fusionfile}" >> "${report}"
else
  echo "No merged fusion file found for ${sample}" >> "${report}"
fi

# copy summary to reports dir
cp "${fusionfile}" "${reportdir}/${sample}.merged_fusions.tsv" 2>/dev/null || true
