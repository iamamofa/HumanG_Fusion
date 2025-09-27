// modules/utils/processes.nf
// Small helper processes used by other modules.
//
// Usage: include { fastqc_run; cutadapt_trim; fq_pair_to_tuple } from './modules/utils/processes.nf'

process fastqc_run {
  tag { sample_id ? "fastqc:${sample_id}" : "fastqc" }
  publishDir "${params.outdir ?: './results'}/qc", mode: 'copy', optional: true
  input:
    tuple val(sample_id), path(r1), path(r2)
  output:
    tuple val(sample_id), path("${sample_id}.R1_fastqc.zip") optional true, path("${sample_id}.R2_fastqc.zip") optional true
  script:
  """
  set -euo pipefail
  mkdir -p qc_${sample_id}
  fastqc -t ${task.cpus} -o qc_${sample_id} ${r1} ${r2} || true
  mv qc_${sample_id}/*_R1*zip ${sample_id}.R1_fastqc.zip || true
  mv qc_${sample_id}/*_R2*zip ${sample_id}.R2_fastqc.zip || true
  """
}

process cutadapt_trim {
  tag { "trim:${sample_id}" }
  publishDir "${params.outdir ?: './results'}/\${sample_id}/fastq_trimmed", mode: 'copy', optional: true
  input:
    tuple val(sample_id), path(r1), path(r2)
  output:
    tuple val(sample_id), path("${sample_id}_R1.trim.fastq.gz"), path("${sample_id}_R2.trim.fastq.gz")
  script:
  """
  set -euo pipefail
  cutadapt -q 20 -m 30 \
    -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
    -o ${sample_id}_R1.trim.fastq.gz -p ${sample_id}_R2.trim.fastq.gz ${r1} ${r2}
  """
}

process fq_pair_to_tuple {
  // Convenience: convert file pair inputs into a tuple for workflows
  cpus 1
  input:
    path files
  output:
    tuple path(files)
  script:
  """
  # passthrough process; Nextflow channels can be used directly but this is here for completeness
  ls -l ${files}
  """
}
