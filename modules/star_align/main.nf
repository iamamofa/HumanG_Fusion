// modules/star_align/main.nf
nextflow.enable.dsl = 2

/*
 STAR 2-pass alignment module.

 Input: tuples (sample_id, r1, r2)  [these should be trimmed or decontaminated]
 Output: tuples (sample_id, sorted_bam, Chimeric.out.junction, Aligned.out.bam)
 Optionally emits metadata with original FASTQ paths (meta)
*/

process star_2pass {
  tag { "star:${sample_id}" }
  cpus { params.threads ?: 8 }
  time '12h'
  publishDir "${params.outdir ?: './results'}/\${sample_id}/align", mode: 'copy', optional: true

  input:
    tuple val(sample_id), path(r1), path(r2)

  output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.Chimeric.out.junction"), path("${sample_id}.Aligned.out.bam")

  script:
  """
  set -euo pipefail
  if [ -z "${params.star_index}" ]; then
    echo "ERROR: params.star_index must be set (STAR genome index path)" >&2
    exit 1
  fi

  mkdir -p star_tmp_${sample_id}
  # STAR first pass (we rely on STAR's 2-pass via --twopassMode Basic)
  STAR --runThreadN ${task.cpus} \
       --genomeDir ${params.star_index} \
       --readFilesIn ${r1} ${r2} \
       --readFilesCommand zcat \
       --outFileNamePrefix ${sample_id}. \
       --outSAMtype BAM Unsorted \
       --chimSegmentMin 12 \
       --chimJunctionOverhangMin 12 \
       --chimOutType Junctions \
       --outSAMattributes NH HI AS nM MD XS \
       --outTmpDir star_tmp_${sample_id} \
       --twopassMode Basic

  # sort and index
  samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam ${sample_id}.Aligned.out.bam
  samtools index ${sample_id}.sorted.bam || true

  # ensure chimeric junction file present (STAR writes Chimeric.out.junction)
  if [ ! -f Chimeric.out.junction ]; then
    echo "WARNING: Chimeric.out.junction not found for ${sample_id}" >&2
  fi

  """
}

workflow star_align {
  take: reads_in

  main:
    reads_in
      .map { tuple -> star_2pass(tuple[0], tuple[1], tuple[2]) }
      .flatten()
      .set { aligned_ch }

  emit:
    aligned_ch
}
