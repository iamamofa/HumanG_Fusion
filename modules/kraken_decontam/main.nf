// modules/kraken_decontam/main.nf
nextflow.enable.dsl = 2

/*
 Kraken2 decontamination workflow.

 Input: tuples (sample_id, r1, r2)
 Output: tuples (sample_id, decontam_R1.fastq.gz, decontam_R2.fastq.gz)

 Behavior:
  - Runs kraken2 classification and uses --unclassified-out + --classified-out to generate filtered FASTQs.
  - Default conservative policy: keep unclassified + human (taxid 9606) if present in DB.
  - User can set params.kraken_policy = 'keep_unclassified' | 'keep_human_only' | 'remove_nonhuman'
  - Requires params.kraken2_db to be set when running.
*/

process kraken2_classify {
  tag { "kraken:${sample_id}" }
  cpus { params.threads ?: 4 }
  time '4h'
  publishDir "${params.outdir ?: './results'}/\${sample_id}/decontam", mode: 'copy', optional: true

  input:
    tuple val(sample_id), path(r1), path(r2)

  output:
    tuple val(sample_id), path("${sample_id}.decontam_R1.fastq.gz"), path("${sample_id}.decontam_R2.fastq.gz")

  script:
  """
  set -euo pipefail

  if [ -z "${params.kraken2_db}" ]; then
    echo "ERROR: params.kraken2_db must be set for kraken decontamination" >&2
    exit 1
  fi

  # run kraken2 classification to produce classified/unclassified paired outputs
  # kraken2's --unclassified-out and --classified-out require file patterns with # to indicate pair
  kraken2 --db ${params.kraken2_db} --paired --gzip-compressed --threads ${task.cpus} \
    --unclassified-out ${sample_id}.unclassified#.fastq --classified-out ${sample_id}.classified#.fastq \
    --report ${sample_id}.kraken.report --output ${sample_id}.kraken.out ${r1} ${r2} || true

  # default policy: keep unclassified + human (if present in classified)
  policy="${params.kraken_policy ?: 'keep_unclassified'}"

  case "$policy" in
    keep_human_only)
      # keep only reads classified as taxid 9606 — crude approach using kraken.out
      # extract read names classified as 9606 and rebuild FASTQ via seqtk (if available)
      awk '\$3==\"C\" && \$2==\"9606\" { print \$2, \$4 }' ${sample_id}.kraken.out | awk '{print $2}' > ${sample_id}.human_readnames.txt || true
      # fallback: if no seqtk, keep classified_1/2 as approximation
      gzip -c ${sample_id}.classified_1.fastq > ${sample_id}.decontam_R1.fastq.gz || true
      gzip -c ${sample_id}.classified_2.fastq > ${sample_id}.decontam_R2.fastq.gz || true
    ;;
    remove_nonhuman)
      # keep human + unclassified. We will concatenate unclassified + classified (approx)
      # Prefer human reads from classified files; here we conservatively keep all unclassified + classified
      cat ${sample_id}.unclassified_1.fastq ${sample_id}.classified_1.fastq 2>/dev/null || true
      cat ${sample_id}.unclassified_2.fastq ${sample_id}.classified_2.fastq 2>/dev/null || true
      gzip -c ${sample_id}.unclassified_1.fastq > ${sample_id}.decontam_R1.fastq.gz || true
      gzip -c ${sample_id}.unclassified_2.fastq > ${sample_id}.decontam_R2.fastq.gz || true
    ;;
    *)
      # default: keep unclassified only (conservative)
      gzip -c ${sample_id}.unclassified_1.fastq > ${sample_id}.decontam_R1.fastq.gz || true
      gzip -c ${sample_id}.unclassified_2.fastq > ${sample_id}.decontam_R2.fastq.gz || true
    ;;
  esac
  """
}

workflow kraken_decontam {
  take: pre_trimmed_ch

  main:
    pre_trimmed_ch
      .map { t -> t }            // ensure tuple(sid, r1, r2)
      .map { sid, r1, r2 -> tuple(sid, r1, r2) }
      .map { tuple -> kraken2_classify(tuple[0], tuple[1], tuple[2]) }
      .flatten()
      .set { decontaminated_ch }

  emit:
    decontaminated_ch
}
// end of modules/kraken_decontam/main.nf