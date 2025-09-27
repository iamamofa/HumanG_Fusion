// modules/fusion_callers/main.nf
nextflow.enable.dsl = 2

/*
 Fusion callers module: runs Arriba, STAR-Fusion, optionally FusionCatcher.

 Input: tuples (sample_id, sorted_bam, chimeric_junction, aligned_bam)
      OR: (sample_id, trimmed_r1, trimmed_r2, sorted_bam, chimeric_junction)
 Output: emits a channel of maps/tuples with caller outputs for downstream postprocessing.

 Requires params.arriba_blacklist (optional), params.genome_fasta, params.gtf, params.star_fusion_ctat_lib for STAR-Fusion.
*/

process run_arriba {
  tag { "arriba:${sample_id}" }
  cpus { params.threads ?: 4 }
  publishDir "${params.outdir ?: './results'}/\${sample_id}/fusion/arriba", mode: 'copy', optional: true

  input:
    tuple val(sample_id), path(chimeric), path(sorted_bam)

  output:
    tuple val(sample_id), path("arriba.tsv")

  script:
  """
  set -euo pipefail
  if [ -z "${params.genome_fasta}" ] || [ -z "${params.gtf}" ]; then
    echo "ERROR: params.genome_fasta and params.gtf must be set for Arriba" >&2
    exit 1
  fi

  # Arriba expects the chimeric junction file produced by STAR in BAM/chimeric format; typical usage:
  arriba -x ${chimeric} -b ${sorted_bam} -o arriba.tsv -O arriba.discarded.tsv \
    -a ${params.genome_fasta} -g ${params.gtf} ${ params.arriba_blacklist ? "-t ${params.arriba_blacklist}" : "" } || true
  """
}

process run_starfusion {
  tag { "starfusion:${sample_id}" }
  cpus { params.threads ?: 8 }
  publishDir "${params.outdir ?: './results'}/\${sample_id}/fusion/starfusion", mode: 'copy', optional: true

  /*
   STAR-Fusion can accept FASTQs or a genome lib dir (ctat resource).
   The simplest approach here is to attempt to run STAR-Fusion with FASTQ inputs if available,
   else run with the BAM (STAR-Fusion typically can re-run STAR internally — depends on version).
  */
  input:
    tuple val(sample_id), path(r1) optional true, path(r2) optional true, path(chimeric) optional true, path(sorted_bam) optional true

  output:
    tuple val(sample_id), path("starfusion.fusion_predictions.tsv") optional true

  script:
  """
  set -euo pipefail
  if [ -z "${params.star_fusion_ctat_lib}" ]; then
    echo "WARNING: params.star_fusion_ctat_lib not set. STAR-Fusion may fail without CTAT resource." >&2
  fi

  # Attempt to call STAR-Fusion with FASTQs if present, else fallback to running with --genome_lib_dir only.
  if [ -f "${r1:-}" ] && [ -f "${r2:-}" ]; then
    STAR-Fusion --left_fq ${r1} --right_fq ${r2} --genome_lib_dir ${params.star_fusion_ctat_lib} --CPU ${task.cpus} --output_dir starfusion_out || true
  else
    # try running STAR-Fusion with BAM/chimeric inputs if supported by the installed version (may not be supported)
    STAR-Fusion --genome_lib_dir ${params.star_fusion_ctat_lib} --CPU ${task.cpus} --output_dir starfusion_out || true
  fi

  # copy expected result
  if [ -d starfusion_out ]; then
    cp -r starfusion_out/* . || true
  fi
  """
}

process run_fusioncatcher {
  tag { "fusioncatcher:${sample_id}" }
  cpus { params.threads ?: 8 }
  when: params.run_fusioncatcher
  publishDir "${params.outdir ?: './results'}/\${sample_id}/fusion/fusioncatcher", mode: 'copy', optional: true

  input:
    tuple val(sample_id), path(r1), path(r2)

  output:
    tuple val(sample_id), path("fusioncatcher_final-list_candidate-fusion-genes.txt") optional true

  script:
  """
  set -euo pipefail
  fusioncatcher -d . -p ${task.cpus} -i ${r1},${r2} || true
  # copy or move results to workspace (names depend on fusioncatcher version)
  if [ -f final-list_candidate-fusion-genes.txt ]; then
    cp final-list_candidate-fusion-genes.txt fusioncatcher_final-list_candidate-fusion-genes.txt || true
  fi
  """
}

workflow fusion_callers {
  take: aligned_ch

  main:
    /*
     aligned_ch is expected to carry tuples. It might contain different shapes depending on upstream wiring.
     We'll normalize expected shapes using a small mapping.
    */

    // map aligned channel to call processes
    aligned_ch.map { tup ->
      // allowed shapes:
      // 1) (sample_id, sorted.bam, Chimeric.out.junction, Aligned.out.bam)
      // 2) (sample_id, trimmed_R1, trimmed_R2, sorted.bam, Chimeric.out.junction) - if metadata passed
      tup
    }
    .subscribe { /* debug hook */ }

    // launch Arriba and STAR-Fusion. We'll try to extract fields defensively.
    caller_ch = Channel.create()

    aligned_ch.subscribe { t ->
      def sid = t[0]
      def paths = t[1..-1]
      // find files by suffix if possible
      def sorted_bam = paths.find { it.name.endsWith('.sorted.bam') || it.name.endsWith('.Aligned.out.bam') || it.name.endsWith('.bam') }
      def chimeric = paths.find { it.name.contains('Chimeric.out') || it.name.endsWith('.junction') || it.name.endsWith('.Chimeric.out.junction') }
      def r1 = paths.find { it.name.toLowerCase().contains('_r1') || it.name.toLowerCase().endsWith('.fastq.gz') } ?: null
      def r2 = paths.find { it.name.toLowerCase().contains('_r2') || it.name.toLowerCase().endsWith('.fastq.gz') } ?: null

      // run arriba if we have chimeric+sorted
      if (chimeric && sorted_bam) {
        caller_ch << run_arriba(sid, chimeric, sorted_bam)
      }

      // run starfusion with FASTQs if present, else attempt with BAM
      caller_ch << run_starfusion(sid, r1, r2, chimeric, sorted_bam)

      // optional fusioncatcher (requires FASTQs)
      if (params.run_fusioncatcher && r1 && r2) {
        caller_ch << run_fusioncatcher(sid, r1, r2)
      }
    }

    // flatten and make a combined output channel of maps per sample with paths
    // The downstream postprocess expects per-sample files; we will rely on publishDir to place caller outputs
    emit:
      caller_ch
}
/*
 STAR alignment module: runs STAR 2-pass alignment.

 Input: tuples (sample_id, r1, r2)  [these should be trimmed or decontaminated]
 Output: tuples (sample_id, sorted_bam, Chimeric.out.junction, Aligned.out.bam)
 Optionally emits metadata with original FASTQ paths (meta)
*/