// modules/postprocess/main.nf
nextflow.enable.dsl = 2

/*
 Postprocess module: collects per-sample caller outputs and runs tools/merge_fusions.py
 Input: channel of caller outputs (tuples or process outputs from fusion_callers)
 Output: per-sample merged TSV reports in params.outdir/reports
*/

process merge_and_annotate {
  tag { "merge:${sample_id}" }
  cpus 2
  time '2h'
  publishDir "${params.outdir ?: './results'}/reports", mode: 'copy', optional: true

  input:
    val sample_id
    path starfusion_file optional true
    path arriba_file optional true
    path fusioncatcher_file optional true

  output:
    path("${sample_id}.merged_fusions.tsv")

  script:
  """
  set -euo pipefail

  # default args for merge_fusions.py, will pass empty values for missing inputs
  STARF="${starfusion_file ?: ''}"
  ARR="${arriba_file ?: ''}"
  FC="${fusioncatcher_file ?: ''}"

  python3 ${workDir}/../tools/merge_fusions.py \
    --starfusion "${STARF}" \
    --arriba "${ARR}" \
    --fusioncatcher "${FC}" \
    --gtf "${params.gtf ?: ''}" \
    --out ${sample_id}.merged_fusions.tsv
  """
}

workflow postprocess {
  take: callers_ch

  /*
   callers_ch: expects a mixed stream of tuples produced by fusion_callers processes.
   For simplicity we will assume the pipeline publishes caller outputs into
   standard locations under ${params.outdir}/<sample>/fusion/<caller>/ and merge per sample by globbing.
  */
  main:
    callers_ch
      .map { it }    // passthrough
      .collect()     // gather emitted tuples (caller processes already wrote results to publishDir)
      .map { list ->
        // build a list of sample IDs by scanning results directory
        def outdir = params.outdir ?: './results'
        def samples = []
        new File(outdir).eachDir { d -> samples << d.name }
        samples
      }
      .flatten()
      .map { sid ->
        // find files produced by callers for this sample
        def base = file("${params.outdir ?: './results'}/${sid}")
        def starf = file("${base}/fusion/starfusion/star-fusion.fusion_predictions.tsv")
        if (!starf.exists()) starf = file("${base}/fusion/starfusion/starfusion.fusion_predictions.tsv")
        def arriba = file("${base}/fusion/arriba/arriba.tsv")
        def fusioncatcher = file("${base}/fusion/fusioncatcher/final-list_candidate-fusion-genes.txt")
        merge_and_annotate(sid, starf.exists() ? starf : null, arriba.exists() ? arriba : null, fusioncatcher.exists() ? fusioncatcher : null)
      }
      .flatten()
      .set { merged_reports }

  emit:
    merged_reports
}
/*
 STAR alignment module: runs STAR 2-pass alignment.

 Input: tuples (sample_id, r1, r2)  [these should be trimmed or decontaminated]
 Output: tuples (sample_id, sorted_bam, Chimeric.out.junction, Aligned.out.bam)
 Optionally emits metadata with original FASTQ paths (meta)
*/