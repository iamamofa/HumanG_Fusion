// modules/preprocess/main.nf
nextflow.enable.dsl = 2

include { fastqc_run; cutadapt_trim } from '../utils/processes.nf'

/*
 Preprocess workflow:
  - Input: Channel of tuples (sample_id, [r1, r2]) OR Channel.fromFilePairs glob
  - Output: Channel of tuples (sample_id, trimmed_r1, trimmed_r2)
*/

workflow preprocess {
  take: samples_ch

  main:
    samples_ch
      .map { it }               // ensure tuple(sid, [r1,r2]) format
      .map { sid, reads -> tuple(sid, reads[0], reads[1]) }
      .set { _reads_pairs }

    // Run pre-trim FastQC (best-effort)
    _reads_pairs
      .map { sid, r1, r2 -> tuple(sid, r1, r2) }
      .ifEmpty { channel.empty() }
      .map { sid, r1, r2 -> tuple(sid, r1, r2) }
      .subscribe { /* for debugging */ }

    // run fastqc (non-blocking; errors suppressed)
    pre_qc = _reads_pairs.map { sid, r1, r2 -> tuple(sid, r1, r2) }.into { a -> a }
    fastqc_run(pre_qc)

    // trimming
    trimmed = _reads_pairs.map { sid, r1, r2 -> tuple(sid, r1, r2) }.map { t -> cutadapt_trim(t) }.flatten()

  emit:
    trimmed
}
// End of modules/preprocess/main.nf