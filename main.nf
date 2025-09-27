#!/usr/bin/env nextflow
/*
 * RNA-Seq fusion detection pipeline
 */

nextflow.enable.dsl=2

include { preprocess }      from './modules/preprocess/main.nf'
include { kraken_decontam } from './modules/kraken_decontam/main.nf'
include { star_align }      from './modules/star_align/main.nf'
include { fusion_callers }  from './modules/fusion_callers/main.nf'
include { postprocess }     from './modules/postprocess/main.nf'

// --- SAMPLE INPUTS ---
workflow {
    samples_ch = params.samples_tsv ? \
        Channel.fromPath(params.samples_tsv).splitCsv(header:true, sep:'\t').map{ row -> tuple(row.sample_id, [file(row.fastq_r1), file(row.fastq_r2)]) } : \
        Channel.fromFilePairs(params.reads ?: './sample_data/*_R1.fastq.gz', flat:true)

    processed      = preprocess(samples_ch)
    decontaminated = params.kraken2_db ? kraken_decontam(processed) : processed
    aligned        = star_align(decontaminated)
    callers        = fusion_callers(aligned)
    merged         = postprocess(callers)
}
