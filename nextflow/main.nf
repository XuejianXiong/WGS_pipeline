nextflow.enable.dsl=2

include { FASTQC; MULTIQC } from './qc_reads.nf'
include { TRIM_FASTQ }      from './trim_fastp.nf'
include { BWA_ALIGN }       from './alignment.nf'

// -------------------------
// Load parameters
// -------------------------
params.samples   = params.samples
params.reference = params.reference
params.outdir    = params.outdir

// -------------------------
// Create channel of tuples
// -------------------------
samples_ch = Channel.fromList(params.samples.entrySet() as List)
    .map { entry -> tuple(entry.key, file(entry.value[0]), file(entry.value[1]), file(params.reference)) }

// -------------------------
// Main workflow
// -------------------------
workflow {

    // Prepare input for FASTQC (drop reference)
    qc_input = samples_ch.map { sid, fq1, fq2, ref -> tuple(sid, fq1, fq2) }

    // Run FASTQC
    fastqc_results = FASTQC(qc_input)

    // Collect FASTQC outputs
    fastqc_files = fastqc_results
        .map { sid, zip1, html1, zip2, html2 -> [zip1, html1, zip2, html2] }
        .flatten()
        .collect()

    // Run MultiQC
    multiqc_report = MULTIQC(fastqc_files)

    // Run TRIM_FASTQ
    trimmed_reads = TRIM_FASTQ(qc_input)

    // Send trimmed reads to alignment
    bwa_input = trimmed_reads.map { sid, fq1, fq2, html, json -> tuple(sid, fq1, fq2, file(params.reference)) }

    aligned_bams = BWA_ALIGN(bwa_input)
}
