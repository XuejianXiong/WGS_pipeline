nextflow.enable.dsl=2

process QC_GATK {
    label 'gatk' 
    tag "$sample_id"
    cpus 3
    memory '4 GB'
    maxForks 2

    input:
    tuple val(sample_id), path(dedup_bam), path(reference)

    output:
    tuple val(sample_id), 
          path("${sample_id}.alignment_metrics.txt"),
          path("${sample_id}.insert_size_metrics.txt"),
          path("${sample_id}.insert_size_histogram.pdf")

    publishDir params.outdir, mode: 'copy'

    script:
    """
    gatk CollectAlignmentSummaryMetrics \
        -R $reference \
        -I $dedup_bam \
        -O ${sample_id}.alignment_metrics.txt

    gatk CollectInsertSizeMetrics \
        -I $dedup_bam \
        -O ${sample_id}.insert_size_metrics.txt \
        -H ${sample_id}.insert_size_histogram.pdf
    """
}
