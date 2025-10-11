nextflow.enable.dsl=2

process DEDUP {
    label 'gatk'
    tag "$sample_id"
    cpus 4
    memory '6 GB'
    maxForks 2

    input:
    tuple val(sample_id), path(sorted_bam), path(sorted_bai)

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.metrics.txt")


    publishDir params.outdir, mode: 'copy'

    script:
    """
    gatk MarkDuplicates \
        -I $sorted_bam \
        -O ${sample_id}.dedup.bam \
        -M ${sample_id}.dedup.metrics.txt \
        --VALIDATION_STRINGENCY LENIENT
    """
}
