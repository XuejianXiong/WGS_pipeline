nextflow.enable.dsl=2

process BWA_ALIGN {
    label 'align'
    tag "$sample_id"
    cpus 4
    memory '6 GB'
    // container 'wgs-align:latest'
    maxForks 2

    input:
    tuple val(sample_id), path(fastq1), path(fastq2), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai")

    // ✅ keep publishDir exactly as before
    publishDir params.outdir, mode: 'copy'

    script:
    """
    bwa index $reference
    bwa mem -t ${task.cpus} $reference $fastq1 $fastq2 \\
        | samtools view -Sb - \\
        | samtools sort -o ${sample_id}.sorted.bam
    samtools index ${sample_id}.sorted.bam
    """
}
