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
    # Index reference (only runs if not already present)
    if [ ! -f "${reference}.bwt" ]; then
        bwa index $reference
    fi

    # Align with read group header for GATK compatibility
    bwa mem -t ${task.cpus} -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \
        $reference $fastq1 $fastq2 \
        | samtools view -Sb - \
        | samtools sort -o ${sample_id}.sorted.bam

    # Index the BAM file
    samtools index ${sample_id}.sorted.bam
    """
}
