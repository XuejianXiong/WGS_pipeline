nextflow.enable.dsl=2

// -------------------------
// FASTQC process
// -------------------------
process FASTQC {
    publishDir "${params.outdir}", mode: 'copy'
    label 'qc'
    tag "$sample_id"
    cpus 2
    memory '4 GB'
    // container 'wgs-qc-reads:latest'
    maxForks 2

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id),
          path("${fastq1.baseName}_fastqc.zip"),
          path("${fastq1.baseName}_fastqc.html"),
          path("${fastq2.baseName}_fastqc.zip"),
          path("${fastq2.baseName}_fastqc.html")

    script:
    """
    fastqc -o . $fastq1 $fastq2
    """
}


// -------------------------
// MULTIQC process
// -------------------------
process MULTIQC {
    publishDir "${params.outdir}", mode: 'copy'
    label 'qc'
    tag "multiqc"
    cpus 1
    memory '2 GB'
    // container 'wgs-qc-reads:latest'

    input:
    path fastqc_files

    output:
    path "multiqc_report.html"
    path "multiqc_data"


    script:
    """
    multiqc $fastqc_files -o .
    """
}
