nextflow.enable.dsl=2

// -------------------------
// FASTP trimming process
// -------------------------
process TRIM_FASTQ {
    publishDir "${params.outdir}", mode: 'copy'
    label 'trim'
    tag "$sample_id"
    cpus 4
    memory '6 GB'
    // container 'wgs-fastp'
    maxForks 2

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id),
          path("${sample_id}_1.trimmed.fastq"),
          path("${sample_id}_2.trimmed.fastq"),
          path("${sample_id}_fastp_report.html"),
          path("${sample_id}_fastp_report.json")


    script:
    """
    fastp \
      -i $fastq1 \
      -I $fastq2 \
      -o ${sample_id}_1.trimmed.fastq \
      -O ${sample_id}_2.trimmed.fastq \
      --detect_adapter_for_pe \
      --thread ${task.cpus} \
      --html ${sample_id}_fastp_report.html \
      --json ${sample_id}_fastp_report.json \
      --report_title "fastp report for $sample_id"
    """
}
