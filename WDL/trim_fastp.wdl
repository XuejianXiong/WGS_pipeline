version 1.0

task trim_fastp {
    input {
        String sample_id
        File fastq1
        File fastq2
        Int threads = 4
    }

    command <<<
        set -euo pipefail

        echo "Running fastp for sample ~{sample_id}"

        fastp \
            -i ~{fastq1} \
            -I ~{fastq2} \
            -o ~{sample_id}_1.trimmed.fastq \
            -O ~{sample_id}_2.trimmed.fastq \
            --detect_adapter_for_pe \
            --thread ~{threads} \
            --html ~{sample_id}_fastp_report.html \
            --json ~{sample_id}_fastp_report.json \
            --report_title "fastp report for ~{sample_id}"

        echo "Creating output archive..."

        tar -czf ~{sample_id}_fastp_outputs.tar.gz \
            ~{sample_id}_1.trimmed.fastq \
            ~{sample_id}_2.trimmed.fastq \
            ~{sample_id}_fastp_report.html \
            ~{sample_id}_fastp_report.json
    >>>

    output {
        File output_tar = "~{sample_id}_fastp_outputs.tar.gz"
    }

    runtime {
        cpu: 6
        memory: "8G"
        docker: "wgs-fastp"
    }
}
