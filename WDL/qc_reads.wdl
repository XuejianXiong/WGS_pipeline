version 1.0

task qc_fastqc {
    input {
        File fastq1
        File fastq2
        String qc_dir
    }

    command <<<
        set -euo pipefail

        echo "qc_dir is: ~{qc_dir}"

        mkdir -p "~{qc_dir}"

        echo "📊 Running FastQC..."
        #export JAVA_OPTS="-Xmx1g"
        fastqc -o "~{qc_dir}" ~{fastq1} ~{fastq2}

        mv "~{qc_dir}/$(basename ~{fastq1} .fastq)_fastqc.zip" fastqc1.zip
        mv "~{qc_dir}/$(basename ~{fastq2} .fastq)_fastqc.zip" fastqc2.zip
    >>>

    output {
        File fastqc_output1 = "fastqc1.zip"
        File fastqc_output2 = "fastqc2.zip"
    }

    runtime {
        cpu: 6
        memory: "8G"
        docker: "wgs-qc_reads"
    }
}

task multiqc_summary {
    input {
        String qc_dir
    }

    command <<<
        echo "📊 Running MultiQC..."
        multiqc "~{qc_dir}" -o "~{qc_dir}"
    >>>

    output {
        File multiqc_report = "~{qc_dir}/multiqc_report.html"
    }

    runtime {
        cpu: 6
        memory: "8G"
        docker: "wgs-qc_reads"
    }
}
