version 1.0

task qc_gatk {
    input {
        File reference
        File dedup_bam
    }

    command {
        gatk CollectAlignmentSummaryMetrics -R ${reference} -I ${dedup_bam} -O alignment_metrics.txt

        gatk CollectInsertSizeMetrics -I ${dedup_bam} -O insert_size_metrics.txt -H insert_size_histogram.pdf
    }

    output {
        File align_txt = "alignment_metrics.txt"
        File insert_txt = "insert_size_metrics.txt"
        File hist_pdf = "insert_size_histogram.pdf"
    }

    runtime {
        cpu: 3
        memory: "4G"
        docker: "wgs-gatk"
    }
}