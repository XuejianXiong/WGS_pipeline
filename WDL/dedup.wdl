version 1.0

task dedup {
    input {
        File aligned_sam
    }

    command {
        gatk SortSam -I ${aligned_sam} -O sorted.bam --SORT_ORDER coordinate

        gatk MarkDuplicates -I sorted.bam -O dedup.bam -M dedup.metrics.txt
    }

    output {
        File bam = "dedup.bam"
        File dedup_report = "dedup.metrics.txt"
    }

    runtime {
        cpu: 6
        memory: "8G"
        docker: "wgs-gatk"
    }
}