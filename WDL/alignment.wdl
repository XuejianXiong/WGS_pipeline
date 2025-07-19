task alignment {
    input {
        File reference
        File fastq1
        File fastq2
    }

    command {
        bwa index ${reference}
        bwa mem -t 6 -R "@RG\tID:sample1\tSM:sample1\tPL:illumina" ${reference} ${fastq1} ${fastq2} > aligned.sam
    }

    output {
        File sam = "aligned.sam"  
    }

    runtime {
        cpu: 6
        memory: "8G"
        docker: "wgs-align"
    }
}