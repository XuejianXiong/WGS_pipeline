version 1.0

task bqsr {
    input {
        File reference
        File reference_fai
        File reference_dict
        File known_variants
        File dedup_bam
    }

    command <<<

        # Compress and index known variants
        bgzip -c ${known_variants} > known_sites.vcf.gz        
        tabix -p vcf known_sites.vcf.gz
        
        # Index the BAM file
        samtools index ${dedup_bam}        

        # Base recalibration
        gatk BaseRecalibrator \
            -I ${dedup_bam} \
            -R ${reference} \
            --known-sites known_sites.vcf.gz \
            -O bqsr_data.table

        # Apply recalibration
        gatk ApplyBQSR \
            -I ${dedup_bam} \
            -R ${reference} \
            --bqsr-recal-file bqsr_data.table \
            -O bqsr.bam
    >>>

    output {
        File bqsr_bam = "bqsr.bam"
        File bqsr_report = "bqsr_data.table"
    }

    runtime {
        cpu: 6
        memory: "8G"
        docker: "wgs-gatk"
    }
}
