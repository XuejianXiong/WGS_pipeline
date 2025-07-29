version 1.0

task bqsr {
    input {
        File reference
        File reference_fai
        File reference_dict
        File known_variants_snps
        File known_variants_indels
        File dedup_bam
    }

    command {

        # Compress and index known variants if needed
        # (it's normally finished in the setup step)
        #bgzip -c ${known_variants_snps} > known_sites.vcf.gz        
        #tabix -p vcf known_sites.vcf.gz

        # Index known variants
        tabix -p vcf ${known_variants_snps}
        tabix -p vcf ${known_variants_indels}
        

        # Index the BAM file
        samtools index ${dedup_bam}        

        # Base recalibration
        gatk BaseRecalibrator \
            -I ${dedup_bam} \
            -R ${reference} \
            --known-sites ${known_variants_snps} \
            --known-sites ${known_variants_indels} \
            -O bqsr_data.table

        # Apply recalibration
        gatk ApplyBQSR \
            -I ${dedup_bam} \
            -R ${reference} \
            --bqsr-recal-file bqsr_data.table \
            -O bqsr.bam
    }

    output {
        File bqsr_bam = "bqsr.bam"
        File bqsr_report = "bqsr_data.table"
    }

    runtime {
        cpu: 3
        memory: "4G"
        docker: "wgs-gatk"
    }
}
