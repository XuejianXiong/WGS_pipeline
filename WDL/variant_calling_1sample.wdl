version 1.0

task variant_calling {
    input {
        File reference
        File reference_fai
        File reference_dict
        File bqsr_bam
    }

    command {

        # Index the sorted BAM (recommended for GATK)
        samtools index ${bqsr_bam}

        # Run GATK HaplotypeCaller
        gatk HaplotypeCaller \
        -R ${reference} \
        -I ${bqsr_bam} \
        -O raw_variants.vcf \
        --native-pair-hmm-threads 6 \
        -L chr22

    }

    output {
        File vcf = "raw_variants.vcf"
    }

    runtime {
        cpu: 6
        memory: "8G"
        docker: "wgs-gatk"
    }
}