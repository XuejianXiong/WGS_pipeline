version 1.0

task select_indels {
    input {
        File reference
        File reference_fai
        File reference_dict
        File raw_vcf
    }

    command {

        # Index for SelectVariants
        tabix -f -p vcf ${raw_vcf}

        # Separate INDELs 
        gatk SelectVariants \
            -R ${reference} \
            -V ${raw_vcf} \
            --select-type-to-include INDEL \
            -O raw_indels.vcf.gz

        # Index for VariantFiltration
        tabix -f -p vcf raw_indels.vcf.gz

        # Flag low-quality variants
        gatk VariantFiltration \
            -R ${reference} \
            -V raw_indels.vcf.gz \
            -O filtered_indels.vcf.gz \
            --filter-name "QD_filter" -filter "QD < 2.0" \
            -filter-name "FS_filter" -filter "FS > 200.0" \
            -filter-name "SOR_filter" -filter "SOR > 10.0" \
            -genotype-filter-expression "DP < 10" \
            -genotype-filter-name "DP_filter" \
            -genotype-filter-expression "GQ < 10" \
            -genotype-filter-name "GQ_filter"

    }

    output {
        File indels_vcf = "filtered_indels.vcf.gz"
    }

    runtime {
        cpu: 6
        memory: "8G"
        docker: "wgs-gatk"
    }
}
