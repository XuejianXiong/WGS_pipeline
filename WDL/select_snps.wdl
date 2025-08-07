version 1.0

task select_snps {
    input {
        File reference
        File reference_fai
        File reference_dict
        File raw_vcf
    }

    command {

        set -euo pipefail

        # Index for SelectVariants
        tabix -f -p vcf ${raw_vcf}

        # Separate SNPs 
        gatk SelectVariants \
            -R ${reference} \
            -V ${raw_vcf} \
            --select-type-to-include SNP \
            -O raw_snps.vcf.gz

        # Index for VariantFiltration
        tabix -f -p vcf raw_snps.vcf.gz

        # Apply Variant Quality Filtering
        gatk VariantFiltration \
            -R ${reference} \
            -V raw_snps.vcf.gz \
            -O filtered_snps.vcf.gz \
            -filter-name "QD_filter" -filter "QD < 2.0" \
            -filter-name "FS_filter" -filter "FS > 60.0" \
            -filter-name "MQ_filter" -filter "MQ < 40.0" \
            -filter-name "SOR_filter" -filter "SOR > 4.0" \
            -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
            -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
            -genotype-filter-expression "DP < 10" \
            -genotype-filter-name "DP_filter" \
            -genotype-filter-expression "GQ < 10" \
            -genotype-filter-name "GQ_filter"

    }

    output {
        File snps_vcf = "filtered_snps.vcf.gz"
    }

    runtime {
        cpu: 6
        memory: "8G"
        docker: "wgs-gatk"
    }
}
