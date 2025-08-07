version 1.0

task bcftools_merge {
    input {
        File filtered_snps
        File filtered_indels
    }

    command {
        set -euo pipefail

        # Index the compressed VCFs
        tabix -f -p vcf ${filtered_snps}
        tabix -f -p vcf ${filtered_indels}

        # Merge SNPs and INDELs
        bcftools concat -a -O v \
            ${filtered_snps} \
            ${filtered_indels} \
            -o merged_variants.vcf
    }

    output {
        File merged_vcf = "merged_variants.vcf"
    }

    runtime {
        cpu: 6
        memory: "8G"
        docker: "wgs-bcftools"
    }
}
