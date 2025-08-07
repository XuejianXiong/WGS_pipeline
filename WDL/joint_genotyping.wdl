version 1.0

task joint_genotyping {
    input {
        File reference
        File reference_fai
        File reference_dict
        File gvcf_1
        File gvcf_2
        File gvcf_3
    }

    command {

        # Index the GVCFs
        gatk IndexFeatureFile -I ${gvcf_1}
        gatk IndexFeatureFile -I ${gvcf_2}
        gatk IndexFeatureFile -I ${gvcf_3}

        gatk CombineGVCFs \
            -R ${reference} \
            --variant ${gvcf_1}\
            --variant ${gvcf_2}\
            --variant ${gvcf_3}\
            -O combined.g.vcf.gz

        # Joint Genotyping
        gatk GenotypeGVCFs \
            -R ${reference} \
            -V combined.g.vcf.gz \
            -O raw_variants.vcf.gz
    }

    output {
        File raw_vcf = "raw_variants.vcf.gz"
    }

    runtime {
        cpu: 6
        memory: "8G"
        docker: "wgs-gatk"
    }
}
