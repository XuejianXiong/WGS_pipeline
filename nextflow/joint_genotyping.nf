nextflow.enable.dsl=2

process JOINT_GENOTYPING {
    label 'gatk'
    tag "joint_genotyping"
    cpus 6
    memory '8 GB'
    maxForks 1

    input:
    tuple path(reference), path(reference_fai), path(reference_dict), path(gvcf_1), path(gvcf_2), path(gvcf_3)

    output:
    tuple val('joint'), path("raw_variants.vcf.gz")

    publishDir params.outdir, mode: 'copy'

    script:
    """
    echo "=== [JOINT_GENOTYPING] Starting joint genotyping for trio ==="

    # Index the GVCFs
    gatk IndexFeatureFile -I ${gvcf_1}
    gatk IndexFeatureFile -I ${gvcf_2}
    gatk IndexFeatureFile -I ${gvcf_3}

    # Combine the trio GVCFs
    gatk CombineGVCFs \\
        -R ${reference} \\
        --variant ${gvcf_1} \\
        --variant ${gvcf_2} \\
        --variant ${gvcf_3} \\
        -O combined.g.vcf.gz

    # Joint Genotyping
    gatk GenotypeGVCFs \\
        -R ${reference} \\
        -V combined.g.vcf.gz \\
        -O raw_variants.vcf.gz

    echo "=== [JOINT_GENOTYPING] Completed joint genotyping ==="
    """
}
