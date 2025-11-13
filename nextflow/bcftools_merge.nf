nextflow.enable.dsl=2

process BCFTOOLS_MERGE {
    publishDir params.outdir, mode: 'copy'
    label 'bcftools'
    tag "$sample_id"
    cpus 6
    memory '8 GB'
    maxForks 1

    input:
    tuple val(sample_id), path(filtered_snps), path(filtered_indels)

    output:
    tuple val(sample_id), path("merged_variants.vcf")


    script:
    """
    echo "=== [BCFTOOLS_MERGE] Starting merge of SNPs and INDELs for ${sample_id} ==="

    # -------------------------------------------------
    # Index the input VCFs
    # -------------------------------------------------
    if [[ ! -f ${filtered_snps}.tbi && ! -f ${filtered_snps}.csi ]]; then
        tabix -f -p vcf ${filtered_snps}
    fi

    if [[ ! -f ${filtered_indels}.tbi && ! -f ${filtered_indels}.csi ]]; then
        tabix -f -p vcf ${filtered_indels}
    fi

    # -------------------------------------------------
    # Merge VCFs using bcftools
    # -------------------------------------------------
    bcftools concat -a -O v \\
        ${filtered_snps} \\
        ${filtered_indels} \\
        -o merged_variants.vcf

    echo "=== [BCFTOOLS_MERGE] Completed successfully for ${sample_id} ==="
    """
}
