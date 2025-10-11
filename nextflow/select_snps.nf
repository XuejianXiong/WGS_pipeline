nextflow.enable.dsl=2

process SELECT_SNPS {
    label 'gatk'
    tag "$sample_id"
    cpus 6
    memory '8 GB'
    maxForks 1

    input:
    tuple val(sample_id), 
          path(reference),
          path(reference_fai),
          path(reference_dict),
          path(raw_vcf)

    output:
    tuple val(sample_id), path("filtered_snps.vcf.gz")


    publishDir params.outdir, mode: 'copy'

    script:
    """
    echo "=== [SELECT_SNPS] Starting SNP selection and filtering for ${sample_id} ==="

    # -------------------------------------------------
    # Step 1: Index input VCF if not already indexed
    # -------------------------------------------------
    if [[ ! -f ${raw_vcf}.tbi && ! -f ${raw_vcf}.csi ]]; then
        echo "[INFO] Indexing raw VCF: ${raw_vcf}"
        tabix -f -p vcf ${raw_vcf}
    else
        echo "[INFO] Raw VCF already indexed"
    fi

    # -------------------------------------------------
    # Step 2: Select SNP variants only
    # -------------------------------------------------
    echo "[STEP 1] Selecting SNP variants..."
    gatk SelectVariants \\
        -R ${reference} \\
        -V ${raw_vcf} \\
        --select-type-to-include SNP \\
        -O raw_snps.vcf.gz

    # -------------------------------------------------
    # Step 3: Index SNP VCF
    # -------------------------------------------------
    tabix -f -p vcf raw_snps.vcf.gz

    # -------------------------------------------------
    # Step 4: Apply hard filters on SNPs
    # -------------------------------------------------
    echo "[STEP 2] Filtering SNPs..."
    gatk VariantFiltration \\
        -R ${reference} \\
        -V raw_snps.vcf.gz \\
        -O filtered_snps.vcf.gz \\
        -filter-name "QD_filter" -filter "QD < 2.0" \\
        -filter-name "FS_filter" -filter "FS > 60.0" \\
        -filter-name "MQ_filter" -filter "MQ < 40.0" \\
        -filter-name "SOR_filter" -filter "SOR > 4.0" \\
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \\
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \\
        -genotype-filter-expression "DP < 10" \\
        -genotype-filter-name "DP_filter" \\
        -genotype-filter-expression "GQ < 10" \\
        -genotype-filter-name "GQ_filter"

    echo "=== [SELECT_SNPS] Completed successfully for ${sample_id} ==="
    """
}
