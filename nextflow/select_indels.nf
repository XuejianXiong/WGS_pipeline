nextflow.enable.dsl=2

process SELECT_INDELS {
    publishDir params.outdir, mode: 'copy'
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
    tuple val(sample_id), path("filtered_indels.vcf.gz")


    script:
    """
    echo "=== [SELECT_INDELS] Starting INDEL selection and filtering for ${sample_id} ==="

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
    # Step 2: Select INDEL variants only
    # -------------------------------------------------
    echo "[STEP 1] Selecting INDEL variants..."
    gatk SelectVariants \\
        -R ${reference} \\
        -V ${raw_vcf} \\
        --select-type-to-include INDEL \\
        -O raw_indels.vcf.gz

    # -------------------------------------------------
    # Step 3: Index INDEL VCF
    # -------------------------------------------------
    tabix -f -p vcf raw_indels.vcf.gz

    # -------------------------------------------------
    # Step 4: Apply hard filters on INDELs
    # -------------------------------------------------
    echo "[STEP 2] Filtering INDELs..."
    gatk VariantFiltration \\
        -R ${reference} \\
        -V raw_indels.vcf.gz \\
        -O filtered_indels.vcf.gz \\
        -filter-name "QD_filter" -filter "QD < 2.0" \\
        -filter-name "FS_filter" -filter "FS > 200.0" \\
        -filter-name "SOR_filter" -filter "SOR > 10.0" \\
        -genotype-filter-expression "DP < 10" \\
        -genotype-filter-name "DP_filter" \\
        -genotype-filter-expression "GQ < 10" \\
        -genotype-filter-name "GQ_filter"

    echo "=== [SELECT_INDELS] Completed successfully for ${sample_id} ==="
    """
}
