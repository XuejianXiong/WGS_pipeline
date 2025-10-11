nextflow.enable.dsl=2

process QC_VARIANT {
    label 'bcftools'
    tag "${sample_id}"
    cpus 6
    memory '8 GB'
    maxForks 1

    input:
    tuple val(sample_id), path(vcf_file)

    output:
    tuple val(sample_id),
          path("variant_stats.txt"),
          path("qc_plots.zip")

    publishDir params.outdir, mode: 'copy'

    script:
    """
    echo "=== [QC_VARIANT] Starting variant QC for sample: ${sample_id} ==="

    # -------------------------------------------------
    # Compress and index VCF
    # -------------------------------------------------
    bgzip -f ${vcf_file}
    tabix -f -p vcf ${vcf_file}.gz

    # -------------------------------------------------
    # Compute variant statistics
    # -------------------------------------------------
    bcftools stats ${vcf_file}.gz > variant_stats.txt

    # -------------------------------------------------
    # Generate plots
    # -------------------------------------------------
    mkdir -p qc_plots
    plot-vcfstats -p qc_plots variant_stats.txt

    # -------------------------------------------------
    # Zip the plots
    # -------------------------------------------------
    zip -r qc_plots.zip qc_plots

    echo "=== [QC_VARIANT] Completed successfully for ${sample_id} ==="
    """
}
