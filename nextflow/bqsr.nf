nextflow.enable.dsl=2

process BQSR {
    label 'gatk'
    tag "$sample_id"
    cpus 3
    memory '4 GB'
    maxForks 2
    // container 'wgs-gatk:latest'

    input:
    tuple val(sample_id),
          path(reference),
          path(reference_fai),
          path(reference_dict),
          path(known_variants_snps),
          path(known_variants_indels),
          path(dedup_bam)

    output:
    tuple val(sample_id),
          path("${sample_id}.bqsr.bam"),
          path("${sample_id}.bqsr_data.table")

    publishDir params.outdir, mode: 'copy'

    script:
    """
    echo "=== [BQSR] Starting Base Quality Score Recalibration for ${sample_id} ==="

    # -------------------------------------
    # Ensure VCFs are compressed and indexed
    # -------------------------------------
    if [[ ! -f ${known_variants_snps}.tbi && ! -f ${known_variants_snps}.csi ]]; then
        echo "[INFO] Indexing SNP VCF: ${known_variants_snps}"
        tabix -p vcf ${known_variants_snps}
    else
        echo "[INFO] SNP VCF already indexed"
    fi

    if [[ ! -f ${known_variants_indels}.tbi && ! -f ${known_variants_indels}.csi ]]; then
        echo "[INFO] Indexing INDEL VCF: ${known_variants_indels}"
        tabix -p vcf ${known_variants_indels}
    else
        echo "[INFO] INDEL VCF already indexed"
    fi

    # Index BAM (required by GATK)
    samtools index ${dedup_bam}
    
    # -------------------------------------
    # Step 1: Base Recalibration Table
    # -------------------------------------
    echo "[STEP 1] Running BaseRecalibrator..."
    gatk BaseRecalibrator \\
        -I ${dedup_bam} \\
        -R ${reference} \\
        --known-sites ${known_variants_snps} \\
        --known-sites ${known_variants_indels} \\
        -O ${sample_id}.bqsr_data.table

    # -------------------------------------
    # Step 2: Apply BQSR Corrections
    # -------------------------------------
    echo "[STEP 2] Applying BQSR recalibration..."
    gatk ApplyBQSR \\
        -I ${dedup_bam} \\
        -R ${reference} \\
        --bqsr-recal-file ${sample_id}.bqsr_data.table \\
        -O ${sample_id}.bqsr.bam

    echo "=== [BQSR] Completed for ${sample_id} ==="
    """
}
