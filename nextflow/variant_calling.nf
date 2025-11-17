nextflow.enable.dsl=2

process VARIANT_CALLING {
    publishDir params.outdir, mode: 'copy'
    label 'gatk'
    tag "$sample_id"
    cpus 6
    memory '8 GB'
    maxForks 2


    input:
    tuple val(sample_id), path(reference), path(reference_fai), path(reference_dict), path(bqsr_bam)

    output:
    tuple val(sample_id), path("${sample_id}.raw_variants.g.vcf.gz")


    script:
    """
    echo "=== [VARIANT_CALLING] Starting HaplotypeCaller for ${sample_id} ==="

    # Index BAM (required by GATK)
    samtools index ${bqsr_bam}

    # Run GATK HaplotypeCaller to generate GVCF
    gatk HaplotypeCaller \\
        -R ${reference} \\
        -I ${bqsr_bam} \\
        -O ${sample_id}.raw_variants.g.vcf.gz \\
        -ERC GVCF \\
        --native-pair-hmm-threads ${task.cpus} \\
        -L chr22

    echo "=== [VARIANT_CALLING] Completed for ${sample_id} ==="
    """
}
