nextflow.enable.dsl=2

// ------------------------------------------------------------
// Import process modules
// ------------------------------------------------------------
include { FASTQC; MULTIQC }  from './qc_reads.nf'
include { TRIM_FASTQ }       from './trim_fastp.nf'
include { BWA_ALIGN }        from './alignment.nf'
include { DEDUP }            from './dedup.nf'
include { QC_GATK }          from './qc_gatk.nf'
include { BQSR }             from './bqsr.nf'
include { VARIANT_CALLING }  from './variant_calling.nf'
include { JOINT_GENOTYPING } from './joint_genotyping.nf'
include { SELECT_SNPS }      from './select_snps.nf'
include { SELECT_INDELS }    from './select_indels.nf'
include { BCFTOOLS_MERGE }   from './bcftools_merge.nf'
include { QC_VARIANT }       from './qc_variant.nf'


// ------------------------------------------------------------
// Load parameters (samples, reference genome, output directory, etc.)
// ------------------------------------------------------------
params.samples   = params.samples
params.reference = params.reference
params.outdir    = params.outdir
params.known_variants_snps   = params.known_variants_snps
params.known_variants_indels = params.known_variants_indels


// ------------------------------------------------------------
// Main workflow definition
// ------------------------------------------------------------
workflow {

    // ------------------------------------------------------------
    // Build channel of input tuples
    //   -> each tuple = (sample_id, fastq1, fastq2, reference)
    // ------------------------------------------------------------
    samples_ch = Channel.fromList(params.samples.entrySet() as List)
        .map { entry -> 
            tuple(entry.key, file(entry.value[0]), file(entry.value[1]), file(params.reference))
        }

    // --------------------------------------------------------
    // FASTQC: quality check raw FASTQs
    // --------------------------------------------------------
    qc_input = samples_ch.map { sid, fq1, fq2, ref -> tuple(sid, fq1, fq2) }
    fastqc_results = FASTQC(qc_input)

    // Collect FASTQC output files for MultiQC
    fastqc_files = fastqc_results
        .map { sid, zip1, html1, zip2, html2 -> [zip1, html1, zip2, html2] }
        .flatten()
        .collect()

    // Run MultiQC (aggregate QC reports)
    multiqc_report = MULTIQC(fastqc_files)


    // --------------------------------------------------------
    // FASTP: trimming (adapter removal, quality filtering)
    // --------------------------------------------------------
    trimmed_reads = TRIM_FASTQ(qc_input)


    // --------------------------------------------------------
    // BWA: align trimmed reads → sorted BAM
    // --------------------------------------------------------
    bwa_input = trimmed_reads.map { sid, fq1, fq2, html, json -> 
        tuple(sid, fq1, fq2, file(params.reference)) 
    }
    
    aligned_bams = BWA_ALIGN(bwa_input)


    // --------------------------------------------------------
    // GATK MarkDuplicates: deduplication
    //   input  -> sorted.bam (from BWA step)
    //   output -> dedup.bam + metrics report
    // --------------------------------------------------------
    dedup_bams = DEDUP(aligned_bams)


    // --------------------------------------------------------
    // GATK QC: CollectAlignmentSummaryMetrics + InsertSizeMetrics
    // --------------------------------------------------------
    qcgatk_input = dedup_bams.map { sid, bam, metrics ->
        tuple(sid, bam, file(params.reference))
    }

    qcgatk_results = QC_GATK(qcgatk_input)


    // --------------------------------------------------------
    // BQSR: Base Quality Score Recalibration
    // --------------------------------------------------------
    bqsr_input = dedup_bams.map { sid, bam, metrics ->
        tuple(
            sid,            
            file(params.reference),
            file(params.reference + '.fai'),
            file(params.reference.replace('.fasta', '').replace('.fa', '') + '.dict'),
            file(params.known_variants_snps),
            file(params.known_variants_indels),
            bam
        )
    }

    bqsr_results = BQSR(bqsr_input)


    // --------------------------------------------------------
    // VARIANT CALLING: GATK HaplotypeCaller → per-sample gVCF
    // --------------------------------------------------------
    variant_calling_input = bqsr_results.map { sid, bqsr_bam, bqsr_table ->
        tuple(
            sid,            
            file(params.reference),
            file(params.reference + '.fai'),
            file(params.reference.replace('.fasta', '').replace('.fa', '') + '.dict'),
            bqsr_bam 
        )
    }

    variant_calling_results = VARIANT_CALLING(variant_calling_input)


    // --------------------------------------------------------
    // JOINT GENOTYPING: Combine trio GVCFs → joint VCF
    // --------------------------------------------------------
    // Inspect the results first
    // variant_calling_results.view { "🔍 VARIANT_CALLING tuple: ${it}" }

    // Convert to list safely
    joint_genotyping_input = variant_calling_results
        .toList()
        .map { gvcf_tuples ->
            // extract just the file paths (second element of each tuple)
            def gvcfs = gvcf_tuples.collect { tuple -> tuple[1] }
            // println "🧬 GVCF list collected: ${gvcfs}"

            if (gvcfs.size() != 3)
                throw new IllegalStateException("Expected exactly 3 GVCFs (trio), but found ${gvcfs.size()}")

            tuple(
                file(params.reference),
                file(params.reference + '.fai'),
                file(params.reference.replace('.fasta', '').replace('.fa', '') + '.dict'),
                file(gvcfs[0]),
                file(gvcfs[1]),
                file(gvcfs[2])
            )
        }

    joint_genotyping_results = JOINT_GENOTYPING(joint_genotyping_input)


    // --------------------------------------------------------
    // SELECT SNPS: Extract and filter SNPs from the joint VCF
    // --------------------------------------------------------
    select_snps_input = joint_genotyping_results.map { sid, joint_vcf ->
        tuple(
            sid,
            file(params.reference),
            file(params.reference + '.fai'),
            file(params.reference.replace('.fasta', '').replace('.fa', '') + '.dict'),
            file(joint_vcf)
        )
    }

    select_snps_results = SELECT_SNPS(select_snps_input)


    // --------------------------------------------------------
    // SELECT Indels: Extract and filter Indels from the joint VCF
    // --------------------------------------------------------
    select_indels_input = joint_genotyping_results.map { sid, joint_vcf ->
        tuple(
            sid,
            file(params.reference),
            file(params.reference + '.fai'),
            file(params.reference.replace('.fasta', '').replace('.fa', '') + '.dict'),
            file(joint_vcf)
        )
    }

    select_indels_results = SELECT_INDELS(select_indels_input)


    // --------------------------------------------------------
    // BCFTOOLS MERGE: Merge filtered SNPs and INDELs per sample
    // --------------------------------------------------------
    bcftools_merge_input = select_snps_results
        .join(select_indels_results)  // join by sample_id
        .map { sid, snps_tuple, indels_tuple ->
            tuple(
                sid,
                snps_tuple, 
                indels_tuple
            )
        }

    merged_vcf_results = BCFTOOLS_MERGE(bcftools_merge_input)


    // --------------------------------------------------------
    // VARIANT QC: compute stats and plots for merged VCF
    // --------------------------------------------------------
    qc_variant_input = merged_vcf_results.map { sid, merged_vcf ->
        tuple(sid, file(merged_vcf))
    }

    qc_variant_results = QC_VARIANT(qc_variant_input)

}
