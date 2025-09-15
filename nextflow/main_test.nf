nextflow.enable.dsl=2

// -------------------------
// Load parameters from input.json
// -------------------------
params.samples   = params.samples
params.reference = params.reference
params.outdir    = params.outdir

// -------------------------
// Create channel of tuples
// -------------------------
samples_ch = Channel.fromList(params.samples.entrySet() as List)
    .map { entry -> 
        tuple(entry.key, file(entry.value[0]), file(entry.value[1]), file(params.reference)) 
    }

// -------------------------
// Alignment process
// -------------------------
process BWA_ALIGN {
    tag "$sample_id"
    cpus 4
    memory '6 GB'
    container 'wgs-align:latest'
    maxForks 2   // at most 2 samples in parallel

    input:
    tuple val(sample_id), path(fastq1), path(fastq2), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai")
    publishDir params.outdir, mode: 'copy'

    script:
    """
    bwa index $reference
    bwa mem -t ${task.cpus} $reference $fastq1 $fastq2 \
        | samtools view -Sb - \
        | samtools sort -o ${sample_id}.sorted.bam
    samtools index ${sample_id}.sorted.bam
    """
}

// -------------------------
// Workflow definition
// -------------------------
workflow {
    BWA_ALIGN(samples_ch)
}
