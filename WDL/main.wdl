version 1.0

import "qc_reads.wdl" as QCReadWDL
import "trim_fastp.wdl" as TrimWDL
import "alignment.wdl" as AlignWDL
import "dedup.wdl" as DedupWDL
import "qc_gatk.wdl" as QCgatkWDL
import "bqsr.wdl" as BqsrWDL
import "variant_calling.wdl" as VarCallWDL


workflow main {
    input {
        File fastq1
        File fastq2
        String qc_dir
        File reference
        File reference_fai
        File reference_dict
        File trim_fastq1
        File trim_fastq2 
        File known_variants_snps
        File known_variants_indels
    }

    #call QCReadWDL.qc_fastqc as my_fastqc {
    #        input: 
    #            fastq1 = fastq1, 
    #            fastq2 = fastq2,
    #            qc_dir = qc_dir
    #}

    #call TrimWDL.trim_fastp as my_fastp {
    #        input: 
    #            sample_id = "SRR062634",
    #            fastq1 = fastq1, 
    #            fastq2 = fastq2,
    #            threads = 4
    #}



    call AlignWDL.alignment as my_align {
        input:
            reference = reference,
            fastq1 = trim_fastq1,
            fastq2 = trim_fastq2
    }
    
    call DedupWDL.dedup as my_dedup {
        input:
            aligned_sam = my_align.sam
            #aligned_sam = aligned_sam
    }

    call QCgatkWDL.qc_gatk as my_qc {
        input:
            reference = reference,
            dedup_bam = my_dedup.bam
            #dedup_bam = dedup_bam
    }      

    call BqsrWDL.bqsr as my_bqsr {
        input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            known_variants_snps = known_variants_snps,
            known_variants_indels = known_variants_indels,
            dedup_bam = my_dedup.bam
            #dedup_bam = dedup_bam
    }

    call VarCallWDL.variant_calling as my_varcall {
        input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            bqsr_bam = my_bqsr.bqsr_bam
            #bqsr_bam = bqsr_bam
    }    

    output {
        #File fastqc1 = my_fastqc.fastqc_output1
        #File fastqc2 = my_fastqc.fastqc_output2
        #File trim_output = my_fastp.output_tar

        File align_output = my_align.sam
        File dedup_report = my_dedup.dedup_report
        File qc_hist_pdf = my_qc.hist_pdf        
        File file_bqsr_bam = my_bqsr.bqsr_bam
        File file_bqsr_report = my_bqsr.bqsr_report
        File file_varcall_vcf = my_varcall.gvcf
    }
}
