version 1.0

import "WDL/qc_reads.wdl" as QCReadWDL


workflow main {
    input {
        File sample11
        File sample12
        File sample21
        File sample22
        File sample31
        File sample32     
        String qc_dir      
    }

    call QCReadWDL.fastqc_1sample as myfastqc1 {
            input: 
                fastq1 = sample11, 
                fastq2 = sample12,
                qc_dir = qc_dir
    }
    call QCReadWDL.fastqc_1sample as myfastqc2 {
            input: 
                fastq1 = sample21, 
                fastq2 = sample22,
                qc_dir = qc_dir
    }
    call QCReadWDL.fastqc_1sample as myfastqc3 {
            input: 
                fastq1 = sample31, 
                fastq2 = sample32,
                qc_dir = qc_dir
    }

    #call QCReadWDL.multiqc_summary as mymultiqc {
    #    input: 
    #        qc_dir = myfastqc3.out_qc_dir
    #}

    output {
        File sample1_fastqc1 = myfastqc1.fastqc_output1
        File sample1_fastqc2 = myfastqc1.fastqc_output2
        File sample2_fastqc1 = myfastqc2.fastqc_output1
        File sample2_fastqc2 = myfastqc2.fastqc_output2
        File sample3_fastqc1 = myfastqc3.fastqc_output1
        File sample3_fastqc2 = myfastqc3.fastqc_output2
        #File report = mymultiqc.multiqc_report
    }
}
