version 1.0

import "qc_reads.wdl" as QCReadWDL
import "trim_fastp.wdl" as TrimWDL
import "alignment.wdl" as AlignWDL


workflow main {
    input {
        File sample11
        File sample12
        File sample21
        File sample22
        File sample31
        File sample32     
        String qc_dir
        File reference
        File reference_fai
        File reference_dict
        File fastq1
        File fastq2      
    }

    call TrimWDL.trim_fastp as myfastp1 {
            input: 
                sample_id = "SRR062634",
                fastq1 = sample11, 
                fastq2 = sample12,
                threads = 4
    }
    call TrimWDL.trim_fastp as myfastp2 {
            input: 
                sample_id = "SRR062635",
                fastq1 = sample21, 
                fastq2 = sample22,
                threads = 4
    }
    call TrimWDL.trim_fastp as myfastp3 {
            input: 
                sample_id = "SRR062637",
                fastq1 = sample31, 
                fastq2 = sample32,
                threads = 4
    }

    output {
        File sample1_trimmed1 = myfastp1.trimmed_fastq1
        File sample1_trimmed2 = myfastp1.trimmed_fastq1
        File sample1_report = myfastp1.report_html
        #FIle sample1_json = myfastp1.report_json

        File sample2_trimmed1 = myfastp2.trimmed_fastq1
        File sample2_trimmed2 = myfastp2.trimmed_fastq1
        File sample2_report = myfastp2.report_html
        #FIle sample2_json = myfastp2.report_json

        File sample3_trimmed1 = myfastp3.trimmed_fastq1
        File sample3_trimmed2 = myfastp3.trimmed_fastq1
        File sample3_report = myfastp3.report_html
        #FIle sample3_json = myfastp3.report_json

    }
}
