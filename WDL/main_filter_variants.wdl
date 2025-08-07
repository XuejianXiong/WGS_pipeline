version 1.0

import "joint_genotyping.wdl" as JointWDL
import "select_snps.wdl" as SNPsWDL
import "select_indels.wdl" as INDELsWDL
import "bcftools_merge.wdl" as BCftoolsWDL
import "qc_variant.wdl" as QCvarWDL


workflow main {
    input {
        File reference
        File reference_fai
        File reference_dict
        File known_variants_snps
        File known_variants_indels
        File gvcf_1
        File gvcf_2
        File gvcf_3
    }

    call JointWDL.joint_genotyping as my_vjoint {
        input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            gvcf_1 = gvcf_1,
            gvcf_2 = gvcf_2,
            gvcf_3 = gvcf_3
    }

    call SNPsWDL.select_snps as my_snps {
        input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            raw_vcf = my_vjoint.raw_vcf
    }
    
    call INDELsWDL.select_indels as my_indels {
        input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            raw_vcf = my_vjoint.raw_vcf
    }    

    call BCftoolsWDL.bcftools_merge as my_vmerge {
        input:
            filtered_snps = my_snps.snps_vcf,
            filtered_indels = my_indels.indels_vcf
    }  

    call QCvarWDL.variant_qc as my_qcvar {
        input:
            vcf_file = my_vmerge.merged_vcf
    } 

    output {
        File file_joint_vcf = my_vjoint.raw_vcf
        File file_snps_vcf = my_snps.snps_vcf
        File file_indels_vcf = my_indels.indels_vcf
        File file_merged_vcf = my_vmerge.merged_vcf   
        File file_variant_stats = my_qcvar.variant_stats
        File file_variant_plots = my_qcvar.zipped_plots            
    }
}
