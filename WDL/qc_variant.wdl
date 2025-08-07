version 1.0

task variant_qc {
    input {
        File vcf_file
    }

    command {
        set -euo pipefail

        # Compress and index VCF
        bgzip -f ${vcf_file}
        tabix -f -p vcf ${vcf_file}.gz

        # Compute stats
        bcftools stats ${vcf_file}.gz > variant_stats.txt

        # Generate plots
        mkdir -p qc_plots
        plot-vcfstats -p qc_plots variant_stats.txt

        # Zip the plots directory
        zip -r qc_plots.zip qc_plots
    }

    output {
        File variant_stats = "variant_stats.txt"
        File zipped_plots = "qc_plots.zip"
    }

    runtime {
        cpu: 6
        memory: "8G"
        docker: "wgs-bcftools"
    }
}
