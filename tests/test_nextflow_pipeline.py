import subprocess
from pathlib import Path
import pytest

PIPELINE_DIR = Path(__file__).resolve().parents[1]

def test_nextflow_dry_run():
    """Validate Nextflow syntax and configuration."""
    result = subprocess.run(
        ["nextflow", "lint", "nextflow/main.nf"],
        cwd=PIPELINE_DIR,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"Nextflow lint failed:\n{result.stderr}"


@pytest.mark.timeout(1200)
def test_nextflow():
    # Run small sample workflow and confirm output files exist.
    outdir = PIPELINE_DIR / "Result"
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "nextflow", "run", "nextflow/main.nf",
        "-params-file", "tests/test_input.json",
        "-profile", "docker",
        "-c", "nextflow.config",
        "--outdir", str(outdir),
    ]
    result = subprocess.run(cmd, cwd=PIPELINE_DIR, capture_output=True, text=True)
    assert result.returncode == 0, f"Pipeline run failed: {result.stderr}"

    # Ensure MultiQC report is generated.
    qc_report = outdir / "multiqc_report.html"
    assert qc_report.exists(), "MultiQC report not found"
    assert qc_report.stat().st_size > 1000, "MultiQC report too small, may be empty"
    print(f"✅ MultiQC report is found: {qc_report}")


    # Ensure trimed fastq files are generated.
    trimmed_files = list(outdir.glob("**/*.trimmed.fastq"))
    assert len(trimmed_files) > 0, "No trimmed FASTQ files found"
    for f in trimmed_files:
        assert f.stat().st_size > 0, f"Trimmed FASTQ file is empty: {f}"
    print(f"✅ Trimmed FASTQ files are found: {len(trimmed_files)} files")


    # Ensure BWA-aligned BAM files are generated.
    bwa_files = list(outdir.glob("**/*.sorted.bam"))
    assert len(bwa_files) > 0, "No BWA-aligned BAM files found"
    for f in bwa_files:
        assert f.stat().st_size > 0, f"BWA-aligned BAM file is empty: {f}"
    print(f"✅ BWA-aligned BAM files are found: {len(bwa_files)} files")


    # Ensure GATK-QC files are generated.
    qc_files = list(outdir.glob("**/*.alignment_metrics.txt"))
    assert len(qc_files) > 0, "No GATK-QC files found"
    for f in qc_files:
        assert f.stat().st_size > 0, f"GATK-QC file is empty: {f}"
    print(f"✅ GATK-QC files are found: {len(qc_files)} files")


    # Ensure GATK-dedup BAM files are generated.
    dedup_files = list(outdir.glob("**/*.dedup.bam"))
    assert len(dedup_files) > 0, "No GATK-dedup BAM files found"
    for f in dedup_files:
        assert f.stat().st_size > 0, f"GATK-dedup BAM file is empty: {f}"
    print(f"✅ GATK-dedup BAM files are found: {len(dedup_files)} files")


    # Ensure GATK-bqsr BAM files are generated.
    bqsr_files = list(outdir.glob("**/*.bqsr.bam"))
    assert len(bqsr_files) > 0, "No GATK-bqsr BAM files found"
    for f in bqsr_files:
        assert f.stat().st_size > 0, f"GATK-bqsr BAM file is empty: {f}"
    print(f"✅ GATK-bqsr BAM files are found: {len(bqsr_files)} files")


    # Ensure raw variant calling files are generated.
    rvc_files = list(outdir.glob("**/*.raw_variants.g.vcf.gz"))
    assert len(rvc_files) > 0, "No raw variant calling files found"
    for f in rvc_files:
        assert f.stat().st_size > 0, f"Raw variant calling file is empty: {f}"
    print(f"✅ Raw variant calling files are found: {len(rvc_files)} files")


    # Check combined raw variant calling file is generated.
    raw_vcf = outdir / "raw_variants.vcf.gz"
    assert raw_vcf.exists(), f"Combined raw variant calling VCF not found: {raw_vcf}"
    assert raw_vcf.stat().st_size > 0, f"Combined raw variant calling VCF file is empty: {raw_vcf}"
    print(f"✅ Combined raw variant calling VCF is found: {raw_vcf}")


    # Check merged variant calling file is generated.
    merged_vcf = outdir / "merged_variants.vcf"
    assert merged_vcf.exists(), f"Merged variant calling VCF not found: {merged_vcf}"
    assert merged_vcf.stat().st_size > 0, f"Merged variant calling VCF file is empty: {merged_vcf}"
    print(f"✅ Merged variant calling VCF is found: {merged_vcf}")


    # Check variant calling QC report is generated.
    vcf_report = outdir / "qc_plots.zip"
    assert vcf_report.exists(), f"Variant calling QC report not found: {vcf_report}"
    assert vcf_report.stat().st_size > 0, f"Variant calling QC report file is empty: {vcf_report}"
    print(f"✅ Variant calling QC report file is found: {vcf_report}")

