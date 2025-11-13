from pathlib import Path

def test_multiqc_report_exists():
    """Ensure MultiQC report is generated."""
    report = Path("Report/multiqc_report.html")
    assert report.exists(), "MultiQC report not found"
    assert report.stat().st_size > 1000, "MultiQC report too small, may be empty"
