import subprocess
import pytest

@pytest.mark.parametrize("image", [
    "bwa:latest",
    "gatk:latest",
    "samtools:latest",
])
def test_docker_image_available(image):
    """Ensure required Docker images are built or pullable."""
    result = subprocess.run(["docker", "image", "inspect", image],
                            capture_output=True, text=True)
    assert result.returncode == 0, f"Docker image not found: {image}"
