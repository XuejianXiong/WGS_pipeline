#!/opt/homebrew/bin/bash
# =============================================================================
# build_all.sh
# -----------------------------------------------------------------------------
# Build all Docker images required for the WGS_pipeline.
# - Ensures base image builds first (dependency for other images)
# - Builds other images in parallel
# - Supports --force to rebuild all images regardless of local cache
# -----------------------------------------------------------------------------
# NOTE: The shebang points to Homebrew-installed Bash on macOS.
#       This ensures compatibility with modern Bash features (e.g., associative arrays)
#       because macOS default Bash is outdated (v3.x).
# -----------------------------------------------------------------------------
# Professional MLOps practices:
# - Parallel builds speed up CI/CD and local testing
# - Strict error handling
# - Clear logging for maintainability
# =============================================================================

set -euo pipefail  # Exit on error, unset variable, or failed pipe

# -----------------------------------------------------------------------------
# Ensure the script runs relative to its own directory
# -----------------------------------------------------------------------------
cd "$(dirname "${BASH_SOURCE[0]}")"
DOCKER_DIR="../Docker"

# -----------------------------------------------------------------------------
# Check for --force flag
# -----------------------------------------------------------------------------
FORCE_REBUILD=false
if [[ "${1:-}" == "--force" ]]; then
  FORCE_REBUILD=true
  echo "⚡ Force rebuild enabled: all images will be rebuilt."
fi

# -----------------------------------------------------------------------------
# Declare Docker images and corresponding Dockerfiles
# -----------------------------------------------------------------------------
declare -A IMAGES=(
  [base]="Dockerfile.base"
  [qc-reads]="Dockerfile.qc_reads"
  [align]="Dockerfile.alignment"
  [fastp]="Dockerfile.fastp"
  [gatk]="Dockerfile.gatk"
  [bcftools]="Dockerfile.bcftools"
)

# -----------------------------------------------------------------------------
# Function to check if a Docker image exists locally
# -----------------------------------------------------------------------------
image_exists() {
  local image_name=$1
  docker image inspect "$image_name" > /dev/null 2>&1
}

# -----------------------------------------------------------------------------
# 1️⃣ Build base image first
# -----------------------------------------------------------------------------
BASE_TAG="wgs-base:latest"
if [ "$FORCE_REBUILD" = false ] && image_exists "$BASE_TAG"; then
  echo "⚡ Skipping base image (already exists locally)"
else
  echo "🚀 Building base image..."
  docker build -t "$BASE_TAG" -f "$DOCKER_DIR/${IMAGES[base]}" "$DOCKER_DIR"
fi

# -----------------------------------------------------------------------------
# 2️⃣ Build remaining images in parallel
# -----------------------------------------------------------------------------
pids=()  # Track background process IDs
for name in qc-reads align fastp gatk bcftools; do
  image_tag="wgs-$name:latest"
  if [ "$FORCE_REBUILD" = false ] && image_exists "$image_tag"; then
    echo "⚡ Skipping $name image (already exists locally)"
  else
    echo "🚀 Building $name image..."
    docker build -t "$image_tag" -f "$DOCKER_DIR/${IMAGES[$name]}" "$DOCKER_DIR" &
    pids+=($!)  # Store PID of background build
  fi
done

# -----------------------------------------------------------------------------
# Wait for all parallel builds to finish
# -----------------------------------------------------------------------------
if [ ${#pids[@]} -gt 0 ]; then
  echo "⏳ Waiting for all Docker builds to complete..."
  wait "${pids[@]}"
fi

# -----------------------------------------------------------------------------
# Completion message
# -----------------------------------------------------------------------------
echo "🎉 All required Docker images are now available!"
