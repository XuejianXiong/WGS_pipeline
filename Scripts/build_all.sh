#!/usr/bin/env bash
set -e

# Build Docker images for each pipeline module
echo "Building base image..."
docker build -t wgs-base:latest -f ../Docker/Dockerfile.base .

echo "Building QC image..."
docker build -t wgs-qc-reads:latest -f ../Docker/Dockerfile.qc_reads .

echo "Building alignment image..."
docker build -t wgs-align:latest -f ../Docker/Dockerfile.alignment .

echo "Building fastp image..."
docker build -t wgs-fastp:latest -f ../Docker/Dockerfile.fastp .

echo "Building gatk image..."
docker build -t wgs-gatk:latest -f ../Docker/Dockerfile.gatk .

echo "Building bcftools image..."
docker build -t wgs-bcftools:latest -f ../Docker/Dockerfile.bcftools .

echo "✅ All Docker images built successfully!"

