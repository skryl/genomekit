#!/usr/bin/env bash
set -euo pipefail

# Input parameters
YOUR_VCF=$1
DBSNP_VCF=$2
REF_GENOME=$3
OUTPUT_DIR=$4

mkdir -p "$OUTPUT_DIR"

echo "Step 1: Creating sites-only dbSNP VCF..."
bcftools view -G "$DBSNP_VCF" -Oz -o "$OUTPUT_DIR/dbsnp_sites_only.vcf.gz"
bcftools index "$OUTPUT_DIR/dbsnp_sites_only.vcf.gz"

echo "Step 2: Filling in your variants where they exist..."
bcftools annotate -a "$YOUR_VCF" -c CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE \
  "$OUTPUT_DIR/dbsnp_sites_only.vcf.gz" \
  -Oz -o "$OUTPUT_DIR/dbsnp_plus_your_variants.vcf.gz"

echo "Step 3: Filling missing positions with reference genotypes..."
bcftools +fill-tags "$OUTPUT_DIR/dbsnp_plus_your_variants.vcf.gz" -Oz -o "$OUTPUT_DIR/complete_genome.vcf.gz" -- -t all

echo "Step 4: Indexing the final VCF..."
bcftools index "$OUTPUT_DIR/complete_genome.vcf.gz"

echo "Complete genome VCF created at $OUTPUT_DIR/complete_genome.vcf.gz"