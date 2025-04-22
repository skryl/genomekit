#!/usr/bin/env bash
# install_dependencies.sh — Install all project dependencies and download reference data
# Usage: bash install_dependencies.sh [--test]

set -euo pipefail

# ----------- Parse flags -----------
TEST_MODE=false
for arg in "$@"; do
  if [[ "$arg" == "--test" ]]; then
    TEST_MODE=true
  fi
done
if $TEST_MODE; then
  echo "[TEST MODE] No actions will be performed. Only log statements will be shown."
fi

# ----------- Ensure data directories exist -----------
mkdir -p data/reference

# ----------- System Dependencies -----------
echo "[1/5] Installing system packages..."
if ! $TEST_MODE; then
    if [[ "$OSTYPE" == "darwin"* ]]; then
        # macOS
        command -v brew >/dev/null || { echo "Homebrew required. Install from https://brew.sh/"; exit 1; }
        brew install coreutils gnu-sed wget curl parallel bcftools htslib samtools || true # Add more if needed
    else
        # Linux
        sudo apt-get update && sudo apt-get install -y wget curl parallel gawk bcftools htslib samtools
    fi
fi

# ----------- GATK (Genome Analysis Toolkit) -----------
echo "[EXTRA] Checking for GATK..."
if ! command -v gatk &>/dev/null; then
    echo "GATK not found. Please download and install GATK manually from: https://github.com/broadinstitute/gatk/releases"
    echo "(Optional, but required for some scripts: cram_to_gvcf.sh)"
else
    echo "GATK found: $(command -v gatk)"
fi

# ----------- OpenCravat (oc) -----------
echo "[EXTRA] Checking for OpenCravat (oc)..."
if ! command -v oc &>/dev/null; then
    echo "Installing OpenCravat (oc) via pip..."
    if ! $TEST_MODE; then python3 -m pip install --user open-cravat; fi
else
    echo "OpenCravat found: $(command -v oc)"
fi

# ----------- OpenCravat modules -----------
echo "[EXTRA] Installing required OpenCravat modules (clinvar, gnomad, cosmic, omim, dbsnp)..."
for mod in clinvar gnomad cosmic omim dbsnp; do
    if oc module ls | grep -E "^$mod\s" | grep -q installed; then
        echo "oc module '$mod' is already installed."
    else
        echo "Installing oc module: $mod ..."
        if ! $TEST_MODE; then oc module install $mod || { echo "Failed to install oc module: $mod"; exit 1; }; fi
    fi
done

# ----------- Extra tools -----------
echo "[EXTRA] Checking for bcftools, tabix, samtools, bgzip, zcat, gunzip, awk, grep, sed, cut, sort, head, tail, parallel..."
for tool in bcftools tabix samtools bgzip zcat gunzip awk grep sed cut sort head tail parallel; do
    command -v $tool >/dev/null || echo "WARNING: $tool not found in PATH!"
done


# ----------- Python Dependencies -----------
echo "[2/5] Installing Python packages..."
if [[ -f requirements.txt ]]; then
    if ! $TEST_MODE; then python3 -m pip install --user -r requirements.txt; fi
fi

# Optional: install tkmacosx for macOS GUI scripts
if [[ "$OSTYPE" == "darwin"* ]]; then
    if ! $TEST_MODE; then python3 -m pip install --user tkmacosx || true; fi
fi

# ----------- Data Directory Setup -----------
echo "[3/5] Ensuring data/ directory structure..."
# Directories are already created at top of script

# ----------- Download Reference Genome (GRCh38) -----------
REF_FASTA="data/reference/hs38d1.fna.gz"
echo "Checking for reference genome: $REF_FASTA"
if [[ -f $REF_FASTA ]]; then
    echo "Reference genome already exists: $REF_FASTA"
else
    echo "[4/5] Downloading human reference genome (GRCh38) as hs38d1.fna.gz..."
    if ! $TEST_MODE; then wget -O "$REF_FASTA" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz"; fi
fi

# ----------- Index reference genome with samtools faidx -----------
REF_FASTA_FAI="${REF_FASTA}.fai"
echo "Checking for reference genome index: $REF_FASTA_FAI"
if [[ -f $REF_FASTA_FAI ]]; then
    echo "Reference genome index already exists: $REF_FASTA_FAI"
else
    echo "Indexing reference genome with samtools faidx..."
    if ! $TEST_MODE; then samtools faidx "$REF_FASTA"; fi
fi

# ----------- Download dbSNP VCF (build 156, hg38) -----------
dbSNP_VCF="data/reference/dbsnp_156_hg38.vcf.gz"
echo "Checking for dbSNP VCF: $dbSNP_VCF"
if [[ -f $dbSNP_VCF ]]; then
    echo "dbSNP VCF already exists: $dbSNP_VCF"
else
    echo "[5/5] Downloading dbSNP (build 156, hg38)..."
    if ! $TEST_MODE; then wget -O "$dbSNP_VCF" "https://ftp.ncbi.nih.gov/snp/archive/b156/VCF/GCF_000001405.40.gz"; fi
fi

# ----------- Index dbSNP VCF with tabix -----------
dbSNP_TBI="${dbSNP_VCF}.tbi"
echo "Checking for dbSNP VCF index: $dbSNP_TBI"
if [[ -f $dbSNP_TBI ]]; then
    echo "dbSNP VCF index already exists: $dbSNP_TBI"
else
    echo "Indexing dbSNP VCF with tabix..."
    if ! $TEST_MODE; then tabix -p vcf "$dbSNP_VCF"; fi
fi

echo "\nAll dependencies and reference data installed/downloaded."
