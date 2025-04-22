#!/usr/bin/env bash

# This script analyzes VCF files using Open-Cravat with all clinical annotators
# installed on the system (CADD, ClinVar, gnomAD, COSMIC, CancerHotspots, OMIM, dbSNP)

# Default output directory
CLINICAL_OUTPUT="clinical_analysis"

# Parse command line arguments
usage() {
  echo "Usage: $0 [-o output_dir] <vcf_file>"
  echo "Options:"
  echo "  -o <dir>  Specify output directory (default: clinical_analysis)"
  echo "Example: $0 -o my_results allSites.rs.vcf.gz"
  exit 1
}

# Process options
while getopts "o:h" opt; do
  case $opt in
    o) CLINICAL_OUTPUT="$OPTARG" ;;
    h) usage ;;
    \?) usage ;;
  esac
done

# Shift to the positional arguments
shift $((OPTIND-1))

# Check if input file was provided as positional argument
if [[ $# -lt 1 ]]; then
  echo "ERROR: No VCF file specified"
  usage
fi

# Get input VCF file from first positional argument
VCF_FILE="$1"

# Check if input file exists
if [[ ! -f "$VCF_FILE" ]]; then
  echo "ERROR: Input VCF file not found: $VCF_FILE" >&2
  exit 1
fi
mkdir -p "$CLINICAL_OUTPUT"

# 1. Run a single analysis with all clinical annotators at once
echo "Step 1/2: Running comprehensive clinical annotation..."
echo "This may take several hours depending on your system..."

echo "Running oc analysis with all clinical annotators..."

# Try without specifying genome reference since hg38 is already in base modules
oc run "$VCF_FILE" \
  -a clinvar gnomad cosmic omim dbsnp \
  -l hg38 \
  -d "$CLINICAL_OUTPUT" \
  -n "analysis" \

# 2. Generate comprehensive clinical report
if [[ -f "$CLINICAL_OUTPUT/analysis.sqlite" ]]; then
  echo "Step 2/2: Generating comprehensive clinical report..."
  oc report "$CLINICAL_OUTPUT/analysis.sqlite" -f html
  echo "Comprehensive clinical analysis complete. Results are in $CLINICAL_OUTPUT directory."

  echo ""
  echo "==============================================="
  echo "All analyses complete!"
  echo "Analysis: $CLINICAL_OUTPUT/analysis.sqlite"
  echo ""
  echo "To analyze clinical variants, run:"
  echo "python analyze_oc_report.py $CLINICAL_OUTPUT/analysis.sqlite"
  echo "==============================================="
else
  echo "ERROR: Comprehensive clinical analysis failed, no output database created."
  exit 1
fi