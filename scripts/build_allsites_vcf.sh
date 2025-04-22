#!/usr/bin/env bash
set -euo pipefail

# Usage: ./build_allsites_vcf.sh sample.vcf.gz dbsnp.vcf.gz GRCh38.p14.fa

# Input parameters
VCF="$1"
DBSNP="$2"
REF="$3"
PREFIX="allSites"

# Check if files exist
if [[ ! -f $VCF ]]; then
  echo "ERROR: VCF not found: $VCF" >&2; exit 1
fi
if [[ ! -f $DBSNP ]]; then
  echo "ERROR: dbSNP VCF not found: $DBSNP" >&2; exit 1
fi
if [[ ! -f $REF ]]; then
  echo "ERROR: Reference FASTA not found: $REF" >&2; exit 1
fi

echo "1) Annotating VCF with dbSNP rsIDs"
bcftools annotate \
  -a "$DBSNP" -c ID \
  -Oz -o "${PREFIX}.rs.vcf.gz" \
  "$VCF"

echo "2) Indexing annotated VCF"
tabix -p vcf "${PREFIX}.rs.vcf.gz"

echo "Done.  You now have:"
echo "  • ${PREFIX}.rs.vcf.gz      (with rsIDs in the ID column)"
echo
echo "Next, run your SNP report:"
echo "  ./check_snps.sh ${PREFIX}.rs.vcf.gz  > my_snp_report.tsv"
 