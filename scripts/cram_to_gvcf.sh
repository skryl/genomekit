#!/usr/bin/env bash
set -euo pipefail

CRAM="/Users/skryl/Documents/backups/genome/Nebula Genomics/NG1Y0BHNRM.mm2.sortdup.bqsr.cram"
CRAI="${CRAM}.crai"
REF="GRCh38.p14.fa"
OUT="/Users/skryl/Documents/dev/projects/genome/NG1Y0BHNRM.mm2.sortdup.bqsr.g.vcf.gz"

# 1. Ensure the CRAM is indexed (creates .crai)
if [[ ! -f $CRAI ]]; then
  echo "Indexing CRAM…"
  samtools index "$CRAM"
else
  echo "CRAM index found: $CRAI"
fi

# 2. Run GATK HaplotypeCaller in GVCF mode
echo "Calling variants into GVCF…"
gatk --java-options "-Xmx16G" HaplotypeCaller \
  -R "$REF" \
  -I "$CRAM" \
  -O "$OUT" \
  -ERC GVCF \
  --EMIT-ALL-SITES \
  --native-pair-hmm-threads 4

# 3. Compress & index the GVCF
echo "Indexing GVCF…"
tabix -p vcf "$OUT"

echo
echo "✅ gVCF built and indexed as $OUT and ${OUT}.tbi"