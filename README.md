# GenomeKit

A comprehensive Python package for genome analysis, including SNP checking, microarray generation, and clinical variant annotation. GenomeKit provides a unified, modular toolkit for processing genomic data using various reference databases and third-party tools.

## Features

- **Microarray Generation:** Convert BAM/CRAM files to various microarray formats (23andMe, Ancestry, FTDNA) with smart checkpointing to resume interrupted operations and avoid redundant processing.

- **SNP Checking:** Query clinically relevant SNPs with parallel processing for VCF files and display color-coded genotype interpretations (protective, risk, carrier status) with detailed explanations.

- **Clinical Variant Analysis:** Integrate OpenCravat for comprehensive variant annotation with data from multiple sources (ClinVar, gnomAD, COSMIC, OMIM) and generate reports of clinically significant variants.


## Quick Start

### 1. Clone the Repository
```bash
git clone git@github.com:skryl/genome.git
cd genome
```

### 2. Install Dependencies and Reference Data
```bash
bash install_dependencies.sh
```
- Use `bash install_dependencies.sh --test` for a dry run (shows actions, does not modify system).

### 3. Run the CLI Tool

You can use the CLI directly without installation using the included `genome` script:

```bash
# Make the script executable (first time only)
chmod +x genome

# View available commands
./genome --help
```

Or install the package and use it globally:

```bash
# Install the package
pip install .

# Run commands
genomekit --help
```
## Command Examples

### Microarray Generation

Generate microarray files from BAM/CRAM files:

```bash
# Generate a specific microarray format
./genome microarray --bam path/to/file.bam --formats CombinedKit--outdir ./output

# Generate multiple formats
./genome microarray --cram path/to/file.cram --formats CombinedKit 23andMe_V5 Ancestry_V2 --outdir ./output

# Use GATK instead of bcftools
./genome microarray --bam path/to/file.bam --formats 23andMe_V5 --outdir ./output --use-gatk
```

### SNP Checking

Check SNPs in a VCF or microarray file:

```bash
# Check all SNPs in a 23andMe microarray file
./genome snps path/to/snp_report.txt

# Check all SNPs in a VCF file
./genome snps path/to/annotated_vcf.vcf.gz

# Check SNPs for a specific section
./genome snps path/to/annotated_vcf.vcf.gz --section metabolism

# Use a different reference directory (must contain snp_catalog.json)
./genome snps path/to/annotated_vcf.vcf.gz --refdir /path/to/reference --debug
```

### OpenCravat Analysis

Run OpenCravat analysis on VCF files and extract clinical variants:

```bash
# Run basic OpenCravat analysis
./genome oc run --vcf path/to/file.vcf --outdir ./clinical_analysis

# Run with specific annotators
./genome oc run --vcf path/to/file.vcf --annotators clinvar gnomad cosmic --outdir ./clinical_analysis

# Find clinically significant variants
./genome oc find-clinvar --outdir ./clinical_analysis --limit 50
```

## Directory Structure
```
├── data/
│   ├── dbsnp/                 # dbSNP VCF files (e.g., dbsnp_156_hg38.vcf.gz)
│   ├── reference/             # Reference genomes (e.g., hs38d1.fna.gz)
├── scripts/                   # Analysis and utility scripts
├── requirements.txt           # Python dependencies
├── install_dependencies.sh    # Automated setup script
└── README.md                  # This file
```

## Requirements
- **System:** bash, wget, bcftools, htslib, samtools, python3, (optionally: Homebrew or apt-get)
- **Python:** See `requirements.txt`
- **Tools:** OpenCravat (oc), GATK (manual install for some workflows)

## Notes
- Reference and dbSNP files are only downloaded if not already present.
- All VCF and FASTA files are indexed for rapid access.
- OpenCravat modules (`clinvar`, `gnomad`, `cosmic`, `omim`, `dbsnp`) are installed automatically.

## Data Processing Pipeline

The analysis pipeline uses the unified CLI tool (via the top-level `genome` script) that supports microarray generation, SNP checking, and OpenCravat analysis.

### 1. Generate Microarray Files from BAM/CRAM

```bash
# Generate microarray files from CRAM (default bcftools approach)
./genome microarray \
  --cram data/nebula/sample.cram \
  --formats CombinedKit \
  --outdir data/output/microarray \
  --refdir data/reference

# Generate multiple formats at once
./genome microarray \
  --bam data/nebula/sample.bam \
  --formats CombinedKit 23andMe_V5 Ancestry_V2 FTDNA_V3 \
  --outdir data/output/microarray

# Alternative: Use GATK HaplotypeCaller for variant calling
./genome microarray \
  --cram data/nebula/sample.cram \
  --formats CombinedKit 23andMe_V5 \
  --outdir data/output/microarray \
  --use-gatk
```

The tool will:
- Process the BAM/CRAM file to call variants
- Generate microarray files in requested formats (e.g., 23andMe_V5, Ancestry_V2)
- Reuse existing intermediate files when possible for faster processing
- Store all intermediate files in the specified output directory's temp subdirectory

### 2. Analyze SNPs from VCF or Microarray Files

```bash
# Analyze all SNP categories from microarray
./genome snps data/output/microarray/CombinedKit.txt

# Analyze all SNP categories from VCF
./genome snps data/output/microarray/sample_annotated.vcf.gz

# Check only specific SNP categories
./genome snps data/output/microarray/CombinedKit.txt --section metabolism
./genome snps data/output/microarray/CombinedKit.txt --section cardiovascular

# Available sections: metabolism, inflammation, cardiovascular, neurotransmitters,
# antioxidant, vitamin, drug, athletic, longevity, sleep, cognitive, nutrition
```

The SNP analysis will:
- Query your VCF or microarray file for clinically relevant SNPs
- Display color-coded results (green/yellow/red) based on genotype interpretation
- Provide detailed explanations for each variant
- Support strand complementarity checking for accurate genotyping

### 3. Run OpenCravat Analysis

```bash
# Run OpenCravat analysis on a VCF file
./genome oc run --vcf data/output/microarray/sample_annotated_vcf.vcf.gz --outdir data/clinical_analysis

# Run with specific annotators
./genome oc run --vcf data/output/microarray/sample_annotated.vcf \
  --annotators clinvar gnomad cosmic omim dbsnp \
  --outdir data/clinical_analysis

# Find clinically significant variants from OpenCravat results
./genome oc find-clinvar --outdir data/clinical_analysis --limit 50

# List all installed OpenCravat annotators
./genome oc list-annotators

# Generate report with all variants and OpenCravat annotations
./genome oc generate-report --outdir data/clinical_analysis --format excel
```

The OpenCravat analysis will:
- Annotate variants in your VCF file with data from various sources (ClinVar, gnomAD, COSMIC, OMIM, dbSNP)
- Generate comprehensive reports in various formats (TSV, Excel, PDF)
- Identify clinically significant variants by pathogenicity and review status
- Allow filtering based on multiple annotation sources

## License
MIT License. See [LICENSE](LICENSE) for details.

## Author
[skryl](https://github.com/skryl)
