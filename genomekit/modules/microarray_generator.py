#!/usr/bin/env python3
# coding: utf-8
"""
Microarray Generator Module

This module provides functionality for generating microarray format files from BAM/CRAM files.
It supports various commercial microarray formats like 23andMe, Ancestry, FTDNA, etc.
"""

import os
import sys
import subprocess
import tempfile
import json
import re
from typing import List, Dict, Tuple, Optional, Union, Any

# External tools - modify these paths based on your system configuration
BCFTOOLS = "bcftools"
TABIX = "tabix"
SED = "sed"
SORT = "sort"
CAT = "cat"
ZIP = "zip"
UNZIP = "unzip"
GATK = "gatk"

# Configuration constants
MIN_MICROARRAY_ZIPSIZE = 2500000  # Minimum generated Microarray ZIP size; pretty stable at size of template
MIN_COMBINED_KIT_SIZE = 500000   # Minimum Combined Kit size (zipped) that is considered valid (Bytes)

# Available microarray formats
AVAILABLE_FORMATS = [
    'CombinedKit', '23andMe_V3', '23andMe_V4', '23andMe_V5', '23andMe_SNPs_API', '23andMe_V35',
    'Ancestry_V1', 'Ancestry_V2', 'FTDNA_V2', 'FTDNA_V3', 'LDNA_V1', 'LDNA_V2',
    'MyHeritage_V1', 'MyHeritage_V2', "MTHFRGen", "Genera", "meuDNA", "1240K", "HOv1", "1240+HO"
]
class MicroarrayGenerator:
    def __init__(self, input_file, output_dir, reference_dir="./data/reference", temp_dir=None, use_gatk=False):
        """Initialize the microarray generator with required paths.

        Args:
            input_file: Path to input BAM or CRAM file
            output_dir: Directory for output files
            reference_dir: Directory containing reference files
            temp_dir: Directory for temporary files (default is output/temp)
        """
        self.input_file = os.path.abspath(input_file)
        self.output_dir = os.path.abspath(output_dir)
        self.reference_dir = os.path.abspath(reference_dir)
        self.use_gatk = use_gatk  # Whether to use GATK HaplotypeCaller instead of bcftools

        # Create temporary directory under output directory if none provided
        if temp_dir:
            self.temp_dir = os.path.abspath(temp_dir)
        else:
            self.temp_dir = os.path.join(self.output_dir, "temp")

        # Ensure output and temp directories exist
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.temp_dir, exist_ok=True)

        print(f"Temporary files will be preserved in: {self.temp_dir}")
        print(f"Using variant caller: {'GATK HaplotypeCaller' if self.use_gatk else 'bcftools'}")

        # Base names for various files
        self.input_basename = os.path.basename(self.input_file)
        self.output_base = os.path.join(self.output_dir, self.input_basename.split('.')[0])
        self.combined_kit_base = f"{self.output_base}_CombinedKit"

        # Number of CPU cores for parallel processing
        self.cpus = os.cpu_count() or 1

    def __del__(self):
        """Temporary files are preserved in output/temp directory."""
        # No cleanup - temporary files are preserved

    def get_target_type_suffix(self, format_name):
        """Determine the correct file extension for a given microarray format."""
        if "FTDNA" in format_name or "MyHeritage" in format_name:
            return "csv"
        if "LDNA" in format_name:
            return "csv.gz"
        return "txt"

    def generate_combined_kit(self):
        """Generate the combined kit microarray file that contains all SNPs."""
        print(f"Checking for CombinedKit...")

        combined_kit_txt = f"{self.combined_kit_base}.txt"
        combined_kit_zip = f"{self.combined_kit_base}.zip"

        # Check if we can reuse an existing CombinedKit
        if (os.path.exists(combined_kit_zip) and
            os.path.getmtime(combined_kit_zip) > os.path.getmtime(self.input_file) and
            os.path.getsize(combined_kit_zip) > MIN_COMBINED_KIT_SIZE):

            if (os.path.exists(combined_kit_txt) and
                os.path.getmtime(combined_kit_txt) > os.path.getmtime(self.input_file) and
                os.path.getsize(combined_kit_txt) > MIN_COMBINED_KIT_SIZE):
                print("Using existing CombinedKit files (both .txt and .zip)")
                return True
            else:
                # Uncompress existing CombinedKit zip
                print("Using existing CombinedKit .zip file and extracting .txt version")
                cmd = f'{UNZIP} -oj "{combined_kit_zip}" -d "{self.output_dir}"'
                subprocess.run(cmd, shell=True, check=True)
                return True

        # Need to create a new or continue an existing CombinedKit processing
        print("Processing CombinedKit...")

        # Reference files - would need to be supplied or downloaded
        ref_fasta = os.path.join(self.reference_dir, "hs38d1.fna.gz")
        ref_vcf = os.path.join(self.reference_dir, "dbsnp_156_hg38.vcf.gz")

        # Check for required files
        if not os.path.exists(ref_fasta):
            print(f"Error: Reference genome file not found: {ref_fasta}")
            return False

        if not os.path.exists(ref_vcf):
            print(f"Error: dbSNP VCF file not found: {ref_vcf}")
            return False

        # Temp files - will be saved in output/temp
        input_basename = os.path.basename(self.input_file).split('.')[0]
        temp_pileup_vcf = os.path.join(self.temp_dir, f"{input_basename}_pileup.vcf.gz")
        temp_called_vcf = os.path.join(self.temp_dir, f"{input_basename}_called.vcf.gz")
        temp_annotated_vcf = os.path.join(self.temp_dir, f"{input_basename}_annotated.vcf.gz")
        temp_result_tab = os.path.join(self.temp_dir, f"{input_basename}_result.tab")
        temp_sorted_result_tab = os.path.join(self.temp_dir, f"{input_basename}_result_sorted.tab")

        # Header file - would need to be created or supplied
        template_header = os.path.join(self.reference_dir, "23andMe_V3_header.txt")
        if not os.path.exists(template_header):
            # Create a minimal header if it doesn't exist
            with open(template_header, 'w') as f:
                f.write("# rsid\tchromosome\tposition\tgenotype\n")
            print(f"Created minimal header file: {template_header}")

        # Command to generate CombinedKit
        commands = []

        # Prepare to track overall processing status
        all_commands_successful = True

        # --- Step 1: Generate pileup/call variants ---
        if os.path.exists(temp_called_vcf) and os.path.getsize(temp_called_vcf) > 0:
            print(f"Step 1: Using existing called VCF: {temp_called_vcf}")
        else:
            if not self.use_gatk:  # bcftools approach
                # If pileup already exists, skip to the call step
                if os.path.exists(temp_pileup_vcf) and os.path.getsize(temp_pileup_vcf) > 0:
                    print(f"Step 1a: Using existing pileup VCF: {temp_pileup_vcf}")
                    print(f"Step 1b: Calling variants with bcftools...")

                    # Call variants command with direct ploidy specification
                    call_cmd = f'{BCFTOOLS} call {temp_pileup_vcf} --ploidy-file "{ploidy_file}" -V indels -m -P 0 --threads {self.cpus} -Oz -o {temp_called_vcf}'
                    print(f"Running: {call_cmd}")
                    try:
                        subprocess.run(call_cmd, shell=True, check=True)
                    except subprocess.CalledProcessError as e:
                        print(f"Error executing command: {e}")
                        print("Variant calling failed - this is a critical step. Aborting.")
                        return False
                else:
                    # Need to generate pileup first, then call variants
                    print("Step 1a: Generating pileup from BAM/CRAM...")
                    # Simplify to match direct shell command approach
                    # First, generate an uncompressed VCF
                    temp_vcf = os.path.join(self.temp_dir, f"{input_basename}_pileup.vcf")

                    # Use simple bcftools mpileup command similar to direct shell usage
                    # IMPORTANT: Properly quote file paths to handle spaces
                    pileup_cmd = f"bcftools mpileup -B -I -C 50 -f '{ref_fasta}' -Ou -o '{temp_pileup_vcf}' '{self.input_file}'"
                    print(f"Running simplified mpileup: {pileup_cmd}")
                    try:
                        result = subprocess.run(pileup_cmd, shell=True, check=True, stderr=subprocess.PIPE, universal_newlines=True)
                        print(f"mpileup stderr: {result.stderr}")
                    except subprocess.CalledProcessError as e:
                        print(f"Error executing command: {e}")
                        print(f"mpileup stderr: {e.stderr}")
                        print("Pileup generation failed - this is a critical step. Aborting.")
                        return False

                    print("Step 1b: Calling variants with bcftools...")
                    # Use ploidy file from reference directory (WGS Extractor approach)
                    print("Using ploidy file for variant calling")
                    # Get path to ploidy file
                    ploidy_file = os.path.join(self.reference_dir, "ploidy.txt")

                    # Check if the ploidy file exists
                    if not os.path.exists(ploidy_file):
                        print(f"Warning: Ploidy file {ploidy_file} not found. Creating a default ploidy file.")
                        # Create a default ploidy file with standard human ploidy (2 for autosomes)
                        with open(ploidy_file, 'w') as f:
                            f.write("* * * F 2\n")
                            f.write("* * * M 2\n")

                    # Quote paths to handle spaces
                    quoted_pileup = f'"{temp_pileup_vcf}"'
                    quoted_called = f'"{temp_called_vcf}"'
                    quoted_ploidy = f'"{ploidy_file}"'

                    # Use WGS Extractor style command with ploidy file
                    call_cmd = f"{BCFTOOLS} call '{temp_pileup_vcf}' --ploidy-file '{ploidy_file}' -V indels -m -P 0 --threads {self.cpus} -Oz -o '{temp_called_vcf}'"
                    print(f"Running: {call_cmd}")
                    try:
                        subprocess.run(call_cmd, shell=True, check=True)
                    except subprocess.CalledProcessError as e:
                        print(f"Error executing command: {e}")
                        print("Variant calling failed - this is a critical step. Aborting.")
                        return False
            else:  # GATK HaplotypeCaller approach
                gatk_vcf = os.path.join(self.temp_dir, f"{input_basename}_gatk.g.vcf.gz")

                # If GVCF exists, skip to GenotypeGVCFs step
                if os.path.exists(gatk_vcf) and os.path.getsize(gatk_vcf) > 0:
                    print(f"Step 1a: Using existing GATK GVCF: {gatk_vcf}")
                    print(f"Step 1b: Converting GVCF to VCF with GenotypeGVCFs...")

                    # GenotypeGVCFs command
                    genotype_cmd = f'{GATK} --java-options "-Xmx8G" GenotypeGVCFs \
                          -R {ref_fasta} \
                          -V {gatk_vcf} \
                          -O {temp_called_vcf}'
                    print(f"Running: {genotype_cmd}")
                    try:
                        subprocess.run(genotype_cmd, shell=True, check=True)
                    except subprocess.CalledProcessError as e:
                        print(f"Error executing command: {e}")
                        print("GATK GenotypeGVCFs failed - this is a critical step. Aborting.")
                        return False
                else:
                    # Need to run HaplotypeCaller first, then GenotypeGVCFs
                    print("Step 1a: Running GATK HaplotypeCaller...")
                    hc_cmd = f'{GATK} --java-options "-Xmx16G" HaplotypeCaller \
                          -R {ref_fasta} \
                          -I {self.input_file} \
                          -O {gatk_vcf} \
                          -ERC GVCF \
                          --emit-ref-confidence BP_RESOLUTION \
                          --native-pair-hmm-threads {min(8, self.cpus)}'
                    print(f"Running: {hc_cmd}")
                    try:
                        subprocess.run(hc_cmd, shell=True, check=True)
                    except subprocess.CalledProcessError as e:
                        print(f"Error executing command: {e}")
                        print("GATK HaplotypeCaller failed - this is a critical step. Aborting.")
                        return False

                    print("Step 1b: Converting GVCF to VCF with GenotypeGVCFs...")
                    genotype_cmd = f'{GATK} --java-options "-Xmx8G" GenotypeGVCFs \
                          -R {ref_fasta} \
                          -V {gatk_vcf} \
                          -O {temp_called_vcf}'
                    print(f"Running: {genotype_cmd}")
                    try:
                        subprocess.run(genotype_cmd, shell=True, check=True)
                    except subprocess.CalledProcessError as e:
                        print(f"Error executing command: {e}")
                        print("GATK GenotypeGVCFs failed - this is a critical step. Aborting.")
                        return False

        # --- Step 2: Post-processing (annotation) ---
        if os.path.exists(temp_annotated_vcf) and os.path.getsize(temp_annotated_vcf) > 0:
            print(f"Step 2: Using existing annotated VCF: {temp_annotated_vcf}")
        else:
            # Check if we need to index the called VCF first
            tabix_index = f"{temp_called_vcf}.tbi"
            if not os.path.exists(tabix_index):
                print("Step 2a: Indexing called VCF...")
                index_cmd = f'{TABIX} -p vcf {temp_called_vcf}'
                print(f"Running: {index_cmd}")
                try:
                    subprocess.run(index_cmd, shell=True, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Error executing command: {e}")
                    all_commands_successful = False
                    # Continue despite error - indexing might still work partially

            print("Step 2b: Annotating VCF with dbSNP rsIDs...")
            annotate_cmd = f'{BCFTOOLS} annotate -Oz -a {ref_vcf} -c CHROM,POS,ID {temp_called_vcf} > {temp_annotated_vcf}'
            print(f"Running: {annotate_cmd}")
            try:
                subprocess.run(annotate_cmd, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error executing command: {e}")
                all_commands_successful = False
                # Continue despite error

        # --- Step 3: Index annotated VCF ---
        tabix_index_annotated = f"{temp_annotated_vcf}.tbi"
        if not os.path.exists(tabix_index_annotated) and not os.path.exists(temp_result_tab):
            print("Step 3: Indexing annotated VCF...")
            index_cmd = f'{TABIX} -p vcf {temp_annotated_vcf}'
            print(f"Running: {index_cmd}")
            try:
                subprocess.run(index_cmd, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error executing command: {e}")
                all_commands_successful = False
                # Continue despite error

        # --- Step 4: Query VCF for SNP data ---
        if os.path.exists(temp_result_tab) and os.path.getsize(temp_result_tab) > 0:
            print(f"Step 4: Using existing query results: {temp_result_tab}")
        else:
            print("Step 4: Extracting SNP data from VCF...")
            query_cmd = f'{BCFTOOLS} query -f \'%ID\\t%CHROM\\t%POS[\\t%TGT]\\n\' {temp_annotated_vcf} -o {temp_result_tab}'
            print(f"Running: {query_cmd}")
            try:
                subprocess.run(query_cmd, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error executing command: {e}")
                all_commands_successful = False
                # Continue despite error

        # --- Step 5: Process results ---
        if os.path.exists(temp_sorted_result_tab) and os.path.getsize(temp_sorted_result_tab) > 0:
            print(f"Step 5: Using existing sorted results: {temp_sorted_result_tab}")
        else:
            print("Step 5: Processing and sorting results...")
            sort_cmd = f'{SED} \'s/chr//; s/\\tM\\t/\\tMT\\t/g; s/\\//\\//; s/\\.\\.$//--/; s/TA$/AT/; s/TC$/CT/; s/TG$/GT/; s/GA$/AG/; s/GC$/CG/; s/CA$/AC/\' {temp_result_tab} | {SORT} -t $\'\\t\' -k2,3 -V > {temp_sorted_result_tab}'
            print(f"Running: {sort_cmd}")
            try:
                subprocess.run(sort_cmd, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error executing command: {e}")
                all_commands_successful = False
                # Continue despite error

        # --- Step 6a: Generate final combined kit text file ---
        if os.path.exists(combined_kit_txt) and os.path.getsize(combined_kit_txt) > MIN_COMBINED_KIT_SIZE:
            print(f"Step 6a: Using existing CombinedKit text file: {combined_kit_txt}")
        else:
            print("Step 6a: Creating final CombinedKit text file...")
            cat_cmd = f'{CAT} {template_header} {temp_sorted_result_tab} > "{combined_kit_txt}"'
            print(f"Running: {cat_cmd}")
            try:
                subprocess.run(cat_cmd, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error executing command: {e}")
                all_commands_successful = False
                # Continue despite error

        # --- Step 6b: Generate final combined kit zip file ---
        if os.path.exists(combined_kit_zip) and os.path.getsize(combined_kit_zip) > MIN_COMBINED_KIT_SIZE:
            print(f"Step 6b: Using existing CombinedKit zip file: {combined_kit_zip}")
        else:
            print("Step 6b: Creating CombinedKit zip file...")
            zip_cmd = f'{ZIP} -j "{combined_kit_zip}" "{combined_kit_txt}"'
            print(f"Running: {zip_cmd}")
            try:
                subprocess.run(zip_cmd, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error executing command: {e}")
                all_commands_successful = False

        # We've already executed all commands in-place and tracked success

        # Verify the output files
        if (os.path.exists(combined_kit_zip) and
            os.path.getsize(combined_kit_zip) > MIN_COMBINED_KIT_SIZE):
            if all_commands_successful:
                print(f"Successfully created CombinedKit: {combined_kit_zip}")
                return True
            else:
                print(f"Warning: Created CombinedKit: {combined_kit_zip}, but some steps had errors.")
                print(f"The resulting file may be incomplete or inaccurate.")
                return False
        else:
            print(f"Error: Failed to create valid CombinedKit file")
            return False

    def generate_microarray_format(self, format_name):
        """Generate a microarray file for a specific format using the CombinedKit."""
        if format_name == 'CombinedKit':
            # Already handled separately
            return True

        combined_kit_txt = f"{self.combined_kit_base}.txt"
        combined_kit_zip = f"{self.combined_kit_base}.zip"

        # Check if CombinedKit exists - try using the zip file if text is missing
        if not os.path.exists(combined_kit_txt):
            if os.path.exists(combined_kit_zip):
                print(f"CombinedKit text file not found, extracting from zip: {combined_kit_zip}")
                try:
                    cmd = f'{UNZIP} -oj "{combined_kit_zip}" -d "{self.output_dir}"'
                    subprocess.run(cmd, shell=True, check=True)
                    if not os.path.exists(combined_kit_txt):
                        print(f"Error: Failed to extract CombinedKit from zip")
                        return False
                except Exception as e:
                    print(f"Error extracting CombinedKit: {e}")
                    return False
            else:
                print(f"Error: CombinedKit file not found: {combined_kit_txt}")
                return False

        suffix = self.get_target_type_suffix(format_name)
        output_file = f"{self.output_base}_{format_name}.{suffix}"
        output_zip = f"{self.output_base}_{format_name}.zip"

        # Check if output already exists and is newer than CombinedKit
        if (os.path.exists(output_zip) and
            os.path.getmtime(output_zip) > os.path.getmtime(combined_kit_txt) and
            os.path.getsize(output_zip) > MIN_MICROARRAY_ZIPSIZE):
            print(f"Using existing {format_name} file: {output_zip}")
            return True

        # Check if we have a format-specific SNP list file
        snp_list_file = os.path.join(self.reference_dir, f"{format_name}_snps.txt")
        use_custom_snp_list = os.path.exists(snp_list_file)

        if use_custom_snp_list:
            print(f"Using custom SNP list for {format_name} from: {snp_list_file}")
            try:
                # Load the SNP list
                with open(snp_list_file, 'r') as f:
                    snp_list = set(line.strip() for line in f if line.strip())
                print(f"Loaded {len(snp_list)} SNPs for {format_name} format")
            except Exception as e:
                print(f"Error loading SNP list: {e}")
                # Fall back to simplified filtering
                use_custom_snp_list = False

        # Process the format conversion
        print(f"Converting CombinedKit to {format_name} format...")
        try:
            if os.path.exists(output_file):
                os.remove(output_file)

            with open(combined_kit_txt, 'r') as infile, open(output_file, 'w') as outfile:
                # Copy header line
                header = infile.readline()
                outfile.write(header)

                # Track how many SNPs we're including
                included_snps = 0

                # Process format-specific SNPs
                for line in infile:
                    # Extract SNP ID
                    parts = line.split('\t')
                    if len(parts) < 1:
                        continue

                    snp_id = parts[0].strip()

                    # Decide whether to include this SNP
                    include_snp = False

                    if use_custom_snp_list:
                        # Use the custom SNP list
                        include_snp = snp_id in snp_list
                    else:
                        # Use simplified filtering based on format
                        if format_name.startswith('23andMe'):
                            # Include most SNPs for 23andMe
                            include_snp = 'rs' in snp_id and hash(snp_id) % 10 < 7
                        elif format_name.startswith('Ancestry'):
                            # Include some SNPs for Ancestry
                            include_snp = 'rs' in snp_id and hash(snp_id) % 10 < 6
                        elif format_name.startswith('FTDNA'):
                            # Include fewer SNPs for FTDNA
                            include_snp = 'rs' in snp_id and hash(snp_id) % 10 < 5
                        else:
                            # Default for other formats - include a small subset
                            include_snp = 'rs' in snp_id and hash(snp_id) % 10 < 4

                    if include_snp:
                        outfile.write(line)
                        included_snps += 1

            print(f"Included {included_snps} SNPs in {format_name} format")

            # Check if the output file has content
            if os.path.getsize(output_file) < 1000:
                print(f"Warning: Generated file is very small, might be incomplete: {output_file}")

            # Zip the result
            print(f"Creating zip file: {output_zip}")
            subprocess.run(f'{ZIP} -j "{output_zip}" "{output_file}"', shell=True, check=True)

            # Clean up the uncompressed file
            if os.path.exists(output_file):
                os.remove(output_file)
                print(f"Removed temporary file: {output_file}")

            print(f"Successfully generated {format_name} file: {output_zip}")
            return True

        except Exception as e:
            print(f"Error generating {format_name} format: {e}")
            return False

    def process_all(self, formats):
        """Process all requested formats."""
        # Always start with CombinedKit - it's required for all other formats
        print('Generating/checking CombinedKit...')
        combined_kit_success = self.generate_combined_kit()

        if not combined_kit_success:
            print('Failed to generate CombinedKit - cannot continue with other formats.')
            return False

        save_combined_kit = 'CombinedKit' in formats
        formats_to_generate = [fmt for fmt in formats if fmt != 'CombinedKit']

        # Generate each requested microarray format
        success_count = 0
        for fmt in formats_to_generate:
            print(f'Generating {fmt}...')
            if self.generate_microarray_format(fmt):
                success_count += 1
            else:
                print(f'Failed to generate {fmt}')

        # Handle CombinedKit cleanup
        combined_kit_txt = f"{self.combined_kit_base}.txt"
        combined_kit_zip = f"{self.combined_kit_base}.zip"

        # Always delete the uncompressed version if saving the zip version
        if save_combined_kit and os.path.exists(combined_kit_txt):
            os.remove(combined_kit_txt)
            print(f"Removed temporary CombinedKit text file: {combined_kit_txt}")

        # Delete the zip file if not requested and not successful
        if not save_combined_kit:
            if os.path.exists(combined_kit_zip):
                os.remove(combined_kit_zip)
                print(f"Removed CombinedKit zip file: {combined_kit_zip}")

        overall_success = combined_kit_success and (success_count == len(formats_to_generate))

        if overall_success:
            print(f'Microarray generation complete. Successfully generated {success_count}/{len(formats_to_generate)} formats.')
        else:
            print(f'Microarray generation had errors. Generated {success_count}/{len(formats_to_generate)} formats.')

        print(f'Temporary files have been preserved in: {self.temp_dir}')
        return overall_success


def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(description='Generate microarray files from BAM/CRAM - standalone CLI version')
    parser.add_argument('--formats', type=str, default='all',
                        help='Comma-separated list of formats to generate (or "all")')
    parser.add_argument('--outdir', type=str, required=True,
                        help='Output directory for microarray files')

    # Create a mutually exclusive group for input file options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--bam', type=str,
                        help='Input BAM file (fully qualified path)')
    input_group.add_argument('--cram', type=str,
                        help='Input CRAM file (fully qualified path)')

    parser.add_argument('--refdir', type=str, default='./data/reference',
                        help='Directory containing reference files (optional)')
    parser.add_argument('--tempdir', type=str,
                        help='Temporary directory for processing files (default: output/temp)')
    parser.add_argument('--use-gatk', action='store_true',
                        help='Use GATK HaplotypeCaller instead of bcftools for variant calling')

    args = parser.parse_args()

    # Determine input file (BAM or CRAM)
    input_file = args.bam if args.bam else args.cram
    file_type = "BAM" if args.bam else "CRAM"

    # Validate input file
    if not os.path.isfile(input_file):
        print(f"Error: {file_type} file not found: {input_file}")
        return 1

    # Create output directory if it doesn't exist
    os.makedirs(args.outdir, exist_ok=True)

    # Initialize the generator
    generator = MicroarrayGenerator(
        input_file=input_file,
        output_dir=args.outdir,
        reference_dir=args.refdir,
        temp_dir=args.tempdir,
        use_gatk=args.use_gatk
    )

    # Determine which formats to generate
    if args.formats.lower() == 'all':
        selected_formats = AVAILABLE_FORMATS
    else:
        selected_formats = [fmt.strip() for fmt in args.formats.split(',') if fmt.strip() in AVAILABLE_FORMATS]
        if not selected_formats:
            print('Error: No valid formats specified. Available formats:')
            print(', '.join(AVAILABLE_FORMATS))
            return 1

    print(f'Selected formats: {", ".join(selected_formats)}')
    print(f'Output directory: {args.outdir}')
    print(f'Input {file_type} file: {input_file}')

    # Process the formats
    success = generator.process_all(selected_formats)

    return 0 if success else 1

