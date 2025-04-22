#!/usr/bin/env python3
# coding: utf8
# Copyright (C) 2018-2020 City Farmer
# Copyright (C) 2020-2022 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

import os
import sys
import argparse
import subprocess
import tempfile
import json
import re
import csv
import sqlite3
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Union, Any

# Import modules from the genomekit package
from genomekit.modules.microarray_generator import (
    MicroarrayGenerator, AVAILABLE_FORMATS,
    MIN_MICROARRAY_ZIPSIZE, MIN_COMBINED_KIT_SIZE,
    BCFTOOLS, TABIX, SED, SORT, CAT, ZIP, UNZIP, GATK
)
from genomekit.modules.opencravat_analyzer import OpenCravatAnalyzer
from genomekit.modules.snp_checker import SNPChecker

"""
Genome Analysis CLI Tool

A unified command-line interface for various genome analysis tasks.
"""

def main():
    # Create main parser
    parser = argparse.ArgumentParser(description="Genome Analysis CLI Tool")
    subparsers = parser.add_subparsers(dest="command", help="Command to run")

    # --- Microarray Format Command ---
    microarray_parser = subparsers.add_parser("microarray", help="Generate microarray files from BAM/CRAM")

    # Group for mutually exclusive input options
    input_group = microarray_parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--bam', type=str,
                            help='Path to BAM file input')
    input_group.add_argument('--cram', type=str,
                            help='Path to CRAM file input')

    microarray_parser.add_argument('--formats', type=str, nargs='+', required=True,
                        choices=AVAILABLE_FORMATS,
                        help='Microarray format(s) to generate')
    microarray_parser.add_argument('--outdir', type=str, required=True,
                        help='Output directory for generated files')
    microarray_parser.add_argument('--refdir', type=str, default="./data/reference",
                        help='Directory containing reference files (optional)')
    microarray_parser.add_argument('--tempdir', type=str,
                        help='Temporary directory for processing files (default: output/temp)')
    microarray_parser.add_argument('--use-gatk', action='store_true',
                        help='Use GATK HaplotypeCaller instead of bcftools for variant calling')

    # --- SNP Check Command ---
    snp_parser = subparsers.add_parser("snps", help="Check SNPs in a VCF or microarray file")
    snp_parser.add_argument('input', type=str,
                           help='Path to annotated VCF file or microarray text file')
    snp_parser.add_argument('--section', type=str, default="all",
                           help='Section to analyze (metabolism, inflammation, cardiovascular, etc.)')
    snp_parser.add_argument('--refdir', type=str, default="./data/reference",
                           help='Directory containing reference files')
    snp_parser.add_argument('--debug', action='store_true',
                           help='Print debug information')
                           
    # --- OpenCravat Command ---
    oc_parser = subparsers.add_parser("oc", help="OpenCravat analysis commands")
    oc_subparsers = oc_parser.add_subparsers(dest="action", help="OpenCravat action to run")
    
    # OpenCravat Run Command
    oc_run_parser = oc_subparsers.add_parser("run", help="Run OpenCravat analysis on a VCF file")
    oc_run_parser.add_argument('--vcf', type=str, required=True,
                              help='Path to VCF file for annotation')
    oc_run_parser.add_argument('--outdir', type=str, default="./clinical_analysis",
                              help='Output directory for OpenCravat results')
    oc_run_parser.add_argument('--annotators', type=str, nargs='+',
                              default=["clinvar", "gnomad", "cosmic", "omim", "dbsnp"],
                              help='Annotators to use (default: clinvar gnomad cosmic omim dbsnp)')
    oc_run_parser.add_argument('--debug', action='store_true',
                              help='Print debug information')
    
    # OpenCravat Find ClinVar Variants Command
    oc_find_parser = oc_subparsers.add_parser("find-clinvar", help="Find high-quality ClinVar variants in OpenCravat results")
    oc_find_parser.add_argument('--db', type=str,
                               help='Path to OpenCravat SQLite database (default: <outdir>/analysis.sqlite)')
    oc_find_parser.add_argument('--csv', type=str,
                               help='Path to output CSV file (default: <outdir>/clinvar_variants.csv)')
    oc_find_parser.add_argument('--outdir', type=str, default="./clinical_analysis",
                               help='Output directory for OpenCravat results')
    oc_find_parser.add_argument('--limit', type=int, default=100,
                               help='Maximum number of variants to return (default: 100)')
    oc_find_parser.add_argument('--debug', action='store_true',
                               help='Print debug information')

    args = parser.parse_args()

    if args.command == "microarray":
        # Handle microarray command
        # Determine input file
        if args.bam:
            input_file = args.bam
            if not os.path.isfile(input_file):
                print(f"Error: BAM file not found: {input_file}")
                return 1
        elif args.cram:
            input_file = args.cram
            if not os.path.isfile(input_file):
                print(f"Error: CRAM file not found: {input_file}")
                return 1

        # Initialize the generator
        generator = MicroarrayGenerator(
            input_file=input_file,
            output_dir=args.outdir,
            reference_dir=args.refdir,
            temp_dir=args.tempdir,
            use_gatk=args.use_gatk
        )

        # Determine which formats to generate
        formats_to_generate = []

        # If CombinedKit is specified, handle that first
        if "CombinedKit" in args.formats:
            formats_to_generate.append("CombinedKit")
            # Also remove CombinedKit from the list if present to avoid duplication
            args.formats = [f for f in args.formats if f != "CombinedKit"]

        # Add all remaining formats
        formats_to_generate.extend(args.formats)

        # Generate each format
        for format_name in formats_to_generate:
            if format_name == "CombinedKit":
                generator.generate_combined_kit()
            else:
                generator.generate_microarray_format(format_name)

        print("Microarray generation completed successfully.")

    elif args.command == "snps":
        # Handle SNP check command
        try:
            checker = SNPChecker(
                input_file=args.input,
                reference_dir=args.refdir,
                debug=args.debug
            )
            checker.run_analysis(section=args.section)
        except Exception as e:
            print(f"Error during SNP analysis: {e}")
            return 1
    
    elif args.command == "oc":
        # Handle OpenCravat subcommand
        if args.action == "run":
            # Run OpenCravat analysis
            if not args.vcf:
                print("Error: VCF file is required for OpenCravat analysis")
                return 1
                
            try:
                analyzer = OpenCravatAnalyzer(
                    output_dir=args.outdir,
                    debug=args.debug
                )
                success = analyzer.run_analysis(
                    vcf_file=args.vcf,
                    annotators=args.annotators
                )
                return 0 if success else 1
            except Exception as e:
                print(f"Error during OpenCravat analysis: {e}")
                return 1
                
        elif args.action == "find-clinvar":
            # Find ClinVar variants in OpenCravat results
            try:
                analyzer = OpenCravatAnalyzer(
                    output_dir=args.outdir,
                    debug=args.debug
                )
                success = analyzer.find_clinvar_variants(
                    db_path=args.db,
                    output_csv=args.csv,
                    limit=args.limit
                )
                return 0 if success else 1
            except Exception as e:
                print(f"Error analyzing ClinVar variants: {e}")
                return 1
        else:
            # Invalid OpenCravat action
            oc_parser.print_help()
            return 1

    else:
        # No command specified
        parser.print_help()
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
