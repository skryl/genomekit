#!/usr/bin/env python3
# coding: utf-8

"""
Example script demonstrating the use of the genomekit SNPChecker module.
"""

import os
import sys
from genomekit.modules.snp_checker import SNPChecker

def main():
    # Example usage of the SNPChecker
    
    # Replace with your actual input file (VCF or microarray file)
    input_file = "../data/examples/sample.vcf"
    
    # Replace with your reference directory
    reference_dir = "../data/reference"
    
    # Check if the input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file not found: {input_file}")
        return 1
    
    # Initialize the SNP checker
    try:
        checker = SNPChecker(
            input_file=input_file,
            reference_dir=reference_dir,
            debug=True  # Set to True for verbose output
        )
        
        # Run the analysis for a specific section
        # Options include: "metabolism", "inflammation", "cardiovascular", etc.
        checker.run_analysis(section="metabolism")
        
        print("\nAnalysis complete!")
        return 0
    
    except Exception as e:
        print(f"Error during SNP analysis: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
