#!/usr/bin/env python3
# coding: utf-8

"""
Example script demonstrating the use of the genomekit OpenCravatAnalyzer module.
"""

import os
import sys
from genomekit.modules.opencravat_analyzer import OpenCravatAnalyzer

def main():
    # Example usage of the OpenCravatAnalyzer
    
    # Replace with your actual input VCF file
    vcf_file = "../data/examples/sample.vcf"
    
    # Replace with your desired output directory
    output_dir = "./clinical_analysis"
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Check if the input file exists
    if not os.path.exists(vcf_file):
        print(f"Error: VCF file not found: {vcf_file}")
        return 1
    
    try:
        # Initialize the OpenCravat analyzer
        analyzer = OpenCravatAnalyzer(
            output_dir=output_dir,
            debug=True  # Set to True for verbose output
        )
        
        # Step 1: Run OpenCravat analysis on the VCF file
        print("Running OpenCravat analysis on VCF file...")
        annotators = ["clinvar", "gnomad", "cosmic", "omim", "dbsnp"] 
        success = analyzer.run_analysis(
            vcf_file=vcf_file,
            annotators=annotators
        )
        
        if not success:
            print("Failed to run OpenCravat analysis.")
            return 1
        
        print("\nOpenCravat analysis complete!")
        
        # Step 2: Find clinically significant variants in the results
        print("\nFinding clinically significant variants...")
        db_path = os.path.join(output_dir, "analysis.sqlite")
        output_csv = os.path.join(output_dir, "clinvar_variants.csv")
        
        success = analyzer.find_clinvar_variants(
            db_path=db_path,
            output_csv=output_csv,
            limit=100  # Maximum number of variants to return
        )
        
        if not success:
            print("Failed to find clinically significant variants.")
            return 1
        
        print(f"\nClinVar variants analysis complete!")
        print(f"Results saved to: {output_csv}")
        
        # Optionally, show the top findings from the CSV
        if os.path.exists(output_csv):
            print("\nTop findings (first 5 lines):")
            with open(output_csv, 'r') as f:
                for i, line in enumerate(f):
                    if i < 5:
                        print(line.strip())
                    else:
                        break
        
        return 0
        
    except Exception as e:
        print(f"Error during OpenCravat analysis: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
