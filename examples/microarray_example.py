#!/usr/bin/env python3
# coding: utf-8

"""
Example script demonstrating the use of the genomekit MicroarrayGenerator module.
"""

import os
import sys
from genomekit.modules.microarray_generator import MicroarrayGenerator, AVAILABLE_FORMATS

def main():
    # Example usage of the MicroarrayGenerator
    
    # Replace with your actual input file (BAM or CRAM)
    input_file = "../data/examples/sample.bam"
    
    # Replace with your desired output directory
    output_dir = "./output"
    
    # Replace with your reference directory
    reference_dir = "../data/reference"
    
    # Check if the input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file not found: {input_file}")
        return 1
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Print available formats
    print("Available microarray formats:")
    for fmt in AVAILABLE_FORMATS:
        print(f"  - {fmt}")
    
    # Initialize the generator
    try:
        generator = MicroarrayGenerator(
            input_file=input_file,
            output_dir=output_dir,
            reference_dir=reference_dir,
            use_gatk=False  # Set to True to use GATK HaplotypeCaller instead of bcftools
        )
        
        # First generate the CombinedKit (contains all SNPs)
        print("Generating CombinedKit file...")
        success = generator.generate_combined_kit()
        if not success:
            print("Failed to generate CombinedKit.")
            return 1
        
        # Then generate a specific format (optional)
        print("\nGenerating 23andMe V5 format...")
        generator.generate_microarray_format("23andMe_V5")
        
        # Generate multiple formats at once (alternative approach)
        formats_to_generate = ["Ancestry_V2", "FTDNA_V3"]
        print(f"\nGenerating multiple formats: {', '.join(formats_to_generate)}")
        for format_name in formats_to_generate:
            print(f"Processing {format_name}...")
            generator.generate_microarray_format(format_name)
        
        print("\nMicroarray generation complete!")
        return 0
    
    except Exception as e:
        print(f"Error during microarray generation: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
