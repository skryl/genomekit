#!/usr/bin/env python3
# coding: utf-8
"""
OpenCravat Analyzer Module

This module provides functionality for running OpenCravat analysis on VCF files 
and processing the results to identify clinically significant variants.

OpenCravat is a powerful annotation tool for genomic variants that integrates 
data from various sources like ClinVar, gnomAD, COSMIC, etc.
"""

import os
import sys
import subprocess
import csv
import sqlite3
from typing import List, Dict, Tuple, Optional, Union, Any

class OpenCravatAnalyzer:
    def __init__(self, output_dir: str = "./clinical_analysis", debug: bool = False):
        """Initialize the OpenCravat analyzer.

        Args:
            output_dir: Directory for OpenCravat output files
            debug: Whether to print debug information
        """
        self.output_dir = os.path.abspath(output_dir)
        self.debug = debug

        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)
    
    def _debug_echo(self, message: str) -> None:
        """Print debug information if debug mode is enabled."""
        if self.debug:
            print(f"DEBUG: {message}", file=sys.stderr)
    
    def run_analysis(self, vcf_file: str, annotators: List[str] = None) -> bool:
        """Run the OpenCravat analysis on a VCF file.

        Args:
            vcf_file: Path to the VCF file to analyze
            annotators: List of annotators to use (default: clinvar, gnomad, cosmic, omim, dbsnp)

        Returns:
            True if analysis was successful, False otherwise
        """
        if not os.path.isfile(vcf_file):
            print(f"Error: Input VCF file not found: {vcf_file}")
            return False
        
        # Set default annotators if none provided
        if not annotators:
            annotators = ["clinvar", "gnomad", "cosmic", "omim", "dbsnp"]
        
        # Create the OpenCravat command
        output_db = os.path.join(self.output_dir, "analysis.sqlite")
        
        print("Running OpenCravat analysis with clinical annotators...")
        print("This may take several hours depending on your system...")
        
        # Prepare the command arguments
        args = [
            "oc", "run", vcf_file,
            "-l", "hg38",
            "-d", self.output_dir,
            "-n", "analysis"
        ]
        
        # Add each annotator separately
        for annotator in annotators:
            args.extend(["-a", annotator])
        
        # Run the command
        try:
            self._debug_echo(f"Running command: {' '.join(args)}")
            subprocess.run(args, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error executing OpenCravat: {e}")
            return False
        
        # Generate HTML report if analysis succeeded
        if os.path.isfile(output_db):
            print("Generating comprehensive clinical report...")
            try:
                subprocess.run(["oc", "report", output_db, "-f", "html"], check=True)
                print(f"OpenCravat analysis complete. Results are in {self.output_dir} directory.")
                return True
            except subprocess.CalledProcessError as e:
                print(f"Error generating report: {e}")
                return False
        else:
            print("Error: OpenCravat analysis failed, no output database created.")
            return False
    
    def find_clinvar_variants(self, db_path: str = None, output_csv: str = None, limit: int = 100) -> bool:
        """Find high-quality ClinVar variants in the OpenCravat results.

        Args:
            db_path: Path to the OpenCravat SQLite database (default: analysis.sqlite in output_dir)
            output_csv: Path to output CSV file (default: clinvar_variants.csv in output_dir)
            limit: Maximum number of variants to return (default: 100)

        Returns:
            True if analysis was successful, False otherwise
        """
        # Use default database path if not specified
        if not db_path:
            db_path = os.path.join(self.output_dir, "analysis.sqlite")
        
        # Use default output CSV path if not specified
        if not output_csv:
            output_csv = os.path.join(self.output_dir, "clinvar_variants.csv")
        
        if not os.path.isfile(db_path):
            print(f"Error: OpenCravat database not found: {db_path}")
            return False
        
        try:
            # Connect to the database
            conn = sqlite3.connect(db_path)
            cursor = conn.cursor()
            
            print(f"Connected to {db_path}")
            
            # Get total variants count
            cursor.execute("SELECT COUNT(*) FROM variant")
            total_variants = cursor.fetchone()[0]
            print(f"Total variants in database: {total_variants:,}")
            
            # Check what columns are available
            cursor.execute("PRAGMA table_info(variant)")
            columns = [col[1] for col in cursor.fetchall()]
            self._debug_echo(f"Available columns: {', '.join(columns)}")
            
            # Check if clinvar columns exist
            has_clinvar = any(col.startswith('clinvar__') for col in columns)
            
            if not has_clinvar:
                print("Warning: ClinVar annotations not found in database. Looking for other informative variants...")
                # Fall back to dbSNP if available
                if 'dbsnp__rsid' in columns:
                    print("Using dbSNP data to find potentially important variants...")
                    base_cols = [col for col in columns if col.startswith('base__') or col.startswith('dbsnp__')]
                    common_cols = ', '.join(base_cols)
                    query = f"SELECT {common_cols} FROM variant WHERE dbsnp__rsid IS NOT NULL LIMIT ?"
                    cursor.execute(query, (limit,))
                    results = cursor.fetchall()
                    
                    # Save results to CSV
                    with open(output_csv, 'w', newline='') as f:
                        writer = csv.writer(f)
                        writer.writerow([col.replace('__', '_') for col in base_cols])
                        for row in results:
                            writer.writerow(row)
                    
                    print(f"Found {len(results)} variants with dbSNP IDs.")
                    print(f"Saved results to {output_csv}")
                    return True
                else:
                    print("No suitable annotation data found for generating a clinical report.")
                    print("Please re-run OpenCravat with the clinvar annotator.")
                    return False
            
            # Run queries to find high-quality variants
            all_results = []
            
            # Build the base column list that we know exists
            base_cols = ["base__chrom", "base__pos", "base__ref_base", "base__alt_base"]
            if "base__hugo" in columns: base_cols.append("base__hugo")
            if "base__so" in columns: base_cols.append("base__so")
            
            # Add optional columns that might exist
            optional_cols = [
                "clinvar__sig", "clinvar__disease_names", "clinvar__rev_stat", 
                "clinvar__id", "dbsnp__rsid", "gnomad__af"
            ]
            
            found_cols = base_cols.copy()
            for col in optional_cols:
                if col in columns:
                    found_cols.append(col)
            
            # Build the SELECT clause
            select_clause = ", ".join(found_cols)
            
            # Try to find highest quality variants first - expert panel reviewed
            print("Searching for highest quality variants (expert panel review)...")
            try:
                # Only add these conditions if the clinvar columns exist
                if "clinvar__sig" in columns and "clinvar__rev_stat" in columns:
                    where_clause = """
                        clinvar__sig IS NOT NULL
                        AND (clinvar__sig LIKE '%pathogenic%' OR clinvar__sig LIKE '%Pathogenic%')
                        AND (clinvar__rev_stat LIKE '%expert panel%' OR clinvar__rev_stat LIKE '%practice guideline%')
                    """
                    
                    query = f"""
                    SELECT {select_clause} 
                    FROM variant
                    WHERE {where_clause}
                    LIMIT ?
                    """
                    
                    cursor.execute(query, (limit // 3,))
                    highest_results = cursor.fetchall()
                    for row in highest_results:
                        # Add quality level indicator
                        all_results.append(row + ("Highest",))
                    
                    print(f"Found {len(highest_results)} variants with expert panel review")
                    
                    # Find variants with multiple submitters
                    print("Searching for high quality variants (multiple submitters, no conflicts)...")
                    where_clause = """
                        clinvar__sig IS NOT NULL
                        AND (clinvar__sig LIKE '%pathogenic%' OR clinvar__sig LIKE '%Pathogenic%')
                        AND clinvar__rev_stat LIKE '%multiple submitters%'
                    """
                    
                    query = f"""
                    SELECT {select_clause} 
                    FROM variant
                    WHERE {where_clause}
                    LIMIT ?
                    """
                    
                    cursor.execute(query, (limit // 3,))
                    high_results = cursor.fetchall()
                    for row in high_results:
                        all_results.append(row + ("High",))
                    
                    print(f"Found {len(high_results)} variants with multiple submitters")
                    
                    # Find variants with single submitter
                    print("Searching for moderate quality variants (single submitter)...")
                    where_clause = """
                        clinvar__sig IS NOT NULL
                        AND (clinvar__sig LIKE '%pathogenic%' OR clinvar__sig LIKE '%Pathogenic%')
                        AND clinvar__rev_stat LIKE '%single submitter%'
                    """
                    
                    query = f"""
                    SELECT {select_clause} 
                    FROM variant
                    WHERE {where_clause}
                    LIMIT ?
                    """
                    
                    cursor.execute(query, (limit // 3,))
                    moderate_results = cursor.fetchall()
                    for row in moderate_results:
                        all_results.append(row + ("Moderate",))
                    
                    print(f"Found {len(moderate_results)} variants with single submitter")
                else:
                    # Just get some variants with rs IDs if clinvar not available
                    query = f"""
                    SELECT {select_clause} 
                    FROM variant
                    WHERE dbsnp__rsid IS NOT NULL
                    LIMIT ?
                    """
                    
                    cursor.execute(query, (limit,))
                    basic_results = cursor.fetchall()
                    for row in basic_results:
                        all_results.append(row + ("Unknown",))
                    
                    print(f"Found {len(basic_results)} variants with dbSNP IDs")
                
                # Write results to CSV
                if all_results:
                    print(f"Writing {len(all_results)} variants to {output_csv}")
                    with open(output_csv, 'w', newline='') as f:
                        writer = csv.writer(f)
                        
                        # Write header - convert column names from database format to readable format
                        header = [col.replace('__', '_') for col in found_cols]
                        header.append("Quality_Level")  # Add quality level column
                        writer.writerow(header)
                        
                        # Write data
                        for row in all_results:
                            writer.writerow(row)
                    
                    print(f"Analysis complete. Found {len(all_results)} variants.")
                    print(f"Results saved to {output_csv}")
                    return True
                else:
                    print("No significant variants found.")
                    return False
                    
            except sqlite3.OperationalError as e:
                print(f"Error analyzing OpenCravat results: {e}")
                # If the specific error is about a missing column, suggest running with correct annotators
                if "no such column" in str(e):
                    print("The database is missing expected columns. Please re-run OpenCravat with the right annotators:")
                    print("./genome oc run --vcf your_vcf_file.vcf --annotators clinvar gnomad cosmic omim dbsnp --outdir output_dir")
                return False
        except Exception as e:
            print(f"Error analyzing OpenCravat results: {e}")
            return False


