#!/usr/bin/env python3
# coding: utf-8
"""
SNP Checker Module

This module provides functionality for checking SNPs in VCF or microarray files 
and reporting their genotypes with interpretations. It supports color-coded output 
to easily identify protective, carrier, and risk genotypes.

The module can work with both VCF files and microarray text files from services 
like 23andMe, Ancestry, etc. It includes capabilities for strand complementarity 
checking to ensure accurate genotype reporting.
"""

import os
import sys
import subprocess
import json
import re
import multiprocessing
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Tuple, Optional, Union, Any

class SNPChecker:
    def __init__(self, input_file: str, reference_dir: str = "./data/reference", debug: bool = False):
        """Initialize the SNP checker with required paths.

        Args:
            input_file: Path to the input file (VCF or microarray text file)
            reference_dir: Directory containing reference files
            debug: Whether to print debug information
        """
        self.input_file = os.path.abspath(input_file)
        self.reference_dir = os.path.abspath(reference_dir)
        self.debug = debug

        # Reference genome path
        self.ref_genome = os.path.join(self.reference_dir, "GRCh38.p14.fa")

        # ANSI color codes for output
        self.colors = {
            "GREEN": "\033[42m\033[30m",  # black text on green background
            "YELLOW": "\033[43m\033[30m", # black text on yellow background
            "RED": "\033[41m\033[37m",    # white text on red background
            "GRAY": "\033[37m",          # gray text (for variant/unknown)
            "NC": "\033[0m"               # reset to normal
        }

        # Determine input file type (VCF or microarray)
        self.is_vcf = self._is_vcf_file()
        self._debug_echo(f"Input file detected as: {'VCF' if self.is_vcf else 'Microarray text file'}")
        
        # If it's a microarray file, preload the data
        self.microarray_data = {}
        if not self.is_vcf:
            self._load_microarray_data()

        # Load SNP catalog from JSON file
        self.catalog = self._load_snp_catalog()

        # Ensure the input file exists
        if not os.path.isfile(self.input_file):
            raise FileNotFoundError(f"Input file not found: {self.input_file}")

    def _is_vcf_file(self) -> bool:
        """Determine if the input file is a VCF file based on content or extension."""
        # Check file extension first
        if self.input_file.endswith(".vcf") or self.input_file.endswith(".vcf.gz"):
            return True
        
        # If not obvious from extension, check content
        try:
            # Try to read the first few lines
            if self.input_file.endswith(".gz"):
                import gzip
                with gzip.open(self.input_file, 'rt') as f:
                    first_line = f.readline().strip()
            else:
                with open(self.input_file, 'r') as f:
                    first_line = f.readline().strip()
            
            # VCF files typically start with ## or have a specific header format
            return first_line.startswith("##") or "#CHROM" in first_line
        except Exception as e:
            self._debug_echo(f"Error checking file type: {e}")
            # Default to assuming it's a VCF if we can't determine
            return True
    
    def _load_microarray_data(self):
        """Load data from a microarray text file into memory."""
        try:
            with open(self.input_file, 'r') as f:
                # Skip header line if present
                first_line = f.readline()
                if not first_line.startswith("rs"):
                    self._debug_echo("Skipping header line in microarray file")
                else:
                    # Process first line if it's not a header
                    self._process_microarray_line(first_line)
                
                # Process the rest of the file
                for line in f:
                    self._process_microarray_line(line)
                    
            self._debug_echo(f"Loaded {len(self.microarray_data)} SNPs from microarray file")
        except Exception as e:
            print(f"Error loading microarray file: {e}")
            sys.exit(1)
    
    def _get_complement(self, base: str) -> str:
        """Return the DNA complement of a base."""
        complements = {
            'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
            'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
            '-': '-', 'N': 'N', 'n': 'n'
        }
        return complements.get(base, base)
        
    def _process_microarray_line(self, line: str):
        """Process a single line from a microarray file."""
        line = line.strip()
        if not line or line.startswith("#"):
            return
        
        # Common microarray formats have tab or space-separated fields
        # Usually: rsid, chromosome, position, genotype
        parts = re.split(r'[\t ,]', line)
        if len(parts) < 3:
            return
        
        # Extract rsid and genotype
        rsid = None
        raw_genotype = None
        
        # Look for rsID pattern in the parts
        for part in parts:
            if part.startswith("rs") and re.match(r'rs\d+', part):
                rsid = part
                break
        
        # Try to find genotype in different positions based on format
        # Format 1: rsid chromosome position genotype (23andMe)
        if len(parts) >= 4 and rsid == parts[0]:
            possible_gt = parts[3]
            if re.match(r'^[ACGTN-]{1,2}$', possible_gt, re.IGNORECASE):
                raw_genotype = possible_gt.upper()
        
        # Format 2: Look for explicit genotype field with slashes or pipes
        if not raw_genotype:
            for possible_gt in parts:
                if re.match(r'^[ACGTN-][/|][ACGTN-]$', possible_gt, re.IGNORECASE):
                    # Extract just the bases
                    bases = re.sub(r'[/|]', '', possible_gt.upper())
                    raw_genotype = bases
                    break
        
        # Format 3: Look for any 1-2 character field that could be a genotype
        if not raw_genotype:
            for possible_gt in parts[1:]:
                if re.match(r'^[ACGTN-]{1,2}$', possible_gt, re.IGNORECASE):
                    raw_genotype = possible_gt.upper()
                    break
        
        # Process and store the raw genotype if found
        if rsid and raw_genotype:
            # Validate raw genotype has valid bases
            if re.match(r'^[ACGTN-]{1,2}$', raw_genotype, re.IGNORECASE):
                # Store the raw genotype in our dictionary (exactly like the shell script)
                # This ensures we're starting with the same raw data
                self.microarray_data[rsid] = raw_genotype
                self._debug_echo(f"Loaded SNP {rsid}: {raw_genotype}")
            else:
                self._debug_echo(f"Invalid genotype format for {rsid}: {raw_genotype}")
    
    def _load_snp_catalog(self) -> Dict:
        """Load the SNP catalog from the JSON file."""
        catalog_path = os.path.join(self.reference_dir, "snp_catalog.json")
        try:
            with open(catalog_path, 'r') as f:
                return json.load(f)
        except Exception as e:
            print(f"Error loading SNP catalog: {e}")
            sys.exit(1)

    def _debug_echo(self, message: str) -> None:
        """Print debug information if debug mode is enabled."""
        if self.debug:
            print(f"DEBUG: {message}", file=sys.stderr)

    def _fetch_genotype(self, chrpos: str, rsid: str) -> Tuple[str, str]:
        """Fetch genotype for a given SNP by rsID or position.

        Args:
            chrpos: Chromosome position (format: chr:pos)
            rsid: rsID of the SNP

        Returns:
            Tuple of (genotype, actual_genotype)
        """
        self._debug_echo(f"Checking {rsid} at position {chrpos}")
        
        # If this is a microarray file, look for the rsID in our preloaded data
        if not self.is_vcf:
            if rsid in self.microarray_data:
                raw_genotype = self.microarray_data[rsid]
                self._debug_echo(f"Found {rsid} in microarray data: {raw_genotype}")
                
                # Validate the raw genotype format (should be 1-2 characters like AA, AG, etc.)
                if re.match(r'^[ACGTN-]{1,2}$', raw_genotype, re.IGNORECASE):
                    # Format the genotype into A/A or A/G format from raw AA or AG format
                    if len(raw_genotype) == 1:
                        formatted_gt = f"{raw_genotype}/{raw_genotype}"
                    elif len(raw_genotype) == 2:
                        formatted_gt = f"{raw_genotype[0]}/{raw_genotype[1]}"
                    else:
                        formatted_gt = "Unknown/Unknown"
                    
                    self._debug_echo(f"Formatted genotype: {formatted_gt}")
                    return formatted_gt, formatted_gt  # For microarray files, these are the same
                else:
                    self._debug_echo(f"Invalid genotype format for {rsid}: {raw_genotype}")
                    return "0/0", "Unknown/Unknown"
            else:
                self._debug_echo(f"Could not find {rsid} in microarray data")
                print(f"Failed to read {rsid} from {self.input_file}: unknown file type or SNP not found")
                return "0/0", "Unknown/Unknown"

        # For VCF files, use bcftools to search
        # First, check if the VCF file contains this SNP by rsID (using exact word match)
        cmd = f"bcftools view -H {self.input_file} | grep -w {rsid} | wc -l"
        self._debug_echo(f"Command: {cmd}")
        snp_count_output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        snp_count = snp_count_output.stdout.strip()
        snp_count = re.sub(r'\s+', '', snp_count)

        self._debug_echo(f"snp_count: {snp_count}")

        if int(snp_count) > 0:
            self._debug_echo(f"Found {rsid} in VCF by ID")
            # Get the full line to extract both genotype and alleles
            cmd = f"bcftools view -H {self.input_file} | grep -w {rsid} | head -1"
            vcf_line_output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            vcf_line = vcf_line_output.stdout.strip()

            # Parse the VCF line
            fields = vcf_line.split("\t")
            if len(fields) >= 10:
                gt_field = fields[9].split(":")
                gt = gt_field[0] if gt_field else "0/0"
                ref = fields[3] if len(fields) > 3 else ""
                alt = fields[4] if len(fields) > 4 else ""

                self._debug_echo(f"ref={ref}, alt={alt}, gt={gt}, vcf_line={vcf_line}")

                # If genotype is missing or empty, set it to 0/0
                if not gt or gt == "./" or gt == "./.":
                    gt = "0/0"
                    actual_gt = f"{ref}/{ref}" if ref else "Unknown/Unknown"
                else:
                    # Convert 0/0, 0/1, 1/1 to actual alleles
                    if gt == "0/0":
                        actual_gt = f"{ref.upper()}/{ref.upper()}" if ref else "Unknown/Unknown"
                    elif gt == "1/1":
                        actual_gt = f"{alt.upper()}/{alt.upper()}" if alt else "Unknown/Unknown"
                    elif gt == "0/1" or gt == "1/0":
                        if ref and alt:
                            actual_gt = f"{ref.upper()}/{alt.upper()}"
                        else:
                            actual_gt = "Unknown/Unknown"
                    else:
                        # For multi-allelic or other complex cases
                        # Try to parse the genotype if possible
                        alleles = gt.split("/")
                        a1, a2 = alleles if len(alleles) == 2 else ["0", "0"]

                        a1_allele = ref if a1 == "0" and ref else (alt if a1 == "1" and alt else "Unknown")
                        a2_allele = ref if a2 == "0" and ref else (alt if a2 == "1" and alt else "Unknown")

                        if a1_allele != "Unknown" and a2_allele != "Unknown":
                            actual_gt = f"{a1_allele}/{a2_allele}"
                        else:
                            actual_gt = "Unknown/Unknown"

                self._debug_echo(f"Genotype: {gt} ({actual_gt})")
                return gt, actual_gt

        # Try by position if not found by rsID
        # Check if we need to convert chr notation
        position_cmd = ""
        if not chrpos.startswith("chr") and re.match(r'^[0-9]+:', chrpos):
            # Convert notation from "22:123456" to "chr22:123456"
            chr_num, pos = chrpos.split(":")
            chr_pos = f"chr{chr_num}:{pos}"
            self._debug_echo(f"Trying with chromosome notation: {chr_pos}")
            position_cmd = f"bcftools view -H -r {chr_pos} {self.input_file} 2>/dev/null | wc -l"
        else:
            position_cmd = f"bcftools view -H -r {chrpos} {self.input_file} 2>/dev/null | wc -l"

        position_count_output = subprocess.run(position_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        position_count = position_count_output.stdout.strip()
        position_count = re.sub(r'\s+', '', position_count)

        self._debug_echo(f"position_count: {position_count}")

        if int(position_count) > 0:
            # Similar process to get the genotype by position
            vcf_cmd = ""
            if not chrpos.startswith("chr") and re.match(r'^[0-9]+:', chrpos):
                chr_num, pos = chrpos.split(":")
                chr_pos = f"chr{chr_num}:{pos}"
                vcf_cmd = f"bcftools view -H -r {chr_pos} {self.input_file} | head -1"
                self._debug_echo(f"Found by position {chr_pos}")
            else:
                vcf_cmd = f"bcftools view -H -r {chrpos} {self.input_file} | head -1"
                self._debug_echo(f"Found by position {chrpos}")

            vcf_line_output = subprocess.run(vcf_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            vcf_line = vcf_line_output.stdout.strip()

            # Parse the VCF line (same as above)
            fields = vcf_line.split("\t")
            if len(fields) >= 10:
                gt_field = fields[9].split(":")
                gt = gt_field[0] if gt_field else "0/0"
                ref = fields[3] if len(fields) > 3 else ""
                alt = fields[4] if len(fields) > 4 else ""

                # Convert the genotype (same logic as above)
                if not gt or gt == "./" or gt == "./.":
                    gt = "0/0"
                    actual_gt = f"{ref}/{ref}" if ref else "Unknown/Unknown"
                else:
                    if gt == "0/0":
                        actual_gt = f"{ref.upper()}/{ref.upper()}" if ref else "Unknown/Unknown"
                    elif gt == "1/1":
                        actual_gt = f"{alt.upper()}/{alt.upper()}" if alt else "Unknown/Unknown"
                    elif gt == "0/1" or gt == "1/0":
                        if ref and alt:
                            actual_gt = f"{ref.upper()}/{alt.upper()}"
                        else:
                            actual_gt = "Unknown/Unknown"
                    else:
                        # Complex cases
                        alleles = gt.split("/")
                        a1, a2 = alleles if len(alleles) == 2 else ["0", "0"]

                        a1_allele = ref if a1 == "0" and ref else (alt if a1 == "1" and alt else "Unknown")
                        a2_allele = ref if a2 == "0" and ref else (alt if a2 == "1" and alt else "Unknown")

                        if a1_allele != "Unknown" and a2_allele != "Unknown":
                            actual_gt = f"{a1_allele}/{a2_allele}"
                        else:
                            actual_gt = "Unknown/Unknown"

                return gt, actual_gt

        # Not found by rsID or position
        return "0/0", "Unknown/Unknown"

    def check_snp(self, snp_data: str) -> None:
        """Check a single SNP and print its status.

        Args:
            snp_data: String with SNP data in format "rsid|chr:pos|protective_gt|risk_gt|interpretation"
        """
        parts = snp_data.split("|")
        if len(parts) < 5:
            print(f"Invalid SNP data format: {snp_data}")
            return

        rsid = parts[0]
        chrpos = parts[1]
        protective_gt = parts[2]
        risk_gt = parts[3]
        interpretation = parts[4]

        # Fetch the genotype
        gt, actual_gt = self._fetch_genotype(chrpos, rsid)
        
        # Prepare different genotype representations for matching
        if "/" in actual_gt and actual_gt != "Unknown/Unknown":
            # Get alleles
            g1, g2 = actual_gt.split("/")
            
            # -- Exactly match the shell script's approach --
            # Standard representation
            user_gt = actual_gt
            # Reversed representation
            rev_gt = f"{g2}/{g1}"
            # Complement representations
            comp_user_gt = f"{self._get_complement(g1)}/{self._get_complement(g2)}"
            comp_rev_gt = f"{self._get_complement(g2)}/{self._get_complement(g1)}"
            
            # Default is to display the user's genotype
            display_gt = user_gt
            status = "VARIANT"  # Default status
            color = self.colors["GRAY"]
            
            # -- First determine which representation to display --
            # Match the logic order in the shell script
            if user_gt == protective_gt or user_gt == risk_gt:
                display_gt = user_gt
            elif rev_gt == protective_gt or rev_gt == risk_gt:
                display_gt = rev_gt
            elif comp_user_gt == protective_gt or comp_user_gt == risk_gt:
                display_gt = comp_user_gt
            elif comp_rev_gt == protective_gt or comp_rev_gt == risk_gt:
                display_gt = comp_rev_gt
            
            # -- Then determine status based on the display genotype --
            if display_gt == protective_gt:
                status = "GOOD"
                color = self.colors["GREEN"]
            elif display_gt == risk_gt:
                status = "RISK"
                color = self.colors["RED"]
            elif g1 != g2:  # Heterozygote
                status = "CARRIER"
                color = self.colors["YELLOW"]
        else:
            # If we couldn't parse the genotype properly
            display_gt = actual_gt
            status = "UNKNOWN"
            color = self.colors["GRAY"]

        # Print result
        print(f"{rsid:<12} {display_gt:<12} {color}{status:<12}{self.colors['NC']} {interpretation}")
        
    def print_section_header(self, section: str) -> None:
        """Print a section header."""
        print("\n{:=^60}\n".format(f" {section} "))

    def _process_snp_list(self, section_name: str, snp_list: list):
        """Process a list of SNPs for a given section.
        
        Args:
            section_name: Name of the section being processed
            snp_list: List of SNPs to process
            
        Returns:
            List of formatted results for each SNP
        """
        results = []
        for snp in snp_list:
            try:
                # Format SNP data string
                snp_data = f"{snp['rs_id']}|{snp['position']}|{snp['ref_genotype']}|{snp['alt_genotype']}|{snp['description']}"
                
                # Get chromosome position and rsid from SNP data
                parts = snp_data.split('|')
                rsid = parts[0]
                chrpos = parts[1]
                ref_gt = parts[2]
                alt_gt = parts[3]
                interpretation = parts[4]
                
                # Get genotype for this SNP
                gt, actual_gt = self._fetch_genotype(chrpos, rsid)
                
                # Determine status based on genotype
                if gt == "Unknown/Unknown":
                    status = "UNKNOWN"
                    status_color = self.colors["GRAY"]
                elif gt == ref_gt:
                    status = "NORMAL"
                    status_color = self.colors["GREEN"]
                elif gt == alt_gt:
                    status = "RISK"
                    status_color = self.colors["RED"]
                else:
                    status = "CARRIER"
                    status_color = self.colors["YELLOW"]
                
                # Format the results as a tuple: (section_name, rsid, genotype, status, status_color, interpretation)
                results.append((section_name, rsid, gt, status, status_color, interpretation))
            except Exception as e:
                # Add error result if processing fails
                self._debug_echo(f"Error processing SNP {snp['rs_id']}: {str(e)}")
                results.append((section_name, snp['rs_id'], "Error", "ERROR", self.colors["GRAY"], str(e)))
                
        return results
    
    def _process_section(self, section_name: str):
        """Process all SNPs in a section.
        
        Args:
            section_name: Name of the section to process
            
        Returns:
            List of results for all SNPs in the section
        """
        if section_name not in self.catalog["categories"]:
            return []
            
        return self._process_snp_list(section_name, self.catalog["categories"][section_name])


        
    def run_analysis(self, section: str = "all"):
        """Run the SNP analysis for the specified section.

        Args:
            section: Section to analyze (all, metabolism, inflammation, etc.)
        """
        try:
            # Print column headers
            print(f"{'RSID':<12} {'Genotype':<12} {'Status':<12} {'Interpretation'}\n")
            
            # Get available sections from the catalog
            available_sections = self.catalog["categories"].keys()
            
            # Check if section is valid
            if section != "all" and section not in available_sections:
                print(f"Error: Unknown section '{section}'. Available sections: {', '.join(available_sections)}")
                return
            
            # Determine which sections to process
            sections_to_process = list(available_sections) if section == "all" else [section]
            
            # For VCF files, use parallel processing with threads
            if self.is_vcf:
                # Store results from each section
                results_lock = threading.Lock()
                section_results = {}
                
                # Define function to process a section in a thread
                def process_section_in_thread(section_name):
                    try:
                        result_list = []
                        for snp in self.catalog["categories"][section_name]:
                            # Format SNP data
                            snp_data = f"{snp['rs_id']}|{snp['position']}|{snp['ref_genotype']}|{snp['alt_genotype']}|{snp['description']}"
                            
                            # Parse SNP data
                            parts = snp_data.split('|')
                            rsid = parts[0]
                            chrpos = parts[1]
                            ref_gt = parts[2]
                            alt_gt = parts[3]
                            interpretation = parts[4]
                            
                            try:
                                # Get genotype for this SNP
                                gt, actual_gt = self._fetch_genotype(chrpos, rsid)
                                
                                # Process the genotype using the same logic as check_snp
                                # Prepare different genotype representations for matching
                                if actual_gt != "Unknown/Unknown":
                                    # Get alleles
                                    g1, g2 = actual_gt.split("/")
                                    
                                    # Standard representation
                                    user_gt = actual_gt
                                    # Reversed representation
                                    rev_gt = f"{g2}/{g1}"
                                    # Complement representations
                                    comp_user_gt = f"{self._get_complement(g1)}/{self._get_complement(g2)}"
                                    comp_rev_gt = f"{self._get_complement(g2)}/{self._get_complement(g1)}"
                                    
                                    # Default display and status
                                    display_gt = user_gt
                                    status = "CARRIER"  # Default status
                                    status_color = self.colors["YELLOW"]
                                    
                                    # Determine which representation to display (same as in check_snp)
                                    if user_gt == ref_gt or user_gt == alt_gt:
                                        display_gt = user_gt
                                    elif rev_gt == ref_gt or rev_gt == alt_gt:
                                        display_gt = rev_gt
                                    elif comp_user_gt == ref_gt or comp_user_gt == alt_gt:
                                        display_gt = comp_user_gt
                                    elif comp_rev_gt == ref_gt or comp_rev_gt == alt_gt:
                                        display_gt = comp_rev_gt
                                    
                                    # Determine status (same logic as in check_snp)
                                    for match_gt in [user_gt, rev_gt, comp_user_gt, comp_rev_gt]:
                                        if match_gt == ref_gt:
                                            status = "GOOD"
                                            status_color = self.colors["GREEN"]
                                            break
                                        elif match_gt == alt_gt:
                                            status = "RISK"
                                            status_color = self.colors["RED"]
                                            break
                                else:
                                    display_gt = "Unknown/Unknown"
                                    status = "UNKNOWN"
                                    status_color = self.colors["GRAY"]
                                
                                # Store result
                                result_list.append((rsid, display_gt, status, status_color, interpretation))
                            except Exception as e:
                                if self.debug:
                                    print(f"Error processing SNP {rsid}: {str(e)}")
                                result_list.append((rsid, "Error", "ERROR", self.colors["GRAY"], str(e)))
                        
                        # Store results for this section
                        with results_lock:
                            section_results[section_name] = result_list
                            
                    except Exception as e:
                        print(f"Error processing section {section_name}: {str(e)}")
                
                # Start a thread for each section
                max_workers = min(len(sections_to_process), multiprocessing.cpu_count() * 2)
                with ThreadPoolExecutor(max_workers=max_workers) as executor:
                    # Submit all sections for processing
                    futures = [executor.submit(process_section_in_thread, section_name) for section_name in sections_to_process]
                    
                    # Wait for all to complete
                    for future in as_completed(futures):
                        try:
                            future.result()  # Check if there were any exceptions
                        except Exception as e:
                            print(f"Error in thread: {str(e)}")
                
                # Print results section by section
                for section_name in sections_to_process:
                    if section_name in section_results and section_results[section_name]:
                        self.print_section_header(section_name.capitalize())
                        for result in section_results[section_name]:
                            # Unpack result: (rsid, genotype, status, status_color, interpretation)
                            rsid, genotype, status, status_color, interpretation = result
                            # Print formatted result
                            print(f"{rsid:<12} {genotype:<12} {status_color}{status:<12}{self.colors['NC']} {interpretation}")
                            
            else:  # For microarray files, use sequential processing (already fast)
                if section == "all":
                    for cat_name in available_sections:
                        self.print_section_header(cat_name.capitalize())
                        for snp in self.catalog["categories"][cat_name]:
                            snp_data = f"{snp['rs_id']}|{snp['position']}|{snp['ref_genotype']}|{snp['alt_genotype']}|{snp['description']}"
                            self.check_snp(snp_data)
                else:
                    self.print_section_header(section.capitalize())
                    for snp in self.catalog["categories"][section]:
                        snp_data = f"{snp['rs_id']}|{snp['position']}|{snp['ref_genotype']}|{snp['alt_genotype']}|{snp['description']}"
                        self.check_snp(snp_data)
        except Exception as e:
            print(f"Error during SNP analysis: {e}")
