#!/usr/bin/env python3
import sqlite3
import csv
import os
import argparse

def find_high_quality_clinvar_variants(db_path, output_csv=None, limit=100, explore=False):
    """Find ClinVar variants with high-quality review status"""
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    print(f"Connected to {db_path}")
    
    # Get total variants count
    cursor.execute("SELECT COUNT(*) FROM variant")
    total_variants = cursor.fetchone()[0]
    print(f"Total variants in database: {total_variants:,}")
    
    # First, explore the ClinVar data to see what's available
    if explore:
        print("\nExploring ClinVar data in the database...")
        
        # Check for variants with any ClinVar data
        cursor.execute("""
            SELECT COUNT(*) FROM variant 
            WHERE clinvar__sig IS NOT NULL
        """)
        clinvar_count = cursor.fetchone()[0]
        print(f"Variants with ClinVar significance: {clinvar_count:,}")
        
        # Sample some ClinVar significance values
        cursor.execute("""
            SELECT DISTINCT clinvar__sig FROM variant 
            WHERE clinvar__sig IS NOT NULL 
            LIMIT 10
        """)
        sig_samples = [row[0] for row in cursor.fetchall()]
        print(f"Sample ClinVar significance values: {sig_samples}")
        
        # Sample some review status values
        cursor.execute("""
            SELECT DISTINCT clinvar__rev_stat FROM variant 
            WHERE clinvar__rev_stat IS NOT NULL 
            LIMIT 10
        """)
        rev_samples = [row[0] for row in cursor.fetchall()]
        print(f"Sample review status values: {rev_samples}")
        
        # Count pathogenic variants
        cursor.execute("""
            SELECT COUNT(*) FROM variant 
            WHERE clinvar__sig LIKE '%pathogenic%' OR clinvar__sig LIKE '%Pathogenic%'
        """)
        pathogenic_count = cursor.fetchone()[0]
        print(f"Variants with pathogenic significance: {pathogenic_count:,}")
    
    # Run three separate queries instead of a combined query
    all_results = []
    
    # Level 1: Highest quality - pathogenic, likely pathogenic, or risk factor with expert review or practice guidelines
    print("Searching for highest quality variants (expert panel review)...")
    query1 = """
    SELECT 
        base__chrom, base__pos, base__ref_base, base__alt_base, 
        base__hugo, base__so, 
        clinvar__sig, clinvar__disease_names, clinvar__rev_stat,
        clinvar__id, dbsnp__rsid, gnomad__af
    FROM variant
    WHERE 
        clinvar__sig IS NOT NULL
        AND (
            clinvar__sig LIKE '%pathogenic%'
            OR clinvar__sig LIKE '%Pathogenic%'
            OR clinvar__sig LIKE '%risk factor%'
            OR clinvar__sig LIKE '%protective%'
        )
        AND (
            clinvar__rev_stat LIKE '%expert panel%'
            OR clinvar__rev_stat LIKE '%practice guideline%'
        )
    LIMIT ?
    """
    cursor.execute(query1, (limit//3 + 1,))
    for row in cursor.fetchall():
        all_results.append(row + ('Highest',))  # Add quality level as last element
    
    highest_count = len(all_results)
    print(f"Found {highest_count} variants with expert panel review")
    
    # Level 2: High quality - pathogenic or clinically significant with multiple submitters, no conflicts
    print("Searching for high quality variants (multiple submitters, no conflicts)...")
    query2 = """
    SELECT 
        base__chrom, base__pos, base__ref_base, base__alt_base, 
        base__hugo, base__so, 
        clinvar__sig, clinvar__disease_names, clinvar__rev_stat,
        clinvar__id, dbsnp__rsid, gnomad__af
    FROM variant
    WHERE 
        clinvar__sig IS NOT NULL
        AND (
            clinvar__sig LIKE '%pathogenic%'
            OR clinvar__sig LIKE '%Pathogenic%'
            OR clinvar__sig LIKE '%risk factor%'
            OR clinvar__sig = 'drug response'
            OR clinvar__sig = 'association'
        )
        AND clinvar__rev_stat LIKE '%multiple submitters, no conflicts%'
    LIMIT ?
    """
    cursor.execute(query2, (limit//3 + 1,))
    for row in cursor.fetchall():
        all_results.append(row + ('High',))  # Add quality level as last element
    
    high_count = len(all_results) - highest_count
    print(f"Found {high_count} variants with multiple submitters, no conflicts")
    
    # Level 3: Moderate quality - pathogenic with single submitter criteria provided
    print("Searching for moderate quality variants (single submitter)...")
    query3 = """
    SELECT 
        base__chrom, base__pos, base__ref_base, base__alt_base, 
        base__hugo, base__so, 
        clinvar__sig, clinvar__disease_names, clinvar__rev_stat,
        clinvar__id, dbsnp__rsid, gnomad__af
    FROM variant
    WHERE 
        clinvar__sig IS NOT NULL
        AND (
            clinvar__sig LIKE '%pathogenic%'
            OR clinvar__sig LIKE '%Pathogenic%'
            OR clinvar__sig LIKE '%drug response%'
        )
        AND clinvar__rev_stat LIKE '%single submitter%'
    LIMIT ?
    """
    cursor.execute(query3, (limit//3 + 1,))
    for row in cursor.fetchall():
        all_results.append(row + ('Moderate',))  # Add quality level as last element
    
    moderate_count = len(all_results) - highest_count - high_count
    print(f"Found {moderate_count} variants with single submitter criteria")
    
    # Level 4: Interesting - pathogenic with conflicting classifications but still of interest
    print("Searching for interesting variants with conflicting classifications...")
    query4 = """
    SELECT 
        base__chrom, base__pos, base__ref_base, base__alt_base, 
        base__hugo, base__so, 
        clinvar__sig, clinvar__disease_names, clinvar__rev_stat,
        clinvar__id, dbsnp__rsid, gnomad__af
    FROM variant
    WHERE 
        clinvar__sig IS NOT NULL
        AND clinvar__sig LIKE '%Conflicting classifications%'
        AND clinvar__rev_stat LIKE '%criteria provided%'
        AND base__hugo IS NOT NULL
    LIMIT ?
    """
    cursor.execute(query4, (limit//3 + 1,))
    for row in cursor.fetchall():
        all_results.append(row + ('Interest',))  # Add quality level as last element
    
    interest_count = len(all_results) - highest_count - high_count - moderate_count
    print(f"Found {interest_count} interesting variants with conflicting classifications")
    
    # Sort results by quality level and gene name
    def sort_key(row):
        quality = row[12]  # Quality level is the last element (index 12)
        gene = row[4] or ""  # Gene name is element 4
        quality_rank = {'Highest': 1, 'High': 2, 'Moderate': 3, 'Interest': 4}.get(quality, 5)
        return (quality_rank, gene)
    
    all_results.sort(key=sort_key)
    
    # Limit to requested number
    results = all_results[:limit]
    
    print(f"\nFound {len(results)} high-quality ClinVar variants")
    
    if results:
        # Print header
        header = "Quality Chrom Pos Ref>Alt Gene Effect ClinVar_Significance Disease Review_Status ClinVar_ID rsID gnomAD_AF"
        print("\n" + header)
        print("-" * len(header))
        
        # Print results
        for row in results:
            chrom, pos, ref, alt, gene, effect, sig, disease, rev_stat, clinvar_id, rsid, af, quality = row
            
            # Format for display
            disease_str = str(disease)[:30] + "..." if disease and len(str(disease)) > 30 else (disease or "Unknown")
            sig_str = str(sig)[:30] if sig else "Unknown"
            rsid_str = rsid or "None"
            af_str = f"{af:.6f}" if af is not None else "Unknown"
            
            print(f"{quality:<8} {chrom} {pos} {ref}>{alt} {gene or 'Unknown'} {effect or 'Unknown'} {sig_str} {disease_str} {rev_stat or 'Unknown'} {clinvar_id or 'None'} {rsid_str} {af_str}")
            
        # Export to CSV if requested
        if output_csv:
            with open(output_csv, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                
                # Write header
                writer.writerow([
                    "Quality_Level", "Chromosome", "Position", "Reference", "Alternate", 
                    "Gene", "Effect", "ClinVar_Significance", "Disease", 
                    "Review_Status", "ClinVar_ID", "rsID", "gnomAD_AF"
                ])
                
                # Write data
                for row in results:
                    writer.writerow(row)
                    
                print(f"\nResults exported to {output_csv}")
    
    conn.close()
    return results

def main():
    parser = argparse.ArgumentParser(description="Find high-quality ClinVar variants in OpenCravat database")
    parser.add_argument("db_path", help="Path to the OpenCravat SQLite database")
    parser.add_argument("--csv", "-c", help="Output CSV file path", default=None)
    parser.add_argument("--limit", "-l", type=int, default=100, help="Limit number of variants to return")
    parser.add_argument("--explore", "-e", action="store_true", help="Explore ClinVar data in the database")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.db_path):
        print(f"Error: Database file {args.db_path} not found")
        return 1
        
    find_high_quality_clinvar_variants(args.db_path, args.csv, args.limit, args.explore)
    return 0

if __name__ == "__main__":
    main()
