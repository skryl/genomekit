#!/usr/bin/env python3
import sqlite3
import argparse
import csv
import sys

def explore_clinvar_data(db_path):
    """Explore ClinVar data in the database to get an overview of the review statuses available"""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Check for total number of variants
    cursor.execute("SELECT COUNT(*) FROM variant")
    total_variants = cursor.fetchone()[0]
    print(f"Total variants in database: {total_variants:,}")
    
    # Check for ClinVar data
    cursor.execute("SELECT COUNT(*) FROM variant WHERE clinvar__sig IS NOT NULL")
    clinvar_count = cursor.fetchone()[0]
    print(f"Variants with ClinVar significance: {clinvar_count:,}")
    
    # Check for review status distribution
    cursor.execute("""
    SELECT clinvar__rev_stat, COUNT(*) 
    FROM variant 
    WHERE clinvar__rev_stat IS NOT NULL 
    GROUP BY clinvar__rev_stat 
    ORDER BY COUNT(*) DESC
    """)
    
    print("\nReview status distribution:")
    for row in cursor.fetchall():
        status, count = row
        print(f"{status}: {count:,}")
    
    # Check for 3-4 star ratings
    cursor.execute("""
    SELECT COUNT(*) 
    FROM variant 
    WHERE (
        clinvar__rev_stat LIKE '%expert panel%' 
        OR clinvar__rev_stat LIKE '%practice guideline%'
    )
    """)
    expert_count = cursor.fetchone()[0]
    print(f"\nVariants with expert panel or practice guideline review (3-4 stars): {expert_count:,}")
    
    conn.close()

def find_expert_clinvar_variants(db_path, output_csv=None):
    """Find all ClinVar variants with expert review (3-4 stars) in the database"""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    print(f"Connected to {db_path}")
    
    # Get column names from the variant table to see what we have available
    cursor.execute("PRAGMA table_info(variant)")
    columns = [row[1] for row in cursor.fetchall()]
    print(f"Available columns in the variant table: {len(columns)}")
    
    # Find all ClinVar variants with expert panel or practice guideline review
    query = """
    SELECT 
        base__chrom, base__pos, base__ref_base, base__alt_base, 
        base__hugo, base__so, 
        clinvar__sig, clinvar__disease_names, clinvar__rev_stat,
        clinvar__id, dbsnp__rsid, gnomad__af, vcfinfo__zygosity
    FROM variant
    WHERE 
        (clinvar__rev_stat LIKE '%expert panel%' 
        OR clinvar__rev_stat LIKE '%practice guideline%')
    ORDER BY 
        CASE
            WHEN clinvar__sig LIKE '%pathogenic%' THEN 1
            WHEN clinvar__sig LIKE '%Pathogenic%' THEN 1
            WHEN clinvar__sig LIKE '%risk factor%' THEN 2
            WHEN clinvar__sig LIKE '%drug response%' THEN 3
            ELSE 4
        END,
        base__chrom, base__pos
    """
    
    cursor.execute(query)
    expert_variants = cursor.fetchall()
    
    print(f"\nFound {len(expert_variants)} ClinVar variants with expert review (3-4 stars)")
    
    # Filter for variants that the user actually carries (non-reference genotype)
    carried_variants = [v for v in expert_variants if v[12] and v[12] not in ['0/0', '0|0', './.', '.|.']]
    print(f"Of these, you carry {len(carried_variants)} variants (heterozygous or homozygous)")
    
    if carried_variants:
        # Print header
        header = "Chrom Pos Ref>Alt Gene Effect ClinVar_Significance Disease Review_Status Stars ClinVar_ID rsID gnomAD_AF Zygosity"
        print("\n" + header)
        print("-" * len(header))
        
        # Print results
        for row in carried_variants:
            chrom, pos, ref, alt, gene, effect, sig, disease, rev_stat, clinvar_id, rsid, af, zygosity = row
            
            # Calculate star rating
            stars = ""
            if rev_stat and "practice guideline" in rev_stat.lower():
                stars = "★★★★"
            elif rev_stat and "expert panel" in rev_stat.lower():
                stars = "★★★"
            
            # Format for display
            disease_str = str(disease)[:30] + "..." if disease and len(str(disease)) > 30 else (disease or "Unknown")
            sig_str = str(sig)[:30] + "..." if sig and len(str(sig)) > 30 else (sig or "Unknown")
            rsid_str = rsid or "None"
            af_str = f"{af:.6f}" if af is not None else "Unknown"
            
            print(f"{chrom} {pos} {ref}>{alt} {gene or 'Unknown'} {effect or 'Unknown'} {sig_str} {disease_str} {rev_stat or 'Unknown'} {stars} {clinvar_id or 'None'} {rsid_str} {af_str} {zygosity or 'Unknown'}")
        
        # Export to CSV if requested
        if output_csv:
            with open(output_csv, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                
                # Write header
                writer.writerow([
                    "Chromosome", "Position", "Reference", "Alternate", 
                    "Gene", "Effect", "ClinVar_Significance", "Disease", 
                    "Review_Status", "Star_Rating", "ClinVar_ID", "rsID", "gnomAD_AF", "Zygosity"
                ])
                
                # Write data
                for row in carried_variants:
                    chrom, pos, ref, alt, gene, effect, sig, disease, rev_stat, clinvar_id, rsid, af, zygosity = row
                    
                    # Calculate star rating
                    star_rating = ""
                    if rev_stat and "practice guideline" in rev_stat.lower():
                        star_rating = "4"
                    elif rev_stat and "expert panel" in rev_stat.lower():
                        star_rating = "3"
                    
                    writer.writerow([
                        chrom, pos, ref, alt, gene, effect, sig, disease, rev_stat, 
                        star_rating, clinvar_id, rsid, af, zygosity
                    ])
                
                print(f"\nResults exported to {output_csv}")
    
    conn.close()

def main():
    parser = argparse.ArgumentParser(description='Find ClinVar variants with expert review in your genome')
    parser.add_argument('db_path', help='Path to the SQLite database')
    parser.add_argument('--csv', help='Output CSV file path')
    parser.add_argument('--explore', action='store_true', help='Explore ClinVar data in the database')
    
    args = parser.parse_args()
    
    if args.explore:
        explore_clinvar_data(args.db_path)
    else:
        find_expert_clinvar_variants(args.db_path, args.csv)

if __name__ == '__main__':
    main()
