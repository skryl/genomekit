#!/usr/bin/env python3
import sqlite3
import csv
import os
import argparse

class ClinicalVariantAnalyzer:
    """Analyze and prioritize variants based on clinical annotations"""

    def __init__(self, db_path):
        """Initialize with path to SQLite database"""
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.tables = self._get_tables()
        print(f"Connected to {db_path}")
        print(f"Available tables: {', '.join(self.tables)}")

        # For tracking total variants
        self.total_variants = self._count_variants()
        print(f"Total variants in database: {self.total_variants:,}")

    def _get_tables(self):
        """Get all tables in the database"""
        cursor = self.conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = [row[0] for row in cursor.fetchall()]
        cursor.close()
        return tables

    def _count_variants(self):
        """Count total variants in the database"""
        cursor = self.conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM variant")
        count = cursor.fetchone()[0]
        cursor.close()
        return count

    def get_table_schema(self, table_name):
        """Get schema for a specific table"""
        cursor = self.conn.cursor()
        cursor.execute(f"PRAGMA table_info({table_name})")
        schema = cursor.fetchall()
        cursor.close()

        print(f"\n==== Schema for {table_name} ====")
        for col in schema:
            print(f"{col[0]}: {col[1]} ({col[2]})")

        return schema

    def find_clinvar_pathogenic(self, limit=100):
        """Find ClinVar pathogenic or likely pathogenic variants"""
        try:
            # Check if variant_annotator table exists
            if "variant_annotator" not in self.tables:
                print("variant_annotator table not found.")
                return []
                
            # First, get the schema of variant_annotator to see which columns exist
            cursor = self.conn.cursor()
            cursor.execute("PRAGMA table_info(variant_annotator)")
            cols = [col[1] for col in cursor.fetchall()]
            
            # Check for ClinVar columns
            clinvar_cols = [col for col in cols if "clinvar" in col.lower()]
            if not clinvar_cols:
                print("No ClinVar columns found in variant_annotator table.")
                print(f"Available columns: {', '.join(cols[:20])}...")
                return []
                
            # Find columns for clinical significance and review status
            sig_col = next((col for col in clinvar_cols if "sig" in col.lower()), None)
            review_col = next((col for col in clinvar_cols if "review" in col.lower()), None)
            star_col = next((col for col in clinvar_cols if "star" in col.lower()), None)
            disease_col = next((col for col in clinvar_cols if "disease" in col.lower() or "phenotype" in col.lower()), None)
            
            print(f"Using ClinVar columns: {sig_col}, {review_col}, {star_col}, {disease_col}")
            
            # Build the query based on available columns
            base_query = """
            SELECT v.base_uid, v.chrom, v.pos, v.ref_base, v.alt_base,
                   v.hugo, v.so, """
                   
            # Add CADD if available
            cadd_col = next((col for col in cols if "cadd" in col.lower() and "phred" in col.lower()), None)
            if cadd_col:
                base_query += f"va.{cadd_col}, "
            else:
                base_query += "NULL as cadd_phred, "
                
            # Add ClinVar columns
            if sig_col:
                base_query += f"va.{sig_col}, "
            else:
                base_query += "NULL as clinvar_sig, "
                
            if review_col:
                base_query += f"va.{review_col}, "
            else:
                base_query += "NULL as review_status, "
                
            if star_col:
                base_query += f"va.{star_col}, "
            else:
                base_query += "NULL as star_rating, "
                
            if disease_col:
                base_query += f"va.{disease_col}"
            else:
                base_query += "NULL as disease"
                
            query = base_query + """
            FROM variant v
            JOIN variant_annotator va ON v.base_uid = va.base_uid
            WHERE 1=1
            """
            
            # Add filter for pathogenic variants if sig_col exists
            if sig_col:
                query += f"""
                AND (va.{sig_col} LIKE '%pathogenic%'
                   OR va.{sig_col} LIKE '%Pathogenic%'
                   OR va.{sig_col} LIKE '%likely_pathogenic%'
                   OR va.{sig_col} LIKE '%Likely_pathogenic%')
                """
                
            # Add filter for high-quality evidence if review_col or star_col exists
            if review_col:
                query += f"""
                AND (va.{review_col} LIKE '%expert%' 
                   OR va.{review_col} LIKE '%panel%'
                   OR va.{review_col} LIKE '%practice%guideline%')
                """
                
            if star_col:
                query += f"""
                AND (va.{star_col} >= 3 OR va.{star_col} = '3' OR va.{star_col} = '4')
                """
                
            # Add ordering and limit
            if cadd_col:
                query += f"ORDER BY va.{cadd_col} DESC "
            query += "LIMIT ?"

            cursor = self.conn.cursor()
            cursor.execute(query, (limit,))
            results = cursor.fetchall()
            cursor.close()

            print(f"\n==== High-Quality ClinVar Pathogenic Variants ({len(results)}) ====\n")

            if results:
                # Print header
                header_parts = ["Chrom", "Position", "Ref>Alt", "Gene", "Effect"]
                if cadd_col:
                    header_parts.append("CADD")
                if sig_col:
                    header_parts.append("ClinVar Sig")
                if review_col:
                    header_parts.append("Review Status")
                if star_col:
                    header_parts.append("Stars")
                if disease_col:
                    header_parts.append("Disease")
                    
                # Format header
                header = ""
                for part in header_parts:
                    width = 30 if part == "Disease" else 20 if part == "Review Status" else 15 if part == "ClinVar Sig" else 6 if part == "CADD" or part == "Stars" else 10 if part == "Position" else 8 if part == "Ref>Alt" else 15
                    header += f"{part:<{width}} "
                print(header)
                print("-" * (len(header) + 10))

                # Print each result
                for row in results:
                    formatted_row = ""
                    # Basic variant info (always present)
                    uid, chrom, pos, ref, alt, gene, effect = row[0:7]
                    formatted_row += f"{chrom:<6} {pos:<10} {ref}>{alt:<6} {gene or 'Unknown':<15} {effect or 'Unknown':<15} "
                    
                    # Additional annotations (may or may not be present)
                    col_index = 7
                    if cadd_col:
                        cadd = row[col_index]
                        formatted_row += f"{cadd or 'N/A':<6} "
                        col_index += 1
                    if sig_col:
                        sig = row[col_index]
                        formatted_row += f"{str(sig)[:15] if sig else 'Unknown':<15} "
                        col_index += 1
                    if review_col:
                        review = row[col_index]
                        formatted_row += f"{str(review)[:20] if review else 'Unknown':<20} "
                        col_index += 1
                    if star_col:
                        stars = row[col_index]
                        formatted_row += f"{stars or 'N/A':<6} "
                        col_index += 1
                    if disease_col:
                        disease = row[col_index]
                        formatted_row += f"{str(disease)[:30] if disease else 'Unknown':<30}"
                        
                    print(formatted_row)

            return results
        except Exception as e:
            print(f"Error finding ClinVar pathogenic variants: {e}")
            import traceback
            traceback.print_exc()
            return []

    def find_rare_high_impact(self, gnomad_af_threshold=0.01, cadd_threshold=20, limit=100):
        """Find rare variants with high predicted impact"""
        try:
            query = """
            SELECT v.base__uid, v.base__chrom, v.base__pos, v.base__ref_base, v.base__alt_base,
                   v.base__hugo, v.base__so, v.cadd__phred,
                   g.gnomad__af, g.gnomad__af_popmax
            FROM variant v
            JOIN gnomad g ON v.base__uid = g.base__uid
            WHERE (g.gnomad__af < ? OR g.gnomad__af IS NULL)
              AND v.cadd__phred > ?
            ORDER BY v.cadd__phred DESC
            LIMIT ?
            """

            cursor = self.conn.cursor()
            cursor.execute(query, (gnomad_af_threshold, cadd_threshold, limit))
            results = cursor.fetchall()
            cursor.close()

            print(f"\n==== Rare High-Impact Variants (AF < {gnomad_af_threshold}, CADD > {cadd_threshold}) ({len(results)}) ====")

            if results:
                # Print header
                print(f"{'Chrom':<6} {'Position':<10} {'Ref>Alt':<8} {'Gene':<15} {'Effect':<20} {'CADD':<6} {'gnomAD AF':<10}")
                print("-" * 90)

                for row in results:
                    uid, chrom, pos, ref, alt, gene, effect, cadd, af, af_popmax = row
                    af_str = f"{af:.6f}" if af is not None else "Unknown"

                    print(f"{chrom:<6} {pos:<10} {ref}>{alt:<6} {gene or 'Unknown':<15} {effect or 'Unknown':<20} {cadd or 'N/A':<6} {af_str:<10}")

            return results
        except Exception as e:
            print(f"Error finding rare high-impact variants: {e}")
            if "gnomad" not in self.tables:
                print("gnomAD table not found. Make sure gnomAD annotator was run.")
            return []

    def find_cancer_variants(self, limit=100):
        """Find variants associated with cancer in COSMIC or CancerHotspots"""
        try:
            # First try query with both COSMIC and CancerHotspots
            query = """
            SELECT v.base__uid, v.base__chrom, v.base__pos, v.base__ref_base, v.base__alt_base,
                   v.base__hugo, v.base__so, v.cadd__phred,
                   c.cosmic__cosmic_id, c.cosmic__primary_site,
                   ch.cancerhotspots__tissue, ch.cancerhotspots__q_value
            FROM variant v
            LEFT JOIN cosmic c ON v.base__uid = c.base__uid
            LEFT JOIN cancerhotspots ch ON v.base__uid = ch.base__uid
            WHERE c.cosmic__cosmic_id IS NOT NULL OR ch.cancerhotspots__q_value IS NOT NULL
            ORDER BY v.cadd__phred DESC
            LIMIT ?
            """

            cursor = self.conn.cursor()
            cursor.execute(query, (limit,))
            results = cursor.fetchall()
            cursor.close()

            # If no results, try just with COSMIC
            if not results and "cosmic" in self.tables:
                query = """
                SELECT v.base__uid, v.base__chrom, v.base__pos, v.base__ref_base, v.base__alt_base,
                       v.base__hugo, v.base__so, v.cadd__phred,
                       c.cosmic__cosmic_id, c.cosmic__primary_site,
                       NULL, NULL
                FROM variant v
                JOIN cosmic c ON v.base__uid = c.base__uid
                ORDER BY v.cadd__phred DESC
                LIMIT ?
                """

                cursor = self.conn.cursor()
                cursor.execute(query, (limit,))
                results = cursor.fetchall()
                cursor.close()

            print(f"\n==== Cancer-Associated Variants ({len(results)}) ====")

            if results:
                # Print header
                print(f"{'Chrom':<6} {'Position':<10} {'Ref>Alt':<8} {'Gene':<15} {'CADD':<6} {'COSMIC ID':<15} {'Cancer Site':<20}")
                print("-" * 90)

                for row in results:
                    uid, chrom, pos, ref, alt, gene, effect, cadd, cosmic_id, site, tissue, q_value = row
                    cancer_site = site or tissue or "Unknown"

                    print(f"{chrom:<6} {pos:<10} {ref}>{alt:<6} {gene or 'Unknown':<15} {cadd or 'N/A':<6} {cosmic_id or 'N/A':<15} {cancer_site[:20]}")

            return results
        except Exception as e:
            print(f"Error finding cancer variants: {e}")
            if "cosmic" not in self.tables and "cancerhotspots" not in self.tables:
                print("Neither COSMIC nor CancerHotspots tables found.")
            return []

    def find_omim_disease_variants(self, limit=100):
        """Find variants in genes associated with OMIM diseases"""
        try:
            query = """
            SELECT v.base__uid, v.base__chrom, v.base__pos, v.base__ref_base, v.base__alt_base,
                   v.base__hugo, v.base__so, v.cadd__phred,
                   o.omim__omim_id, o.omim__phenotypes, o.omim__inheritance
            FROM variant v
            JOIN omim o ON v.base__uid = o.base__uid
            WHERE o.omim__phenotypes IS NOT NULL
            ORDER BY v.cadd__phred DESC
            LIMIT ?
            """

            cursor = self.conn.cursor()
            cursor.execute(query, (limit,))
            results = cursor.fetchall()
            cursor.close()

            print(f"\n==== OMIM Disease-Associated Variants ({len(results)}) ====")

            if results:
                # Print header
                print(f"{'Chrom':<6} {'Position':<10} {'Ref>Alt':<8} {'Gene':<15} {'CADD':<6} {'OMIM ID':<10} {'Inheritance':<12} {'Phenotype':<30}")
                print("-" * 100)

                for row in results:
                    uid, chrom, pos, ref, alt, gene, effect, cadd, omim_id, phenotypes, inheritance = row
                    phenotype_str = str(phenotypes)[:30] if phenotypes else "Unknown"
                    inheritance_str = inheritance or "Unknown"

                    print(f"{chrom:<6} {pos:<10} {ref}>{alt:<6} {gene or 'Unknown':<15} {cadd or 'N/A':<6} {omim_id or 'N/A':<10} {inheritance_str:<12} {phenotype_str}")

            return results
        except Exception as e:
            print(f"Error finding OMIM disease variants: {e}")
            if "omim" not in self.tables:
                print("OMIM table not found. Make sure OMIM annotator was run.")
            return []

    def export_medical_prioritized_variants(self, output_file, limit=1000):
        """Export prioritized list of medically significant variants"""
        def score_variant(var_data):
            """Score variant based on clinical significance"""
            score = 100  # Base score

            # CADD score contribution (0-100 points)
            cadd = var_data[7]
            if cadd is not None:
                try:
                    cadd_value = float(cadd)
                    if cadd_value >= 30:  # Extremely damaging
                        score += 100
                    elif cadd_value >= 25:
                        score += 80
                    elif cadd_value >= 20:
                        score += 60
                    elif cadd_value >= 15:
                        score += 40
                    elif cadd_value < 10:
                        score -= 50
                except (ValueError, TypeError):
                    pass

            # ClinVar contribution (-50 to +100 points)
            clinvar_sig = var_data[8] if len(var_data) > 8 else None
            if clinvar_sig:
                if "pathogenic" in str(clinvar_sig).lower() and "likely" not in str(clinvar_sig).lower():
                    score += 100
                elif "likely_pathogenic" in str(clinvar_sig).lower():
                    score += 50
                elif "uncertain_significance" in str(clinvar_sig).lower() or "vus" in str(clinvar_sig).lower():
                    score += 30
                elif "benign" in str(clinvar_sig).lower() and "likely" not in str(clinvar_sig).lower():
                    score -= 50
                elif "likely_benign" in str(clinvar_sig).lower():
                    score -= 30

            # gnomAD frequency contribution (-50 to 0 points)
            if len(var_data) > 10 and var_data[10] is not None:
                try:
                    af = float(var_data[10])
                    if af < 0.0001:  # Very rare
                        score += 50
                    elif af < 0.001:
                        score += 30
                    elif af < 0.01:
                        score += 10
                    elif af > 0.05:
                        score -= 20
                except (ValueError, TypeError):
                    pass

            # COSMIC/Cancer contribution (0-50 points)
            if len(var_data) > 12 and var_data[12] is not None:
                score += 50  # Present in COSMIC

            # OMIM contribution (0-50 points)
            if len(var_data) > 14 and var_data[14] is not None:
                score += 50  # Associated with OMIM disease

            return score

        try:
            # Get schema for variant_annotator to find relevant columns
            cursor = self.conn.cursor()
            cursor.execute("PRAGMA table_info(variant_annotator)")
            va_cols = [col[1] for col in cursor.fetchall()]

            # Find clinical annotation columns
            cadd_col = next((col for col in va_cols if "cadd" in col.lower() and "phred" in col.lower()), None)
            clinvar_sig_col = next((col for col in va_cols if "clinvar" in col.lower() and "sig" in col.lower()), None)
            clinvar_disease_col = next((col for col in va_cols if "clinvar" in col.lower() and ("disease" in col.lower() or "phenotype" in col.lower())), None)
            gnomad_af_col = next((col for col in va_cols if "gnomad" in col.lower() and "af" in col.lower() and "popmax" not in col.lower()), None)
            gnomad_popmax_col = next((col for col in va_cols if "gnomad" in col.lower() and "popmax" in col.lower()), None)
            cosmic_id_col = next((col for col in va_cols if "cosmic" in col.lower() and "id" in col.lower()), None)
            cosmic_site_col = next((col for col in va_cols if "cosmic" in col.lower() and ("site" in col.lower() or "tissue" in col.lower())), None)
            omim_id_col = next((col for col in va_cols if "omim" in col.lower() and "id" in col.lower()), None)
            omim_phenotype_col = next((col for col in va_cols if "omim" in col.lower() and "phenotype" in col.lower()), None)
            omim_inheritance_col = next((col for col in va_cols if "omim" in col.lower() and "inheritance" in col.lower()), None)
            
            print(f"Using clinical columns from variant_annotator:")
            print(f"CADD: {cadd_col}")
            print(f"ClinVar: {clinvar_sig_col}, {clinvar_disease_col}")
            print(f"gnomAD: {gnomad_af_col}, {gnomad_popmax_col}")
            print(f"COSMIC: {cosmic_id_col}, {cosmic_site_col}")
            print(f"OMIM: {omim_id_col}, {omim_phenotype_col}, {omim_inheritance_col}")
            
            # Build SQL query dynamically based on available columns
            select_parts = ["v.base_uid", "v.chrom", "v.pos", "v.ref_base", "v.alt_base", "v.hugo", "v.so"]
            
            if cadd_col:
                select_parts.append(f"va.{cadd_col}")
            else:
                select_parts.append("NULL as cadd_phred")
                
            if clinvar_sig_col:
                select_parts.append(f"va.{clinvar_sig_col}")
            else:
                select_parts.append("NULL as clinvar_sig")
                
            if clinvar_disease_col:
                select_parts.append(f"va.{clinvar_disease_col}")
            else:
                select_parts.append("NULL as clinvar_disease")
                
            if gnomad_af_col:
                select_parts.append(f"va.{gnomad_af_col}")
            else:
                select_parts.append("NULL as gnomad_af")
                
            if gnomad_popmax_col:
                select_parts.append(f"va.{gnomad_popmax_col}")
            else:
                select_parts.append("NULL as gnomad_af_popmax")
                
            if cosmic_id_col:
                select_parts.append(f"va.{cosmic_id_col}")
            else:
                select_parts.append("NULL as cosmic_id")
                
            if cosmic_site_col:
                select_parts.append(f"va.{cosmic_site_col}")
            else:
                select_parts.append("NULL as cosmic_site")
                
            if omim_id_col:
                select_parts.append(f"va.{omim_id_col}")
            else:
                select_parts.append("NULL as omim_id")
                
            if omim_phenotype_col:
                select_parts.append(f"va.{omim_phenotype_col}")
            else:
                select_parts.append("NULL as omim_phenotype")
                
            if omim_inheritance_col:
                select_parts.append(f"va.{omim_inheritance_col}")
            else:
                select_parts.append("NULL as omim_inheritance")
            
            # Create the full query
            query = f"""
            SELECT {', '.join(select_parts)}
            FROM variant v
            JOIN variant_annotator va ON v.base_uid = va.base_uid
            WHERE 1=1
            """
            # Only include variants with clinical significance
            filter_conditions = []
            
            if clinvar_sig_col:
                filter_conditions.append(f"va.{clinvar_sig_col} IS NOT NULL")
                filter_conditions.append(f"va.{clinvar_sig_col} LIKE '%pathogenic%'")
                
            if gnomad_af_col:
                filter_conditions.append(f"va.{gnomad_af_col} < 0.01")
                
            if cosmic_id_col:
                filter_conditions.append(f"va.{cosmic_id_col} IS NOT NULL")
                
            if omim_id_col:
                filter_conditions.append(f"va.{omim_id_col} IS NOT NULL")
                
            # Exclude non-interesting variant types
            filter_conditions.append("(v.so != 'synonymous_variant' AND v.so != 'intron_variant')")
            
            # Add the WHERE clause if we have conditions
            if filter_conditions:
                query += f"AND ({' OR '.join(filter_conditions)})\n"
                
            query += f"LIMIT {limit}"
            
            print("\nExecuting query:\n", query)
            
            cursor = self.conn.cursor()
            cursor.execute(query)
            all_results = cursor.fetchall()
            cursor.close()

            print(f"\nFound {len(all_results)} potentially medically significant variants")

            # Score and sort variants
            scored_variants = [(row, score_variant(row)) for row in all_results]
            scored_variants.sort(key=lambda x: x[1], reverse=True)

            top_variants = scored_variants[:limit]  # Take top variants for export

            # Export to CSV
            with open(output_file, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)

                # Write header
                writer.writerow([
                    "PriorityScore", "Chromosome", "Position", "Ref", "Alt", "Gene", "Effect", "CADD",
                    "ClinVar", "ClinVar Disease", "gnomAD AF", "gnomAD AF PopMax", "COSMIC ID", "Cancer Site", "OMIM ID", "OMIM Phenotype", "Inheritance"
                ])

                # Write data
                for var, score in top_variants:
                    # Extract the fields (variable number based on what columns were available)
                    row_data = [round(score, 1)]
                    
                    # Basic variant info (always present)
                    for i in range(7):
                        row_data.append(var[i] if i < len(var) else None)
                    
                    # Additional annotations (may not all be present)
                    for i in range(7, 17):
                        row_data.append(var[i] if i < len(var) else None)
                        
                    writer.writerow(row_data)

            print(f"\nExported prioritized variants to {output_file}")
            print("\nTop 20 highest-priority variants:")
            print(f"{'Score':<7} {'Chrom':<6} {'Position':<10} {'Ref>Alt':<8} {'Gene':<15} {'ClinVar':<15} {'OMIM':<8} {'COSMIC':<8}")
            print("-" * 85)

            for var, score in top_variants[:20]:
                # Basic variant info
                uid, chrom, pos, ref, alt, gene, effect = [var[i] if i < len(var) else None for i in range(7)]
                
                # Get clinical annotations if available
                clinvar_sig = var[8] if len(var) > 8 and var[8] else "-"
                cosmic_id = "Yes" if len(var) > 12 and var[12] else "-"
                omim_id = "Yes" if len(var) > 14 and var[14] else "-"

                # Shorten ClinVar significance for display
                if clinvar_sig != "-":
                    clinvar_sig = str(clinvar_sig).lower()
                    if "pathogenic" in clinvar_sig and "likely" not in clinvar_sig:
                        clinvar_sig = "Pathogenic"
                    elif "likely_pathogenic" in clinvar_sig:
                        clinvar_sig = "Likely Path"
                    elif "benign" in clinvar_sig and "likely" not in clinvar_sig:
                        clinvar_sig = "Benign"
                    elif "uncertain" in clinvar_sig:
                        clinvar_sig = "VUS"

                print(f"{score:<7.1f} {chrom:<6} {pos:<10} {ref}>{alt:<6} {gene or '-':<15} {clinvar_sig[:15]:<15} {omim_id:<8} {cosmic_id:<8}")

            return top_variants
        except Exception as e:
            print(f"Error exporting prioritized variants: {e}")
            import traceback
            traceback.print_exc()
            return []

    def close(self):
        """Close database connection"""
        if self.conn:
            self.conn.close()
            print(f"Closed connection to {self.db_path}")


def main():
    parser = argparse.ArgumentParser(description="Clinical variant analysis from OpenCravat results")
    parser.add_argument("db_path", help="Path to the OpenCravat SQLite database")
    parser.add_argument("--out", "-o", help="Output directory for reports", default="clinical_reports")
    parser.add_argument("--clinvar", "-c", action="store_true", help="Find high-quality ClinVar pathogenic variants")
    parser.add_argument("--all", "-A", action="store_true", help="Run all analyses")
    parser.add_argument("--limit", "-l", type=int, default=100, help="Limit results for each query")
    parser.add_argument("--output-csv", "-C", action="store_true", help="Output results to CSV files")
    parser.add_argument("--explore", "-e", action="store_true", help="Explore database schema")
    
    args = parser.parse_args()
    
    # Check if the database file exists
    if not os.path.isfile(args.db_path):
        print(f"Error: Database file '{args.db_path}' not found.")
        return 1
        
    # Create output directory if it doesn't exist
    if args.output_csv and not os.path.exists(args.out):
        os.makedirs(args.out)
        
    # Create analyzer instance
    analyzer = ClinicalVariantAnalyzer(args.db_path)
    
    # If explore mode, show schema for main tables
    if args.explore:
        main_tables = ["variant", "variant_annotator", "gene", "gene_annotator"]
        for table in main_tables:
            if table in analyzer.tables:
                analyzer.get_table_schema(table)
            else:
                print(f"Table '{table}' not found in database.")
        return 0
        
    # Run the analyses
    if args.all or args.clinvar:
        results = analyzer.find_clinvar_pathogenic(limit=args.limit)
        if args.output_csv and results:
            analyzer._write_to_csv(results, os.path.join(args.out, "clinvar_pathogenic.csv"))
            
    return 0




    def _get_tables(self):
        """Get all tables in the database"""
        cursor = self.conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = [row[0] for row in cursor.fetchall()]
        cursor.close()
        return tables

    def _count_variants(self):
        """Count total variants in the database"""
        cursor = self.conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM variant")
        count = cursor.fetchone()[0]
        cursor.close()
        return count

    def get_table_schema(self, table_name):
        """Get schema for a specific table"""
        cursor = self.conn.cursor()
        cursor.execute(f"PRAGMA table_info({table_name})")
        schema = cursor.fetchall()
        cursor.close()

        print(f"\n==== Schema for {table_name} ====")
        for col in schema:
            print(f"{col[0]}: {col[1]} ({col[2]})")

        return schema

    def find_clinvar_pathogenic(self, limit=100):
        """Find ClinVar pathogenic or likely pathogenic variants with high-quality evidence"""
        try:
            # Check if variant_annotator table exists
            if "variant_annotator" not in self.tables:
                print("variant_annotator table not found.")
                return []
                
            # First, get the schema of variant_annotator to see which columns exist
            cursor = self.conn.cursor()
            cursor.execute("PRAGMA table_info(variant_annotator)")
            cols = [col[1] for col in cursor.fetchall()]
            
            # Check for ClinVar columns
            clinvar_cols = [col for col in cols if "clinvar" in col.lower()]
            if not clinvar_cols:
                print("No ClinVar columns found in variant_annotator table.")
                print(f"Available columns: {', '.join(cols[:20])}...")
                return []
                
            # Find columns for clinical significance and review status
            sig_col = next((col for col in clinvar_cols if "sig" in col.lower()), None)
            review_col = next((col for col in clinvar_cols if "review" in col.lower()), None)
            star_col = next((col for col in clinvar_cols if "star" in col.lower()), None)
            disease_col = next((col for col in clinvar_cols if "disease" in col.lower() or "phenotype" in col.lower()), None)
            
            print(f"Using ClinVar columns: {sig_col}, {review_col}, {star_col}, {disease_col}")
            
            # Build the query based on available columns
            base_query = """
            SELECT v.base_uid, v.chrom, v.pos, v.ref_base, v.alt_base,
                   v.hugo, v.so, """
                   
            # Add CADD if available
            cadd_col = next((col for col in cols if "cadd" in col.lower() and "phred" in col.lower()), None)
            if cadd_col:
                base_query += f"va.{cadd_col}, "
            else:
                base_query += "NULL as cadd_phred, "
                
            # Add ClinVar columns
            if sig_col:
                base_query += f"va.{sig_col}, "
            else:
                base_query += "NULL as clinvar_sig, "
                
            if review_col:
                base_query += f"va.{review_col}, "
            else:
                base_query += "NULL as review_status, "
                
            if star_col:
                base_query += f"va.{star_col}, "
            else:
                base_query += "NULL as star_rating, "
                
            if disease_col:
                base_query += f"va.{disease_col}"
            else:
                base_query += "NULL as disease"
                
            query = base_query + """
            FROM variant v
            JOIN variant_annotator va ON v.base_uid = va.base_uid
            WHERE 1=1
            """
            
            # Add filter for pathogenic variants if sig_col exists
            if sig_col:
                query += f"""
                AND (va.{sig_col} LIKE '%pathogenic%'
                   OR va.{sig_col} LIKE '%Pathogenic%'
                   OR va.{sig_col} LIKE '%likely_pathogenic%'
                   OR va.{sig_col} LIKE '%Likely_pathogenic%')
                """
                
            # Add filter for high-quality evidence if review_col or star_col exists
            if review_col:
                query += f"""
                AND (va.{review_col} LIKE '%expert%' 
                   OR va.{review_col} LIKE '%panel%'
                   OR va.{review_col} LIKE '%practice%guideline%')
                """
                
            if star_col:
                query += f"""
                AND (va.{star_col} >= 3 OR va.{star_col} = '3' OR va.{star_col} = '4')
                """
                
            # Add ordering and limit
            if cadd_col:
                query += f"ORDER BY va.{cadd_col} DESC "
            query += "LIMIT ?"

            cursor = self.conn.cursor()
            cursor.execute(query, (limit,))
            results = cursor.fetchall()
            cursor.close()

            print(f"\n==== High-Quality ClinVar Pathogenic Variants ({len(results)}) ====\n")

            if results:
                # Print header
                header_parts = ["Chrom", "Position", "Ref>Alt", "Gene", "Effect"]
                if cadd_col:
                    header_parts.append("CADD")
                if sig_col:
                    header_parts.append("ClinVar Sig")
                if review_col:
                    header_parts.append("Review Status")
                if star_col:
                    header_parts.append("Stars")
                if disease_col:
                    header_parts.append("Disease")
                    
                # Format header
                header = ""
                for part in header_parts:
                    width = 30 if part == "Disease" else 20 if part == "Review Status" else 15 if part == "ClinVar Sig" else 6 if part == "CADD" or part == "Stars" else 10 if part == "Position" else 8 if part == "Ref>Alt" else 15
                    header += f"{part:<{width}} "
                print(header)
                print("-" * (len(header) + 10))

                # Print each result
                for row in results:
                    formatted_row = ""
                    # Basic variant info (always present)
                    uid, chrom, pos, ref, alt, gene, effect = row[0:7]
                    formatted_row += f"{chrom:<6} {pos:<10} {ref}>{alt:<6} {gene or 'Unknown':<15} {effect or 'Unknown':<15} "
                    
                    # Additional annotations (may or may not be present)
                    col_index = 7
                    if cadd_col:
                        cadd = row[col_index]
                        formatted_row += f"{cadd or 'N/A':<6} "
                        col_index += 1
                    if sig_col:
                        sig = row[col_index]
                        formatted_row += f"{str(sig)[:15] if sig else 'Unknown':<15} "
                        col_index += 1
                    if review_col:
                        review = row[col_index]
                        formatted_row += f"{str(review)[:20] if review else 'Unknown':<20} "
                        col_index += 1
                    if star_col:
                        stars = row[col_index]
                        formatted_row += f"{stars or 'N/A':<6} "
                        col_index += 1
                    if disease_col:
                        disease = row[col_index]
                        formatted_row += f"{str(disease)[:30] if disease else 'Unknown':<30}"
                        
                    print(formatted_row)

            return results
        except Exception as e:
            print(f"Error finding ClinVar pathogenic variants: {e}")
            import traceback
            traceback.print_exc()
            return []
            
    def _write_to_csv(self, results, filename):
        """Write query results to a CSV file"""
        try:
            with open(filename, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                # Write header - use generic column names as the structure may vary
                header = ["UID", "Chrom", "Position", "Ref", "Alt", "Gene", "Effect"]
                # Add extra columns based on result length
                for i in range(7, len(results[0]) if results else 7):
                    header.append(f"Column{i+1}")
                writer.writerow(header)
                
                # Write data
                for row in results:
                    writer.writerow(row)
                print(f"Results written to {filename}")
        except Exception as e:
            print(f"Error writing to CSV: {e}")

    def close(self):
        """Close database connection"""
        if self.conn:
            self.conn.close()


def main():
    parser = argparse.ArgumentParser(description="Clinical variant analysis from OpenCravat results")
    parser.add_argument("db_path", help="Path to the OpenCravat SQLite database")
    parser.add_argument("--out", "-o", help="Output directory for reports", default="clinical_reports")
    parser.add_argument("--clinvar", "-c", action="store_true", help="Find high-quality ClinVar pathogenic variants")
    parser.add_argument("--all", "-A", action="store_true", help="Run all analyses")
    parser.add_argument("--limit", "-l", type=int, default=100, help="Limit results for each query")
    parser.add_argument("--output-csv", "-C", action="store_true", help="Output results to CSV files")
    parser.add_argument("--explore", "-e", action="store_true", help="Explore database schema")
    
    args = parser.parse_args()
    
    # Check if the database file exists
    if not os.path.isfile(args.db_path):
        print(f"Error: Database file '{args.db_path}' not found.")
        return 1
        
    # Create output directory if it doesn't exist
    if args.output_csv and not os.path.exists(args.out):
        os.makedirs(args.out)
        
    # Create analyzer instance
    analyzer = ClinicalVariantAnalyzer(args.db_path)
    
    # If explore mode, show schema for main tables
    if args.explore:
        main_tables = ["variant", "variant_annotator", "gene", "gene_annotator"]
        for table in main_tables:
            if table in analyzer.tables:
                analyzer.get_table_schema(table)
            else:
                print(f"Table '{table}' not found in database.")
        return 0
        
    # Run the analyses
    if args.all or args.clinvar:
        results = analyzer.find_clinvar_pathogenic(limit=args.limit)
        if args.output_csv and results:
            analyzer._write_to_csv(results, os.path.join(args.out, "clinvar_pathogenic.csv"))
            
    return 0


if __name__ == "__main__":
    main()