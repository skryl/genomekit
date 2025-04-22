#!/usr/bin/env python3
import sqlite3
import sys

def check_variant_columns(db_path):
    """Check what columns are available in the variant table"""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Get column names from the variant table
    cursor.execute("PRAGMA table_info(variant)")
    columns = cursor.fetchall()
    
    print(f"Total columns in variant table: {len(columns)}")
    print("\nColumns that might contain genotype information:")
    
    # Print columns that might contain genotype information
    for col in columns:
        col_id, col_name, col_type, _, _, _ = col
        if "genotype" in col_name.lower() or "zygos" in col_name.lower() or "allele" in col_name.lower():
            print(f"- {col_name} ({col_type})")
    
    print("\nSample of base__ columns:")
    for col in columns:
        col_id, col_name, col_type, _, _, _ = col
        if col_name.startswith("base__"):
            print(f"- {col_name} ({col_type})")
    
    conn.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python check_variant_columns.py <db_path>")
        sys.exit(1)
    
    check_variant_columns(sys.argv[1])
