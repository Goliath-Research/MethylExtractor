#!/usr/bin/env python3
"""
Simple script to convert MethylExtractor CSV output to Parquet format with Zstd compression.
This demonstrates the Parquet conversion functionality.
"""

import pandas as pd
import sys
import os

def csv_to_parquet(csv_file):
    """Convert CSV file to Parquet with Zstd compression."""
    if not os.path.exists(csv_file):
        print(f"Error: {csv_file} not found")
        return False
    
    parquet_file = csv_file.replace('.csv', '.parquet')
    
    try:
        # Read CSV
        df = pd.read_csv(csv_file)
        print(f"Read {len(df)} rows from {csv_file}")
        
        # Write Parquet with Zstd compression
        df.to_parquet(parquet_file, compression='zstd', index=False)
        
        # Get file sizes
        csv_size = os.path.getsize(csv_file)
        parquet_size = os.path.getsize(parquet_file)
        ratio = parquet_size / csv_size
        
        print(f"Converted to {parquet_file}")
        print(f"Compression ratio: {ratio:.1f}")
        return True
        
    except Exception as e:
        print(f"Error converting {csv_file}: {e}")
        return False

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 test_parquet.py <csv_file>")
        sys.exit(1)
    
    csv_file = sys.argv[1]
    success = csv_to_parquet(csv_file)
    sys.exit(0 if success else 1)
