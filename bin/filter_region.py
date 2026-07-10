import pandas as pd
import argparse
import os
from compile_results import collapse_dsb_signatures

def filter_region(input_csv, raw_out, formatted_out, max_dist=10):
    """
    Filter for variants in a specific genomic region (currently mmpL5-Rv0678).
    Produces both raw filtered and formatted final reports.
    """
    if not os.path.exists(input_csv):
        print(f"Error: Input file {input_csv} not found.")
        return

    print(f"Loading {input_csv}...")
    df = pd.read_csv(input_csv)
    
    # --- Cleaning ---
    initial_count = len(df)
    df.drop_duplicates(inplace=True)
    
    # Dropping simulation tests
    df = df[~df['SampleID'].str.contains('simulated_test', na=False, case=False)]
    print(f"Loaded {initial_count} records, cleaned to {len(df)} records.")

    # --- Filtering Logic ---
    chrom_target = 'NC_000962.3'
    # Window: mmpL5 (Rv0676c) to Rv0678
    start_pos = 775586
    end_pos = 779487
    
    print(f"Filtering for mmpL5-Rv0678 coordinates ({chrom_target}: {start_pos}-{end_pos})...")
    df['pos'] = pd.to_numeric(df['pos'], errors='coerce')
    region_df = df[(df['chrom'] == chrom_target) & (df['pos'] >= start_pos) & (df['pos'] <= end_pos)]
    
    print(f"Found {len(region_df)} variants in mmpL5-Rv0678 region.")

    # Ensure output directories exist
    os.makedirs(os.path.dirname(raw_out), exist_ok=True)
    os.makedirs(os.path.dirname(formatted_out), exist_ok=True)

    # Save Raw Filtered Result
    region_df.to_csv(raw_out, index=False)
    print(f"Saved raw filtered results to {raw_out}")
    
    # Create and Save Formatted Final Report
    print(f"Creating formatted mmpL5-Rv0678 report...")
    final_report = collapse_dsb_signatures(region_df, max_dist=max_dist)
    final_report.to_csv(formatted_out, index=False)
    print(f"Saved formatted mmpL5-Rv0678 report to {formatted_out}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter results for specific genomic region.")
    parser.add_argument("--input", required=True, help="Input compiled CSV")
    parser.add_argument("--raw-out", required=True, help="Output path for raw filtered CSV")
    parser.add_argument("--formatted-out", required=True, help="Output path for formatted report")
    parser.add_argument("--max-dist", type=int, default=10, help="Max distance to collapse junctions")
    
    args = parser.parse_args()
    
    filter_region(args.input, args.raw_out, args.formatted_out, max_dist=args.max_dist)
