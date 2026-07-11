#!/usr/bin/env python3
import pandas as pd
import argparse
import os
from compile_results import collapse_dsb_signatures

def filter_region(input_csv, raw_out, formatted_out, chrom, start, end, region_name, max_dist=10):
    """
    Filter for variants in a specific genomic region.
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
    print(f"Filtering for {region_name} coordinates ({chrom}: {start}-{end})...")
    df['pos'] = pd.to_numeric(df['pos'], errors='coerce')
    region_df = df[(df['chrom'] == chrom) & (df['pos'] >= start) & (df['pos'] <= end)]
    
    print(f"Found {len(region_df)} variants in {region_name} region.")

    # Ensure output directories exist
    os.makedirs(os.path.dirname(raw_out), exist_ok=True)
    os.makedirs(os.path.dirname(formatted_out), exist_ok=True)

    # Save Raw Filtered Result
    region_df.to_csv(raw_out, index=False)
    print(f"Saved raw filtered results to {raw_out}")
    
    # Create and Save Formatted Final Report
    print(f"Creating formatted {region_name} report...")
    final_report = collapse_dsb_signatures(region_df, max_dist=max_dist)
    final_report.to_csv(formatted_out, index=False)
    print(f"Saved formatted {region_name} report to {formatted_out}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter results for specific genomic region.")
    parser.add_argument("--input", required=True, help="Input compiled CSV")
    parser.add_argument("--raw-out", required=True, help="Output path for raw filtered CSV")
    parser.add_argument("--formatted-out", required=True, help="Output path for formatted report")
    parser.add_argument("--chrom", default="NC_000962.3", help="Target chromosome/sequence")
    parser.add_argument("--start", type=int, default=775586, help="Start coordinate")
    parser.add_argument("--end", type=int, default=779487, help="End coordinate")
    parser.add_argument("--region-name", default="mmpL5-Rv0678", help="Name of the region/genes")
    parser.add_argument("--max-dist", type=int, default=10, help="Max distance to collapse junctions")
    
    args = parser.parse_args()
    raw_out_path = os.path.abspath(args.raw_out)
    formatted_out_path = os.path.abspath(args.formatted_out)
    
    filter_region(args.input, raw_out_path, formatted_out_path, args.chrom, args.start, args.end, args.region_name, max_dist=args.max_dist)
