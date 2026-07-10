import pandas as pd
import os
import glob
import argparse
import re

def clean_gene_name(name):
    """Removes empty parentheses '()' from gene names."""
    if pd.isna(name):
        return name
    return re.sub(r'\s*\(\)', '', str(name)).strip()

def collapse_dsb_signatures(df, max_dist=10):
    """
    Collapses nearby structural variant junctions (DSB signatures) into single rows.
    Implements a greedy clustering approach to group all signatures within max_dist.
    """
    if df.empty:
        return df

    # Ensure necessary columns exist for sorting and grouping
    sort_cols = ['SampleID', 'chrom', 'pos']
    df = df.sort_values(sort_cols)
    
    # Rename columns to start with
    df = df.rename(columns={'SampleID': 'Sample', 'chrom': 'Chrom'})
    
    collapsed_rows = []
    
    # Group by sample and chromosome to find neighbors
    for (sample, chrom), group in df.groupby(['Sample', 'Chrom']):
        group = group.reset_index(drop=True)
        i = 0
        while i < len(group):
            # Start a new cluster
            cluster = [group.iloc[i]]
            j = i + 1
            # Greedily add neighbors within max_dist of the LAST added member
            while j < len(group) and (group.iloc[j]['pos'] - cluster[-1]['pos']) <= max_dist:
                cluster.append(group.iloc[j])
                j += 1
            
            # --- New Clustering Logic ---
            
            # 1. Evidence Aggregation (Sum all evidence in cluster)
            total_pe_cluster = sum([c['pe'] for c in cluster])
            total_sr_cluster = sum([c['sr'] for c in cluster])
            
            total_dv_cluster = sum([c.get('dv', 0) for c in cluster if pd.notna(c.get('dv'))])
            total_rv_cluster = sum([c.get('rv', 0) for c in cluster if pd.notna(c.get('rv'))])
            total_dr_cluster = sum([c.get('dr', 0) for c in cluster if pd.notna(c.get('dr'))])
            total_rr_cluster = sum([c.get('rr', 0) for c in cluster if pd.notna(c.get('rr'))])
            
            cluster_vaf_variant = total_dv_cluster + total_rv_cluster
            cluster_vaf_ref = total_dr_cluster + total_rr_cluster
            cluster_vaf_total = cluster_vaf_variant + cluster_vaf_ref
            cluster_vaf = (cluster_vaf_variant / cluster_vaf_total) if cluster_vaf_total > 0 else 0.0
            cluster_hetero = 0.05 <= cluster_vaf <= 0.35
            
            # 2. Peak Selection (Find top 2 coordinates by individual evidence)
            # Group cluster by position to handle multiple entries at same pos if any
            pos_evidence = []
            for pos in sorted(list(set([c['pos'] for c in cluster]))):
                pos_entries = [c for c in cluster if c['pos'] == pos]
                sum_pe = sum([e['pe'] for e in pos_entries])
                sum_sr = sum([e['sr'] for e in pos_entries])
                is_prec = any([e['is_precise'] for e in pos_entries])
                pos_evidence.append({
                    'pos': pos,
                    'total': sum_pe + sum_sr,
                    'is_precise': is_prec
                })
            
            # Sort by precision (True first), then total evidence descending, then by position ascending
            pos_evidence.sort(key=lambda x: (not x['is_precise'], -x['total'], x['pos']))
            
            top_peaks = pos_evidence[:2]
            # Sort peaks by coordinate for A/B reporting consistency (A < B)
            top_peaks.sort(key=lambda x: x['pos'])
            
            pos_a = top_peaks[0]['pos']
            prec_a = top_peaks[0]['is_precise']
            
            pos_b = None
            prec_b = None
            if len(top_peaks) > 1:
                pos_b = top_peaks[1]['pos']
                prec_b = top_peaks[1]['is_precise']

            # 3. Metadata Aggregation
            orientations = sorted(list(set([str(c['orientation']) for c in cluster])))
            orientation_str = "/".join(orientations) if len(orientations) > 1 else orientations[0]
            
            methods = []
            for c in cluster:
                m = str(c['detection_method'])
                for sub_m in m.split('; '):
                    if sub_m not in methods:
                        methods.append(sub_m)
            method_str = "; ".join(methods)
            
            overlapping = "; ".join(sorted(list(set([clean_gene_name(c['overlapping_genes']) for c in cluster if not pd.isna(c['overlapping_genes'])]))))
            # For nearest gene, we pick the one closest to our primary Peak A
            nearest_entries = sorted(cluster, key=lambda x: abs(x['pos'] - pos_a))
            nearest = nearest_entries[0]['nearest_gene']
            dist_to_nearest = min([c['dist_to_nearest'] for c in cluster])

            row = {
                'Sample': sample,
                'Chrom': chrom,
                'DSB Breakpoint A': pos_a,
                'DSB Breakpoint B': pos_b,
                'Orientation': orientation_str,
                'Overlapping Genes': overlapping,
                'Nearest Gene': clean_gene_name(nearest),
                'Dist. to Nearest': dist_to_nearest,
                'PE A': total_pe_cluster,
                'PE B': total_pe_cluster if pos_b is not None else None,
                'SR A': total_sr_cluster,
                'SR B': total_sr_cluster if pos_b is not None else None,
                'DV': total_dv_cluster,
                'RV': total_rv_cluster,
                'DR': total_dr_cluster,
                'RR': total_rr_cluster,
                'VAF': cluster_vaf,
                'Is Heteroresistant': cluster_hetero,
                'Precision A': prec_a,
                'Precision B': prec_b,
                'Detection Method': method_str
            }
            
            collapsed_rows.append(row)
            i = j # Move to the next unclustered signature
                
    result_df = pd.DataFrame(collapsed_rows)
    
    # Final Column Order
    target_order = [
        'Sample', 'Chrom', 'DSB Breakpoint A', 'DSB Breakpoint B', 'Orientation',
        'Overlapping Genes', 'Nearest Gene', 'Dist. to Nearest', 
        'PE A', 'PE B', 'SR A', 'SR B', 
        'DV', 'RV', 'DR', 'RR', 'VAF', 'Is Heteroresistant',
        'Precision A', 'Precision B', 'Detection Method'
    ]
    
    # Ensure all columns exist (even if empty) for robustness
    for col in target_order:
        if col not in result_df.columns:
            result_df[col] = None
            
    return result_df[target_order]

def compile_results(results_dir, raw_out, formatted_out, max_dist=10):
    """
    Finds all *_is6110_candidates.csv files, merges them,
    and produces both raw compiled and formatted final reports.
    """
    search_pattern = os.path.join(results_dir, "*", "*_is6110_candidates.csv")
    candidate_files = glob.glob(search_pattern)
    
    if not candidate_files:
        print(f"No candidate files found in {results_dir}")
        return

    print(f"Found {len(candidate_files)} candidate files. Compiling...")
    
    all_dfs = []
    for f in candidate_files:
        filename = os.path.basename(f)
        sample_id = filename.replace("_is6110_candidates.csv", "")
        
        try:
            df = pd.read_csv(f)
            if not df.empty:
                df.insert(0, 'SampleID', sample_id)
                all_dfs.append(df)
        except Exception as e:
            print(f"Error reading {f}: {e}")

    if not all_dfs:
        print("No data found in any of the candidate files.")
        return

    merged_df = pd.concat(all_dfs, ignore_index=True)
    
    # Ensure output directories exist
    os.makedirs(os.path.dirname(raw_out), exist_ok=True)
    os.makedirs(os.path.dirname(formatted_out), exist_ok=True)
    
    # Save Raw Merged Result
    merged_df.to_csv(raw_out, index=False)
    print(f"Saved raw compiled results to {raw_out}")
    
    # Create and Save Formatted Final Report
    print(f"Creating formatted final report...")
    final_report = collapse_dsb_signatures(merged_df, max_dist=max_dist)
    final_report.to_csv(formatted_out, index=False)
    print(f"Saved formatted final report to {formatted_out}")
    
    # Save compressed raw results
    compressed_raw = raw_out + ".gz"
    merged_df.to_csv(compressed_raw, index=False, compression='gzip')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compile IS6110 candidate CSVs.")
    parser.add_argument("--results", default="results", help="Directory containing sample results")
    parser.add_argument("--raw-out", required=True, help="Output path for raw compiled CSV")
    parser.add_argument("--formatted-out", required=True, help="Output path for formatted final report")
    parser.add_argument("--max-dist", type=int, default=10, help="Max distance to collapse junctions")
    
    args = parser.parse_args()
    
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    results_path = os.path.join(base_dir, args.results)
    
    compile_results(results_path, args.raw_out, args.formatted_out, max_dist=args.max_dist)
