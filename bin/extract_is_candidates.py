import argparse
import sys
import pandas as pd
from Bio import SeqIO

def get_reference_is6110_regions(record):
    """
    Extracts coordinates of existing IS6110 elements from the H37Rv GenBank record.
    Returns a list of (start, end) tuples.
    """
    is_regions = []
    for feature in record.features:
        if 'mobile_element_type' in feature.qualifiers:
            mobile_type = feature.qualifiers['mobile_element_type'][0]
            if 'IS6110' in mobile_type:
                # GenBank features use 0-based indexing for start, 1-based for end in Python Bio
                # But we just need the range to check overlaps
                start = int(feature.location.start)
                end = int(feature.location.end)
                is_regions.append((start, end))
    return is_regions

def is_overlap(pos, regions, padding=50):
    """
    Checks if a position is within any of the regions (with padding).
    """
    for start, end in regions:
        if (start - padding) <= pos <= (end + padding):
            return True
    return False

def parse_vcf_for_is6110(vcf_path, is_regions, decoy_name="IS6110", min_len=1300, max_len=1450):
    """
    Parses a VCF file for:
    1. INS variants with length ~1350bp.
    2. TRA/BND variants where one end is in a known IS6110 region.
    3. TRA/BND variants where one end is on the Decoy IS6110 contig.
    
    Returns a list of dictionaries with hit info.
    """
    hits = []
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            
            chrom = parts[0]
            pos = int(parts[1])
            info = parts[7]
            alt = parts[4]
            
            # Parse FORMAT and sample data for genotype read support
            dv, rv, dr, rr = 0, 0, 0, 0
            if len(parts) >= 10:
                format_keys = parts[8].split(':')
                sample_vals = parts[9].split(':')
                format_dict = dict(zip(format_keys, sample_vals))
                
                # In some VCFs missing values might not cast well, use 0 as fallback
                try: dv = int(format_dict.get('DV', 0))
                except ValueError: dv = 0
                try: rv = int(format_dict.get('RV', 0))
                except ValueError: rv = 0
                try: dr = int(format_dict.get('DR', 0))
                except ValueError: dr = 0
                try: rr = int(format_dict.get('RR', 0))
                except ValueError: rr = 0
                
            # Calculate VAF
            variant_support = dv + rv
            reference_support = dr + rr
            total_support = variant_support + reference_support
            vaf = 0.0
            if total_support > 0:
                vaf = variant_support / total_support
                
            # Determine Heteroresistance
            is_heteroresistant = 0.05 <= vaf <= 0.9
            
            # Parse INFO field
            info_dict = {}
            for item in info.split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    info_dict[k] = v
                else:
                    info_dict[item] = True
            
            sv_type = info_dict.get('SVTYPE', '')
            pe_support = int(info_dict.get('PE', 0))
            sr_support = int(info_dict.get('SR', 0))
            is_precise = 'PRECISE' in info_dict
            
            # Determine orientation based on connection type (CT)
            ct = info_dict.get('CT', '')
            if ct in ['3to3', '5to5']:
                orientation = 'Reverse'
            elif ct in ['3to5', '5to3']:
                orientation = 'Forward'
            else:
                orientation = 'Unknown'
            
            # --- Strategy 1: Direct Insertion (INS) ---
            if sv_type == 'INS':
                sv_len = 0
                if 'INSLEN' in info_dict:
                    sv_len = int(info_dict['INSLEN'])
                elif 'SVLEN' in info_dict:
                    sv_len = abs(int(info_dict['SVLEN']))
                
                if min_len <= sv_len <= max_len:
                    hits.append({
                        'chrom': chrom,
                        'pos': pos,
                        'sv_type': 'INS',
                        'length': sv_len,
                        'pe': pe_support,
                        'sr': sr_support,
                        'dv': dv,
                        'rv': rv,
                        'dr': dr,
                        'rr': rr,
                        'vaf': vaf,
                        'is_heteroresistant': is_heteroresistant,
                        'is_precise': is_precise,
                        'orientation': orientation,
                        'info': info,
                        'detection_method': 'Direct Size Match'
                    })

            # --- Strategy 2 & 3: Structural Variants (TRA/DEL/DUP/INV) ---
            elif sv_type in ['TRA', 'BND', 'DEL', 'DUP', 'INV']:
                # Get the other end coordinate
                chrom2 = info_dict.get('CHR2', chrom)
                # For TRA/BND, Delly usually provides POS2. For others, END.
                if 'POS2' in info_dict:
                    pos2 = int(info_dict['POS2'])
                else:
                    pos2 = int(info_dict.get('END', pos))
                
                # --- Strategy 3: Decoy Linkage (State-of-the-Art) ---
                # Check if one of the chromosomes is the decoy
                is_decoy_event = False
                insertion_chrom = None
                insertion_pos = None
                
                if chrom == decoy_name and chrom2 != decoy_name:
                    # Decoy is REF, Genome is ALT
                    is_decoy_event = True
                    insertion_chrom = chrom2
                    insertion_pos = pos2
                elif chrom2 == decoy_name and chrom != decoy_name:
                    # Genome is REF, Decoy is ALT
                    is_decoy_event = True
                    insertion_chrom = chrom
                    insertion_pos = pos
                
                if is_decoy_event:
                    # It's a link to the decoy. Now check if it's a known copy or new.
                    # Note: Known copies were masked, so reads mapped to decoy instead. 
                    # The breakpoint on the genome should be near the masked region.
                    if is_overlap(insertion_pos, is_regions):
                        hits.append({
                            'chrom': insertion_chrom,
                            'pos': insertion_pos,
                            'sv_type': 'TRA_Decoy_RefCopy',
                            'length': 0,
                            'pe': pe_support,
                            'sr': sr_support,
                            'dv': dv,
                            'rv': rv,
                            'dr': dr,
                            'rr': rr,
                            'vaf': vaf,
                            'is_heteroresistant': is_heteroresistant,
                            'is_precise': is_precise,
                            'orientation': orientation,
                            'info': info,
                            'detection_method': 'Decoy (Reference Copy)'
                        })
                    else:
                        hits.append({
                            'chrom': insertion_chrom,
                            'pos': insertion_pos,
                            'sv_type': 'TRA_Decoy_NewInsertion',
                            'length': 0,
                            'pe': pe_support,
                            'sr': sr_support,
                            'dv': dv,
                            'rv': rv,
                            'dr': dr,
                            'rr': rr,
                            'vaf': vaf,
                            'is_heteroresistant': is_heteroresistant,
                            'is_precise': is_precise,
                            'orientation': orientation,
                            'info': info,
                            'detection_method': 'Decoy (New Insertion)'
                        })
                    continue # Skip Strategy 2 check if Strategy 3 matches

                # --- Strategy 2: Link to Reference IS6110 (Legacy/Unmasked) ---
                # This logic is for when using the unmasked reference.
                if chrom2 == chrom:
                    pos1_is_is6110 = is_overlap(pos, is_regions)
                    pos2_is_is6110 = is_overlap(pos2, is_regions)
                    
                    if pos1_is_is6110 and not pos2_is_is6110:
                        hits.append({
                            'chrom': chrom,
                            'pos': pos2,
                            'sv_type': f'{sv_type}_IS_Link',
                            'length': 0, 
                            'pe': pe_support,
                            'sr': sr_support,
                            'dv': dv,
                            'rv': rv,
                            'dr': dr,
                            'rr': rr,
                            'vaf': vaf,
                            'is_heteroresistant': is_heteroresistant,
                            'is_precise': is_precise,
                            'orientation': orientation,
                            'info': info,
                            'detection_method': 'Linked to Reference IS6110'
                        })
                    elif pos2_is_is6110 and not pos1_is_is6110:
                        hits.append({
                            'chrom': chrom,
                            'pos': pos,
                            'sv_type': f'{sv_type}_IS_Link',
                            'length': 0, 
                            'pe': pe_support,
                            'sr': sr_support,
                            'dv': dv,
                            'rv': rv,
                            'dr': dr,
                            'rr': rr,
                            'vaf': vaf,
                            'is_heteroresistant': is_heteroresistant,
                            'is_precise': is_precise,
                            'orientation': orientation,
                            'info': info,
                            'detection_method': 'Linked to Reference IS6110'
                        })

    return hits

def apply_diagnostic_filters(hits, max_dist=10, min_pe_isolated=10, min_sr_isolated=5, min_total_clustered=15):
    """
    Applies the "Diagnostic-Grade" filtering logic.
    Rule 1: Isolated variants must have PE >= min_pe_isolated AND SR >= min_sr_isolated.
    Rule 2: Clustered variants (within max_dist bp) pass if the POOLED cluster evidence (Sum of PE + SR) >= min_total_clustered.
    """
    if not hits:
        return []

    # Sort hits by chromosome and position for efficient cluster identification
    hits.sort(key=lambda x: (x['chrom'], x['pos']))
    
    passed_indices = set()
    i = 0
    while i < len(hits):
        # Identify a cluster (Greedy approach similar to compilation)
        cluster_indices = [i]
        j = i + 1
        while j < len(hits) and hits[j]['chrom'] == hits[i]['chrom'] and (hits[j]['pos'] - hits[cluster_indices[-1]]['pos']) <= max_dist:
            cluster_indices.append(j)
            j = j + 1
        
        # Now evaluate the cluster
        if len(cluster_indices) == 1:
            # --- Rule 1: Isolated Variant Gate ---
            lead = hits[i]
            if lead['pe'] >= min_pe_isolated and lead['sr'] >= min_sr_isolated:
                passed_indices.add(i)
        else:
            # --- Rule 2: Clustered Variant Gate (Pooled Evidence) ---
            cluster_total_support = sum([(hits[idx]['pe'] + hits[idx]['sr']) for idx in cluster_indices])
            
            if cluster_total_support >= min_total_clustered:
                # All members of the cluster pass if the pooled evidence is sufficient
                for idx in cluster_indices:
                    passed_indices.add(idx)
        
        i = j # Move to the next unclustered signature
                
    # Return passed hits in sorted order
    passed_hits = [hits[idx] for idx in sorted(list(passed_indices))]
    return passed_hits

def annotate_hits(hits, record):
    """
    Annotates hits with gene information from a GenBank record.
    """
    annotated_hits = []
    
    for hit in hits:
        chrom = hit['chrom']
        pos = hit['pos']
        
        # Skip annotation if the hit is on the decoy itself (shouldn't happen with current logic but safe to check)
        if chrom == "IS6110": 
            hit['overlapping_genes'] = "Decoy"
            hit['nearest_gene'] = "Decoy"
            hit['dist_to_nearest'] = 0
            annotated_hits.append(hit)
            continue

        overlapping_features = []
        nearest_gene = None
        min_dist = float('inf')
        
        for feature in record.features:
            if feature.type not in ['CDS', 'gene', 'rRNA', 'tRNA']:
                continue
            
            start = int(feature.location.start)
            end = int(feature.location.end)
            locus_tag = feature.qualifiers.get('locus_tag', [''])[0]
            gene_name = feature.qualifiers.get('gene', [''])[0]
            product = feature.qualifiers.get('product', [''])[0]
            
            # Check overlap
            if start <= pos <= end:
                overlapping_features.append(f"{locus_tag} ({gene_name}): {product}")
            
            # Distance
            dist = min(abs(start - pos), abs(end - pos))
            if dist < min_dist:
                min_dist = dist
                nearest_gene = f"{locus_tag} ({gene_name})"

        hit['overlapping_genes'] = "; ".join(overlapping_features) if overlapping_features else "Intergenic"
        hit['nearest_gene'] = nearest_gene
        hit['dist_to_nearest'] = min_dist
        annotated_hits.append(hit)
        
    return annotated_hits

def main():
    parser = argparse.ArgumentParser(description="Filter Delly VCF for IS6110 insertions.")
    parser.add_argument("--vcf", required=True, help="Path to Delly VCF file")
    parser.add_argument("--gbk", required=True, help="Path to Reference GenBank file")
    parser.add_argument("--out", required=True, help="Output CSV path")
    parser.add_argument("--decoy", default="IS6110", help="Name of the IS6110 decoy contig (if used)")
    
    # Diagnostic Filter Thresholds
    parser.add_argument("--max-dist", type=int, default=10, help="Max distance between partners (default: 10)")
    parser.add_argument("--min-pe-iso", type=int, default=10, help="Min PE for isolated hits (default: 10)")
    parser.add_argument("--min-sr-iso", type=int, default=5, help="Min SR for isolated hits (default: 5)")
    parser.add_argument("--min-total-clust", type=int, default=15, help="Min PE+SR for clustered hits (default: 15)")
    
    args = parser.parse_args()
    
    print(f"Loading reference features from {args.gbk}...")
    try:
        record = SeqIO.read(args.gbk, "genbank")
    except Exception as e:
        print(f"Error loading GenBank: {e}")
        sys.exit(1)
        
    is_regions = get_reference_is6110_regions(record)
    print(f"Loaded {len(is_regions)} existing IS6110 regions from reference.")
    
    print(f"Parsing {args.vcf}...")
    hits = parse_vcf_for_is6110(args.vcf, is_regions, decoy_name=args.decoy)
    print(f"Found {len(hits)} candidate IS6110 events.")
    
    print(f"Applying Diagnostic-Grade filters (max_dist={args.max_dist}, isolated: PE>={args.min_pe_iso}, SR>={args.min_sr_iso}, clustered: total>={args.min_total_clust})...")
    filtered_hits = apply_diagnostic_filters(
        hits, 
        max_dist=args.max_dist, 
        min_pe_isolated=args.min_pe_iso, 
        min_sr_isolated=args.min_sr_iso, 
        min_total_clustered=args.min_total_clust
    )
    print(f"Diagnostic filtering complete. {len(filtered_hits)} hits passed.")
    
    print(f"Annotating hits...")
    annotated = annotate_hits(filtered_hits, record)
    
    df = pd.DataFrame(annotated)
    if not df.empty:
        cols = ['chrom', 'pos', 'detection_method', 'sv_type', 'pe', 'sr', 'dv', 'rv', 'dr', 'rr', 'vaf', 'is_heteroresistant', 'is_precise', 'orientation', 'overlapping_genes', 'nearest_gene', 'dist_to_nearest', 'info']
        # Reorder if columns exist
        existing_cols = [c for c in cols if c in df.columns]
        df = df[existing_cols]
    
    df.to_csv(args.out, index=False)
    print(f"Saved results to {args.out}")

if __name__ == "__main__":
    main()