#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def get_is6110_coordinates(gbk_path):
    """Parses GenBank to find IS6110 coordinates."""
    is_regions = []
    try:
        # GBK might contain multiple records, usually H37Rv is one.
        records = list(SeqIO.parse(gbk_path, "genbank"))
    except Exception as e:
        print(f"Error reading GenBank: {e}")
        return []

    for record in records:
        for feature in record.features:
            if 'mobile_element_type' in feature.qualifiers:
                mobile_type = feature.qualifiers['mobile_element_type'][0]
                if 'IS6110' in mobile_type:
                    # 0-based start, 1-based end (BioPython location)
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    is_regions.append((start, end))
    return is_regions

def mask_sequence(seq_record, regions):
    """Masks the sequence with Ns at specified regions."""
    try:
        # For Biopython < 1.79
        mutable_seq = seq_record.seq.tomutable()
    except AttributeError:
        # For Biopython >= 1.79
        from Bio.Seq import MutableSeq
        mutable_seq = MutableSeq(str(seq_record.seq))
        
    for start, end in regions:
        # Check bounds
        s = max(0, start)
        e = min(len(mutable_seq), end)
        length = e - s
        if length > 0:
            mutable_seq[s:e] = "N" * length
            
    try:
        # For Biopython < 1.79
        seq_record.seq = mutable_seq.toseq()
    except AttributeError:
        # For Biopython >= 1.79
        from Bio.Seq import Seq
        seq_record.seq = Seq(str(mutable_seq))
        
    return seq_record

def main():
    parser = argparse.ArgumentParser(description="Create a masked reference genome with a decoy IS6110 contig.")
    parser.add_argument("--gbk", required=True, help="Reference GenBank file (for IS coordinates)")
    parser.add_argument("--fasta", required=True, help="Reference FASTA file (to mask)")
    parser.add_argument("--is6110", required=True, help="IS6110 sequence FASTA (decoy)")
    parser.add_argument("--out", required=True, help="Output FASTA file")
    
    args = parser.parse_args()
    
    # 1. Get Coordinates
    print(f"Reading coordinates from {args.gbk}...")
    coords = get_is6110_coordinates(args.gbk)
    print(f"Found {len(coords)} IS6110 regions to mask.")
    
    # 2. Read Reference
    print(f"Reading reference from {args.fasta}...")
    ref_records = list(SeqIO.parse(args.fasta, "fasta"))
    if not ref_records:
        print("No records found in reference FASTA.")
        exit(1)
        
    # 3. Mask Reference (Assuming first record is the chromosome)
    # If there are plasmids, we might need logic, but usually H37Rv ref is one seq.
    # We apply masking to the main chromosome (NC_000962.3)
    # The coordinates from GBK usually correspond to the first record if it's the main one.
    
    # For safety, let's assume coordinates apply to the first record found (H37Rv).
    print(f"Masking {len(coords)} regions in {ref_records[0].id}...")
    ref_records[0] = mask_sequence(ref_records[0], coords)
    
    # 4. Read Decoy
    print(f"Reading IS6110 decoy from {args.is6110}...")
    decoy_records = list(SeqIO.parse(args.is6110, "fasta"))
    if not decoy_records:
        print("No records found in IS6110 FASTA.")
        exit(1)
        
    # Rename decoy to something simple if needed, e.g. 'IS6110'
    decoy_record = decoy_records[0]
    decoy_record.id = "IS6110"
    decoy_record.description = "IS6110 Decoy Contig"
    
    # 5. Write Output
    output_records = ref_records + [decoy_record]
    print(f"Writing masked reference + decoy to {args.out}...")
    SeqIO.write(output_records, args.out, "fasta")
    print("Done.")

if __name__ == "__main__":
    main()
