#!/usr/bin/env python3
import pandas as pd
import re
import sys

def extract_te_copy(metadata):
    """Extract the TE_copy from the metadata using the ID= field."""
    match = re.search(r'ID=([^;]+)', metadata)
    return match.group(1) if match else None

def extract_te_consensus(metadata):
    """Extract the TE_consensus from the metadata using the Target= field."""
    match = re.search(r'Target=([^;\s]+)', metadata)
    return match.group(1) if match else None

def main(input_file, output_file):
    # Read the file (assuming tab-delimited)
    try:
        df = pd.read_csv(input_file, sep="\t")
    except Exception as e:
        sys.exit(f"Error reading file: {e}")

    # Create new columns extracting the required fields from the metadata column
    df['TE_copy'] = df['metadata'].apply(extract_te_copy)
    df['TE_consensus'] = df['metadata'].apply(extract_te_consensus)
    
    # Optionally, you can rearrange columns if desired.
    # For example, to have TE_copy and TE_consensus first:
    cols = ['TE_copy', 'TE_consensus'] + [col for col in df.columns if col not in ['TE_copy', 'TE_consensus']]
    df = df[cols]

    # Write to a CSV file
    df.to_csv(output_file, index=False)
    print(f"CSV file saved as {output_file}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Generate a CSV file with TE_copy and TE_consensus extracted from a GFF file."
    )
    parser.add_argument('-i', '--input', type=str, required=True, help="Input GFF file")
    parser.add_argument('-o', '--output', type=str, required=True, help="Output CSV file")
    
    args = parser.parse_args()
    main(args.input, args.output)