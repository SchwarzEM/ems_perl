#!/usr/bin/env python

import argparse
import pandas as pd
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Filter bedtools -wa -wb overlap table to keep only same-reading-frame CDS overlaps."
    )
    parser.add_argument(
        "input_file",
        help="Input file: tab-delimited overlap table from bedtools (no header)"
    )
    parser.add_argument(
        "output_file",
        help="Output file name for same-frame overlaps"
    )
    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file

    # ====================== Column names ======================
    col_names = [
        'chrA', 'startA', 'endA', 'nameA', 'scoreA', 'strandA',
        'sourceA', 'typeA', 'phaseA', 'attrA',
        'chrB', 'startB', 'endB', 'nameB', 'scoreB', 'strandB',
        'sourceB', 'typeB', 'phaseB', 'attrB'
    ]

    def get_effective(start, end, strand, phase):
        """Return genomic position of codon start boundary."""
        try:
            phase = int(phase)
            start = int(start)
            end = int(end)
        except (ValueError, TypeError):
            return None
        if strand == '+':
            return start + phase
        else:
            return end - phase

    # ====================== Read and process ======================
    df = pd.read_csv(
        input_file,
        sep='\t',
        header=None,
        names=col_names,
        dtype=str,
        engine='python'
    )

    df['eff_A'] = df.apply(
        lambda row: get_effective(row['startA'], row['endA'], row['strandA'], row['phaseA']),
        axis=1
    )
    df['eff_B'] = df.apply(
        lambda row: get_effective(row['startB'], row['endB'], row['strandB'], row['phaseB']),
        axis=1
    )

    df['same_frame'] = (df['eff_A'] % 3 == df['eff_B'] % 3) & df['eff_A'].notna()
    df['strand_match'] = df['strandA'] == df['strandB']

    same_frame_df = df[df['same_frame'] & df['strand_match']]

    # ====================== Write output ======================
    same_frame_df[col_names].to_csv(
        output_file,
        sep='\t',
        header=False,
        index=False
    )

if __name__ == "__main__":
    main()
