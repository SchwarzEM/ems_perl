#!/usr/bin/env python

import sys

# 1. Check for command-line argument immediately
if len(sys.argv) < 2:
    print("Error: Missing input file argument.", file=sys.stderr)
    print("Usage: add_bh.fdrs_to_pvals_09jul2026a.py <path_to_table.txt>", file=sys.stderr)
    sys.exit(1)

from statsmodels.stats.multitest import multipletests

input_file = sys.argv[1]
genes = []
p_values = []

# Core file reading logic
try:
    with open(input_file, "r") as f:
        # Read and sanitize the header string safely
        first_line = f.readline().strip()

        if first_line:
            print("Pangene\tp-value\tFDR_BH")

        # Process all remaining rows sequentially
        for line in f:
            clean_line = line.strip()
            if not clean_line:
                continue

            parts = clean_line.split("\t")
            if len(parts) < 2:
                continue  # Skip corrupt or uneven lines

            gene_id = parts[0]
            p_val_str = parts[1]

            try:
                p_val = float(p_val_str)
                genes.append(gene_id)
                p_values.append(p_val)
            except ValueError:
                print(f"Skipping row with invalid p-value format: '{clean_line}'", file=sys.stderr)

except FileNotFoundError:
    print(f"Error: The file '{input_file}' was not found.", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"Error reading file structure: {e}", file=sys.stderr)
    sys.exit(1)

# Multiple testing correction and sorting block
if p_values:
    try:
        # Safely compute the BH FDR calculations using original sequence order
        mt_results = multipletests(p_values, alpha=0.001, method="fdr_bh")
        q_values = mt_results[1]

        # Combine arrays into a structured list of lists
        combined_data = []
        for i in range(len(genes)):
            combined_data.append([genes[i], p_values[i], q_values[i]])

        # Sort the data dynamically based on the p-value column (index 1)
        # x[1] references the float p-value item in each sub-list
        sorted_data = sorted(combined_data, key=lambda x: x[1])

        # 3. Print unified sorted rows cleanly to STDOUT
        for row in sorted_data:
            print(row[0], row[1], row[2], sep="\t")

    except Exception as e:
        print(f"Error computing Benjamini-Hochberg FDR adjustments: {e}", file=sys.stderr)
        sys.exit(1)
else:
    print("Warning: No valid p-values were found to process.", file=sys.stderr)
