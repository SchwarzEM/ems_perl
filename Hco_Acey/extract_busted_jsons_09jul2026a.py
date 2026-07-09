#!/usr/bin/env python

import json
import re
import sys
from pathlib import Path

if len(sys.argv) < 2:
    print("Error: Missing input file argument.", file=sys.stderr)
    print("Usage: extract_busted_jsons_09jul2026a.py <path_to_JSON_file_list.txt>", file=sys.stderr)
    sys.exit(1)

input_list_file = sys.argv[1]

# 1. Specify the header text before the loop starts
header_text = "Pangene\tp-value"

with open(input_list_file, "r") as list_f:
    for line in list_f:
        clean_path_str = line.strip()

        if not clean_path_str:
            continue

        json_path = Path(clean_path_str)

        try:
            with open(json_path) as f:
                data = json.load(f)

            # 2. Conditional check: if header_text has content, print and clear it
            if header_text:
                print(header_text)
                header_text = ""  # Reset to empty string so it never prints again

            # Extract the values
            match = re.search(r"pangene_\d+", json_path.name)
            pangene_name = match.group(0) if match else "None"
            p_value = data.get("test results", {}).get("p-value", "None")

            # Print data to STDOUT
            print(pangene_name, p_value, sep="\t")

        except (FileNotFoundError, json.JSONDecodeError) as e:
            print(f"Error processing {json_path}: {e}", file=sys.stderr)
