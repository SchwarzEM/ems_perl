#!/usr/bin/env python

import json
import re
import sys
from pathlib import Path

if len(sys.argv) < 2:
    print("Error: Missing input file argument.", file=sys.stderr)
    print("Usage: extract_absrel_jsons_09jul2026a.py <path_to_JSON_file_list.txt>", file=sys.stderr)
    sys.exit(1)

input_list_file = sys.argv[1]

header_text = "Pangene\tBranch\tLRT\tUncorrected_P-value\tCorrected_P-value"

with open(input_list_file, "r") as list_f:
    for line in list_f:
        clean_path_str = line.strip()

        if not clean_path_str:
            continue

        json_path = Path(clean_path_str)

        try:
            with open(json_path) as f:
                data = json.load(f)

            if header_text:
                print(header_text)
                header_text = ""

            # Extract pangene name from filename
            match = re.search(r"pangene_\d+", json_path.name)
            pangene_name = match.group(0) if match else json_path.stem

            # === 1. Find branches marked as "test" (foreground) ===
            tested_section = data.get("tested", {})
            test_branches = set()

            for _, branches_dict in tested_section.items():
                if isinstance(branches_dict, dict):
                    for branch_name, status in branches_dict.items():
                        if str(status).lower() == "test":
                            test_branches.add(str(branch_name))

            if not test_branches:
                continue

            # === 2. Get branch data (handles both flat and nested "0" structures) ===
            ba = data.get("branch attributes", {})
            if "0" in ba and isinstance(ba["0"], dict):
                branch_attrs = ba["0"]
            else:
                branch_attrs = ba

            # === 3. Extract data only for test/foreground branches ===
            for branch_key, attrs in branch_attrs.items():
                if not isinstance(attrs, dict) or branch_key == "attributes":
                    continue

                raw_original = attrs.get("original name", branch_key)
                branch_id = str(branch_key)
                original_id = str(raw_original) if raw_original is not None else branch_id

                if branch_id not in test_branches and original_id not in test_branches:
                    continue

                lrt      = attrs.get("LRT", "None")
                uncorr_p = attrs.get("Uncorrected P-value", "None")
                corr_p   = attrs.get("Corrected P-value", "None")

                display_name = original_id if original_id != branch_id else branch_id

                print(pangene_name, display_name, lrt, uncorr_p, corr_p, sep="\t")

        except Exception as e:
            print(f"Error processing {json_path}: {e}", file=sys.stderr)
