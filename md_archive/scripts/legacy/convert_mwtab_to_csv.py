#!/usr/bin/env python3
"""
Utility script to fetch mwTab content from Metabolomics Workbench
and convert the tabular section to CSV format.

Usage:
    python scripts/convert_mwtab_to_csv.py --study-id ST004396
    python scripts/convert_mwtab_to_csv.py --study-id ST004396 --output-dir data/mwtab_export
"""

import argparse
import csv
import json
import sys
from pathlib import Path

# Add parent directory to path to import from scripts
sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.harvest_mw_studies import fetch_mw_mwtab


def parse_json_mwtab(mwtab_text: str):
    """
    Parse JSON-formatted mwTab to extract metabolite data.

    Handles mwTab JSON structure with MS_METABOLITE_DATA section.

    Args:
        mwtab_text: Raw mwTab text content (JSON format)

    Returns:
        List of lists representing rows (first row is header), or empty list
    """
    rows = []

    # Try to parse JSON (may have multiple JSON objects or extra content)
    text = mwtab_text.strip()
    first_brace = text.find("{")

    if first_brace == -1:
        return rows

    # Try to parse first complete JSON object
    for end_pos in range(len(text), first_brace, -1):
        try:
            json_str = text[first_brace:end_pos]
            data = json.loads(json_str)

            # Look for metabolite data sections
            metabolite_sections = [
                "MS_METABOLITE_DATA",
                "GC_METABOLITE_DATA",
                "LC_METABOLITE_DATA",
                "METABOLITE_DATA",
            ]

            for section_key in metabolite_sections:
                if section_key not in data:
                    continue

                section = data[section_key]
                if not isinstance(section, dict):
                    continue

                # Extract data array
                data_array = section.get("Data", [])
                if not isinstance(data_array, list) or len(data_array) == 0:
                    continue

                # First row is a dict - extract keys (Metabolite + sample IDs)
                if isinstance(data_array[0], dict):
                    # Get all unique keys across all rows (for header)
                    all_keys = set()
                    for row_dict in data_array:
                        if isinstance(row_dict, dict):
                            all_keys.update(row_dict.keys())

                    # Use 'Metabolite' as first column, then sample IDs
                    header = ["Metabolite"]
                    sample_keys = sorted([k for k in all_keys if k != "Metabolite"])
                    header.extend(sample_keys)
                    rows.append(header)

                    # Convert each dict to a row
                    for row_dict in data_array:
                        if not isinstance(row_dict, dict):
                            continue
                        row = [
                            str(row_dict.get("Metabolite", "")),
                            *[str(row_dict.get(key, "")) for key in sample_keys],
                        ]
                        rows.append(row)

                    if rows:
                        return rows

                # Alternative format: Data might be list of lists
                elif isinstance(data_array[0], list):
                    # First list is header
                    if len(data_array) > 0:
                        rows.append(data_array[0])  # header
                        rows.extend(data_array[1:])  # data rows
                        return rows

            # If we found valid JSON but no metabolite data, break
            break

        except (json.JSONDecodeError, ValueError, KeyError):
            continue

    return rows


def parse_mwtab_to_rows(mwtab_text: str):
    """
    Parse mwTab content to extract tabular data.

    Handles multiple formats:
    1. Tab-separated text (legacy format) - lines with tabs
    2. JSON format (newer format) - extracts metabolite data from JSON structure

    Args:
        mwtab_text: Raw mwTab text content

    Returns:
        List of lists representing rows (first row is header)
    """
    rows = []

    # First, try to find tab-separated sections (most common format)
    header_found = False
    for line in mwtab_text.splitlines():
        line = line.rstrip("\n\r")

        # Skip comment/metadata lines
        if line.startswith("#") or not line:
            continue

        # Only consider tab-separated lines
        if "\t" not in line:
            continue

        parts = line.split("\t")
        if not header_found:
            header_found = True
            rows.append(parts)
        else:
            rows.append(parts)

    # If we found tab-separated data, return it
    if rows:
        return rows

    # Otherwise, try to parse as JSON format
    rows = parse_json_mwtab(mwtab_text)
    if rows:
        return rows

    # If still no data, return empty (will be handled by caller)
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Download mwTab for a Metabolomics Workbench study and convert the tabular section to CSV."
    )
    parser.add_argument(
        "--study-id",
        required=True,
        help="Metabolomics Workbench study ID (e.g. ST004396).",
    )
    parser.add_argument(
        "--output-dir",
        default="data/mwtab",
        help="Directory to save the CSV file (default: data/mwtab).",
    )
    args = parser.parse_args()

    study_id = args.study_id.upper().strip()
    print(f"Fetching mwTab for study {study_id}...")

    try:
        mwtab_text = fetch_mw_mwtab(study_id)
    except Exception as e:
        print(f"Error fetching mwTab for {study_id}: {e}")
        return 1

    if not mwtab_text or len(mwtab_text.strip()) == 0:
        print(f"No mwTab content returned for study {study_id}")
        return 1

    print(f"Parsing mwTab content ({len(mwtab_text)} characters)...")
    rows = parse_mwtab_to_rows(mwtab_text)

    if not rows:
        print(f"⚠️  No tabular rows parsed from mwTab for study {study_id}")
        print(
            f"   This study's mwTab format may be JSON-only without tab-separated data sections."
        )
        print(
            f"   The mwTab content ({len(mwtab_text)} chars) is available but doesn't contain"
        )
        print(f"   tab-separated metabolite data tables in the expected format.")
        print(
            f"\n   Tip: Try a different study ID, or check the mwTab content manually:"
        )
        print(
            f"   curl 'https://www.metabolomicsworkbench.org/rest/study/study_id/{study_id}/mwtab'"
        )
        return 1

    if len(rows) == 1:
        print(f"Warning: Only header row found, no data rows for study {study_id}")
    else:
        print(f"Parsed {len(rows)} rows ({len(rows) - 1} data rows)")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{study_id}.csv"

    print(f"Writing CSV to {output_path}...")
    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        for row in rows:
            writer.writerow(row)

    print(f"✅ Saved CSV for {study_id} to {output_path}")
    print(f"   File size: {output_path.stat().st_size} bytes")
    print(f"   Rows: {len(rows)} ({len(rows[0])} columns)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
