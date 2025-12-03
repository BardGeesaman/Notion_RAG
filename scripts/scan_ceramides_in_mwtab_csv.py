#!/usr/bin/env python3
"""
Utility script to scan mwTab CSV files for ceramide-like metabolites.

Scans metabolite names for ceramide patterns and produces a summary
of which ceramide species appear and how frequently.

Usage:
    python scripts/scan_ceramides_in_mwtab_csv.py --study-id ST004396
    python scripts/scan_ceramides_in_mwtab_csv.py --csv-path data/mwtab/ST004396.csv
    python scripts/scan_ceramides_in_mwtab_csv.py --study-id ST004396 --output-report report.txt
"""

import argparse
import csv
from collections import Counter
from pathlib import Path


def is_ceramide_like(name: str) -> bool:
    """
    Simple ceramide pattern detector, case-insensitive.

    Examples that should match:
    - Cer(d18:1/16:0)
    - Ceramide C16:0
    - CER 24:1

    Args:
        name: Metabolite name to check

    Returns:
        True if name matches ceramide patterns, False otherwise
    """
    if not name:
        return False

    lower = name.lower()

    patterns = [
        "cer(d",  # Cer(d18:1/16:0)
        "ceramide",  # Ceramide C16:0
        " cer ",  # space cer space
        " cer(",  # space cer(
        "cer ",  # Cer 24:1
        "cer-c",  # some naming variants
    ]

    for p in patterns:
        if p in lower:
            return True

    return False


def resolve_csv_path(study_id: str | None, csv_path: str | None) -> Path | None:
    """
    Determine which CSV path to use:
    - If csv_path is provided, use it.
    - Else if study_id is provided, use data/mwtab/<study_id>.csv.
    - Else return None.

    Args:
        study_id: Optional MW study ID
        csv_path: Optional explicit CSV file path

    Returns:
        Path object if valid CSV found, None otherwise
    """
    if csv_path:
        p = Path(csv_path)
        if p.is_file():
            return p
        print(f"CSV path does not exist: {csv_path}")
        return None

    if study_id:
        p = Path("data/mwtab") / f"{study_id.upper().strip()}.csv"
        if p.is_file():
            return p
        print(f"CSV for study_id {study_id} not found at {p}")
        return None

    return None


def scan_ceramides(csv_path: Path):
    """
    Scan the CSV file for ceramide-like metabolite names.

    Heuristic:
    - Look for columns named like 'metabolite_name', 'Metabolite', 'METABOLITE_NAME',
      'refmet_name', 'REFMET_NAME', etc.
    - If none match, assume the first column is metabolite name.
    - Ceramide-like patterns: contains 'cer(d', 'ceramide', 'cer ', ' cer(' (case-insensitive).

    Args:
        csv_path: Path to CSV file to scan

    Returns:
        Tuple of (Counter of ceramide names, total row count)
    """
    name_candidates = [
        "metabolite_name",
        "METABOLITE_NAME",
        "Metabolite",
        "METABOLITE",
        "refmet_name",
        "REFMET_NAME",
        "compound",
        "COMPOUND",
    ]

    cer_hits = Counter()
    total_rows = 0

    with csv_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.reader(f)
        try:
            header = next(reader)
        except StopIteration:
            print("CSV appears to be empty.")
            return cer_hits, total_rows

        # Find name column index
        header_lower = [h.lower() for h in header]
        name_idx = None
        for candidate in name_candidates:
            candidate_lower = candidate.lower()
            if candidate_lower in header_lower:
                name_idx = header_lower.index(candidate_lower)
                break

        if name_idx is None:
            # Fallback: assume first column is metabolite name
            name_idx = 0
            print(
                f"Warning: No metabolite name column found, using first column: '{header[0] if header else 'unknown'}'"
            )

        for row in reader:
            if not row or len(row) <= name_idx:
                continue
            total_rows += 1
            metab_name = row[name_idx].strip()
            if is_ceramide_like(metab_name):
                cer_hits[metab_name] += 1

    return cer_hits, total_rows


def print_summary(cer_hits, total_rows: int):
    """
    Print summary of ceramide hits to stdout.

    Args:
        cer_hits: Counter of ceramide metabolite names
        total_rows: Total number of data rows scanned
    """
    print(f"Scanned {total_rows} data rows.")
    if not cer_hits:
        print("No ceramide-like metabolites found.")
        return

    print(f"\nFound {len(cer_hits)} unique ceramide-like metabolite(s):")
    print("Ceramide-like hits:")
    for name, count in cer_hits.most_common():
        print(f"  {name} : {count}")


def write_report(path_str: str, cer_hits, total_rows: int, csv_path: Path):
    """
    Write a text report of ceramide hits to a file.

    Args:
        path_str: Path to output report file
        cer_hits: Counter of ceramide metabolite names
        total_rows: Total number of data rows scanned
        csv_path: Path to source CSV file
    """
    path = Path(path_str)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as f:
        f.write(f"Source CSV: {csv_path}\n")
        f.write(f"Scanned rows: {total_rows}\n\n")
        if not cer_hits:
            f.write("No ceramide-like metabolites found.\n")
            return
        f.write(f"Found {len(cer_hits)} unique ceramide-like metabolite(s):\n")
        f.write("Ceramide-like hits:\n")
        for name, count in cer_hits.most_common():
            f.write(f"{name}\t{count}\n")

    print(f"\n✅ Report saved to {path}")


def main():
    parser = argparse.ArgumentParser(
        description="Scan mwTab CSV for ceramide-like metabolites and summarize hits."
    )
    parser.add_argument(
        "--study-id",
        required=False,
        help="Metabolomics Workbench study ID (e.g. ST004396). "
        "If provided without --csv-path, CSV is assumed at data/mwtab/<study_id>.csv",
    )
    parser.add_argument(
        "--csv-path",
        required=False,
        help="Path to an mwTab CSV file. If not provided, study_id is used.",
    )
    parser.add_argument(
        "--output-report",
        required=False,
        help="Optional path to save a text summary report.",
    )
    args = parser.parse_args()

    if not args.study_id and not args.csv_path:
        parser.error("Either --study-id or --csv-path must be provided.")

    csv_path = resolve_csv_path(args.study_id, args.csv_path)
    if csv_path is None:
        return 1

    print(f"Scanning CSV: {csv_path}")
    ceramide_hits, total_rows = scan_ceramides(csv_path)
    print_summary(ceramide_hits, total_rows)

    if args.output_report:
        write_report(args.output_report, ceramide_hits, total_rows, csv_path)

    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
