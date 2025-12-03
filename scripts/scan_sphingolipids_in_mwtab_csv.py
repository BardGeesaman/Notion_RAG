#!/usr/bin/env python3
"""
Hybrid sphingolipid scanner for mwTab CSV files.

Scans metabolite names using rule-based patterns and optional RefMet mapping
to classify sphingolipids into classes (Ceramide, SM, HexCer, LacCer, etc.)

Usage:
    python scripts/scan_sphingolipids_in_mwtab_csv.py --study-id ST004396
    python scripts/scan_sphingolipids_in_mwtab_csv.py --csv-path data/mwtab/ST004396.csv
    python scripts/scan_sphingolipids_in_mwtab_csv.py --csv-path data/mwtab/ST004396.csv --refmet-map data/refmet_classes.tsv
"""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


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


def load_refmet_map(tsv_path: Path) -> dict[str, str]:
    """
    Load a RefMet mapping TSV with columns:
        refmet_id, refmet_name, lipid_class

    Returns a dict keyed by both refmet_id and refmet_name (lower-cased)
    for flexible lookups.

    Args:
        tsv_path: Path to RefMet mapping TSV file

    Returns:
        Dictionary mapping RefMet IDs/names (lowercase) to lipid class labels
    """
    refmet_map = {}

    with tsv_path.open("r", encoding="utf-8") as f:
        # Allow tab or comma; simplest is to split on tab.
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                continue

            ref_id, ref_name, lipid_class = parts[0], parts[1], parts[2]

            if ref_id:
                refmet_map[ref_id.lower()] = lipid_class
            if ref_name:
                refmet_map[ref_name.lower()] = lipid_class

    return refmet_map


def classify_sphingolipid_by_name(name: str) -> str | None:
    """
    Rule-based classification of sphingolipid species into classes:

    Ceramide, Dihydroceramide, HexCer, LacCer, SM, Sphingosine, S1P,
    Deoxysphingolipids, Ceramide-1-phosphate, Gangliosides, etc.

    Args:
        name: Metabolite name to classify

    Returns:
        Class label string or None if not recognized as a sphingolipid
    """
    if not name:
        return None

    lower = name.lower()

    # Order matters: more specific patterns first.

    # Ceramide-1-phosphate
    if "ceramide-1-phosphate" in lower or "c1p" in lower:
        return "Ceramide-1-phosphate"

    # Dihydroceramide
    if "dihydroceramide" in lower or "dhcer" in lower:
        return "Dihydroceramide"

    # Hexosylceramide (HexCer, GlcCer, GalCer)
    if "hexcer" in lower or "glccer" in lower or "galcer" in lower:
        return "Hexosylceramide"

    # Lactosylceramide
    if "lacc" in lower or "laccer" in lower or "lactosylceramide" in lower:
        return "Lactosylceramide"

    # Sphingomyelin
    if "sphingomyelin" in lower or lower.startswith("sm(") or lower.startswith("sm "):
        return "Sphingomyelin"

    # Sphingosine-1-phosphate
    if "sphingosine-1-phosphate" in lower or "s1p" in lower:
        return "Sphingosine-1-phosphate"

    # Sphingosine (but not S1P)
    if "sphingosine" in lower and "phosphate" not in lower:
        return "Sphingosine"

    # Deoxysphingolipids
    if "deoxysphingolipid" in lower or "deoxycer" in lower or "m18:" in lower:
        return "Deoxysphingolipid"

    # Gangliosides
    if any(
        g in lower for g in ["gm1", "gm2", "gm3", "gd1", "gd2", "gd3", "gt1", "gt2"]
    ):
        return "Ganglioside"

    # General ceramide patterns
    if "ceramide" in lower or "cer(d" in lower or lower.startswith("cer "):
        return "Ceramide"

    return None


def scan_sphingolipids(
    csv_path: Path,
    refmet_map: dict[str, str],
    refmet_id_column: str | None = None,
    refmet_name_column: str | None = None,
):
    """
    Scan CSV rows and assign sphingolipid classes using:
    - Rule-based name patterns (first)
    - Optional RefMet mapping (if refmet_map provided)

    Args:
        csv_path: Path to CSV file to scan
        refmet_map: Dictionary mapping RefMet IDs/names to lipid classes
        refmet_id_column: Optional CSV column name for RefMet ID
        refmet_name_column: Optional CSV column name for RefMet name

    Returns:
        Tuple of (class_counts, species_counts, total_rows):
        - class_counts: Counter(class_label -> count)
        - species_counts: dict[class_label -> Counter(species_name -> count)]
        - total_rows: int
    """
    class_counts = Counter()
    species_counts: dict[str, Counter] = defaultdict(Counter)
    total_rows = 0

    with csv_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.reader(f)
        try:
            header = next(reader)
        except StopIteration:
            print("CSV appears to be empty.")
            return class_counts, species_counts, total_rows

        header_lower = [h.lower() for h in header]

        # Guess a metabolite name column
        name_col_candidates = [
            "metabolite_name",
            "metabolite",
            "refmet_name",
            "compound",
        ]
        name_idx = None
        for c in name_col_candidates:
            if c in header_lower:
                name_idx = header_lower.index(c)
                break
        if name_idx is None:
            name_idx = 0  # fallback
            print(
                f"Warning: No metabolite name column found, using first column: '{header[0] if header else 'unknown'}'"
            )

        # Optional RefMet ID / name columns
        refmet_id_idx = None
        refmet_name_idx = None
        if refmet_id_column and refmet_id_column.lower() in header_lower:
            refmet_id_idx = header_lower.index(refmet_id_column.lower())
        if refmet_name_column and refmet_name_column.lower() in header_lower:
            refmet_name_idx = header_lower.index(refmet_name_column.lower())

        for row in reader:
            if not row or len(row) <= name_idx:
                continue
            total_rows += 1
            metab_name = row[name_idx].strip()

            # 1) Rule-based classification by name
            cls = classify_sphingolipid_by_name(metab_name)

            # 2) If no class and we have a RefMet map + ID/name, try that
            if not cls and refmet_map:
                ref_ids = []
                if refmet_id_idx is not None and len(row) > refmet_id_idx:
                    ref_ids.append(row[refmet_id_idx].strip().lower())
                if refmet_name_idx is not None and len(row) > refmet_name_idx:
                    ref_ids.append(row[refmet_name_idx].strip().lower())

                for ref in ref_ids:
                    if ref in refmet_map:
                        cls = refmet_map[ref]
                        break

            if not cls:
                continue  # skip non-sphingolipid rows

            class_counts[cls] += 1
            species_counts[cls][metab_name] += 1

    return class_counts, species_counts, total_rows


def print_summary(class_counts, species_counts, total_rows: int):
    """
    Print summary of sphingolipid hits to stdout.

    Args:
        class_counts: Counter of lipid class labels
        species_counts: Dict of class_label -> Counter of species names
        total_rows: Total number of data rows scanned
    """
    print(f"Scanned {total_rows} data rows.")
    if not class_counts:
        print("No sphingolipid-like metabolites found.")
        return

    print(f"\nFound {len(class_counts)} unique sphingolipid class(es):")
    print("\nClass-level counts:")
    for cls, count in class_counts.most_common():
        print(f"  {cls}: {count}")

    print("\nTop species per class (up to 10):")
    for cls, counter in species_counts.items():
        print(f"\n[{cls}]")
        for name, count in counter.most_common(10):
            print(f"  {name} : {count}")


def write_report(
    path_str: str, class_counts, species_counts, total_rows: int, csv_path: Path
):
    """
    Write a text report of sphingolipid hits to a file.

    Args:
        path_str: Path to output report file
        class_counts: Counter of lipid class labels
        species_counts: Dict of class_label -> Counter of species names
        total_rows: Total number of data rows scanned
        csv_path: Path to source CSV file
    """
    path = Path(path_str)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as f:
        f.write(f"Source CSV: {csv_path}\n")
        f.write(f"Scanned rows: {total_rows}\n\n")
        if not class_counts:
            f.write("No sphingolipid-like metabolites found.\n")
            return

        f.write(f"Found {len(class_counts)} unique sphingolipid class(es):\n")
        f.write("Class-level counts:\n")
        for cls, count in class_counts.most_common():
            f.write(f"{cls}\t{count}\n")

        f.write("\nSpecies detail (top 10 per class):\n")
        for cls, counter in species_counts.items():
            f.write(f"\n[{cls}]\n")
            for name, count in counter.most_common(10):
                f.write(f"{name}\t{count}\n")

    print(f"\n✅ Report saved to {path}")


def main():
    parser = argparse.ArgumentParser(
        description="Scan mwTab CSV for sphingolipid-like metabolites (ceramide, SM, HexCer, LacCer, etc.) and summarize hits, using rule-based and optional RefMet mapping."
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
    parser.add_argument(
        "--refmet-map",
        required=False,
        help=(
            "Optional path to a RefMet mapping TSV file with columns: refmet_id,refmet_name,lipid_class. "
            "If provided and the CSV has a RefMet ID or name column, this will refine class assignments."
        ),
    )
    parser.add_argument(
        "--refmet-id-column",
        required=False,
        help="Optional CSV column name for RefMet ID (e.g. REFMET_ID).",
    )
    parser.add_argument(
        "--refmet-name-column",
        required=False,
        help="Optional CSV column name for RefMet name (e.g. REFMET_NAME).",
    )
    args = parser.parse_args()

    if not args.study_id and not args.csv_path:
        parser.error("Either --study-id or --csv-path must be provided.")

    csv_path = resolve_csv_path(args.study_id, args.csv_path)
    if csv_path is None:
        return 1

    refmet_map = {}
    if args.refmet_map:
        refmet_map_path = Path(args.refmet_map)
        if not refmet_map_path.is_file():
            print(f"Warning: RefMet map file not found: {args.refmet_map}")
        else:
            print(f"Loading RefMet map from {refmet_map_path}...")
            refmet_map = load_refmet_map(refmet_map_path)
            print(f"Loaded {len(refmet_map)} RefMet mappings")

    print(f"Scanning CSV: {csv_path}")
    class_counts, species_counts, total_rows = scan_sphingolipids(
        csv_path,
        refmet_map=refmet_map,
        refmet_id_column=args.refmet_id_column,
        refmet_name_column=args.refmet_name_column,
    )
    print_summary(class_counts, species_counts, total_rows)

    if args.output_report:
        write_report(
            args.output_report, class_counts, species_counts, total_rows, csv_path
        )

    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
