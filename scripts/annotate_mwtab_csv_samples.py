#!/usr/bin/env python3
"""
Annotate mwTab CSV sample columns with metadata from mwTab JSON.

Enriches CSV column headers with human-readable sample metadata like:
- Group (FXS, TD, etc.)
- Disease
- Sample source/type
- Other relevant factors

Usage:
    python scripts/annotate_mwtab_csv_samples.py --study-id ST004396
    python scripts/annotate_mwtab_csv_samples.py --csv-path data/mwtab/ST004396.csv
"""

import argparse
import csv
import json
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.harvest_mw_studies import fetch_mw_mwtab


def resolve_csv_path(study_id: str | None, csv_path: str | None) -> Path | None:
    """
    Determine which CSV path to use.

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


def load_mwtab_json(study_id: str | None, mwtab_json_path: str | None) -> dict | None:
    """
    Load mwTab JSON data from file or API.

    Args:
        study_id: Optional MW study ID (used if mwtab_json_path not provided)
        mwtab_json_path: Optional path to mwTab JSON file

    Returns:
        Parsed JSON dict or None if error
    """
    if mwtab_json_path:
        p = Path(mwtab_json_path)
        if not p.is_file():
            print(f"mwTab JSON path does not exist: {mwtab_json_path}")
            return None
        with p.open("r", encoding="utf-8") as f:
            try:
                return json.load(f)
            except json.JSONDecodeError as e:
                print(f"Failed to parse mwTab JSON file: {e}")
                return None

    if study_id:
        text = fetch_mw_mwtab(study_id.upper().strip())
        if not text:
            print(f"No mwTab content returned for study {study_id}")
            return None

        # Try to parse JSON (handle potential multiple JSON objects)
        text_stripped = text.strip()
        first_brace = text_stripped.find("{")

        if first_brace == -1:
            print(f"mwTab text for {study_id} doesn't appear to be JSON")
            return None

        # Try to parse first complete JSON object
        for end_pos in range(len(text_stripped), first_brace, -1):
            try:
                json_str = text_stripped[first_brace:end_pos]
                return json.loads(json_str)
            except (json.JSONDecodeError, ValueError):
                continue

        print(f"Failed to parse mwTab JSON for {study_id}")
        return None

    return None


def build_sample_metadata_map(mwtab_data: dict) -> dict[str, str]:
    """
    Build a mapping from sample ID (e.g., 'FXS1') to a human-readable annotation string.

    Example output:
      'FXS1' -> 'FXS1 (Group=FXS, Source=LCL, Disease=Fragile X syndrome)'

    Args:
        mwtab_data: Parsed mwTab JSON data

    Returns:
        Dictionary mapping sample_id -> annotated_label
    """
    sample_meta = {}

    # Look for SUBJECT_SAMPLE_FACTORS section
    # It can be either a direct list or a dict with "Data" key
    factors_section = mwtab_data.get("SUBJECT_SAMPLE_FACTORS")

    if not factors_section:
        # Try alternative key names
        factors_section = mwtab_data.get("SUBJECT_SAMPLE_FACTORS_DATA")

    # Handle both list and dict structures
    data_list = None
    if isinstance(factors_section, list):
        data_list = factors_section
    elif isinstance(factors_section, dict):
        data_list = factors_section.get("Data") or factors_section.get("data")

    if not isinstance(data_list, list):
        return sample_meta

    for entry in data_list:
        if not isinstance(entry, dict):
            continue

        # Extract sample ID
        sample_id = (
            entry.get("Sample ID")
            or entry.get("Sample_ID")
            or entry.get("sample_id")
            or entry.get("SampleID")
        )

        if not sample_id:
            continue

        # Extract factors
        factors = (
            entry.get("Factors") or entry.get("FACTOR") or entry.get("factor") or {}
        )

        if not isinstance(factors, dict):
            factors = {}

        # Extract key fields for the label
        disease = (
            factors.get("Disease")
            or factors.get("disease")
            or factors.get("DISEASE")
            or ""
        )

        source = (
            factors.get("Sample source")
            or factors.get("Sample Source")
            or factors.get("sample_source")
            or factors.get("Source")
            or ""
        )

        genotype = (
            factors.get("Genotype")
            or factors.get("genotype")
            or factors.get("GENOTYPE")
            or ""
        )

        group = None

        # Heuristic: derive group from disease, genotype, or sample_id prefix
        disease_lower = disease.lower() if disease else ""
        genotype_lower = genotype.lower() if genotype else ""
        sample_id_upper = sample_id.upper()

        if (
            "fragile" in disease_lower
            or "fxs" in genotype_lower
            or "fxs" in disease_lower
        ):
            group = "FXS"
        elif (
            "control" in disease_lower
            or "typical" in disease_lower
            or "td" in genotype_lower
        ):
            group = "TD"
        elif "case" in disease_lower:
            group = "Case"
        elif "control" in sample_id_upper:
            group = "Control"

        # Fallback: prefix of sample_id
        if group is None:
            if sample_id_upper.startswith("FXS"):
                group = "FXS"
            elif sample_id_upper.startswith("TD"):
                group = "TD"
            elif sample_id_upper.startswith("CASE"):
                group = "Case"
            elif sample_id_upper.startswith("CONTROL"):
                group = "Control"

        # Build label parts
        label_parts = []

        if group:
            label_parts.append(f"Group={group}")

        if genotype:
            label_parts.append(f"Genotype={genotype}")

        if disease:
            # Truncate long disease names
            disease_short = disease[:30] + "..." if len(disease) > 30 else disease
            label_parts.append(f"Disease={disease_short}")

        if source:
            # Truncate long source descriptions
            source_short = source[:20] + "..." if len(source) > 20 else source
            label_parts.append(f"Source={source_short}")

        # Combine into annotation
        label_suffix = ", ".join(label_parts) if label_parts else ""
        if label_suffix:
            sample_meta[sample_id] = f"{sample_id} ({label_suffix})"
        else:
            sample_meta[sample_id] = sample_id

    return sample_meta


def annotate_csv_headers(
    csv_path: Path, sample_meta_map: dict[str, str], output_csv: str | None
):
    """
    Annotate CSV column headers with sample metadata.

    Args:
        csv_path: Path to input CSV file
        sample_meta_map: Dictionary mapping sample_id -> annotated_label
        output_csv: Optional output path (defaults to <input>_annotated.csv)
    """
    if output_csv:
        out_path = Path(output_csv)
    else:
        out_path = csv_path.with_name(csv_path.stem + "_annotated.csv")

    with csv_path.open("r", encoding="utf-8", newline="") as f_in, out_path.open(
        "w", encoding="utf-8", newline=""
    ) as f_out:

        reader = csv.reader(f_in)
        writer = csv.writer(f_out)

        try:
            header = next(reader)
        except StopIteration:
            print("CSV is empty; nothing to annotate.")
            return

        new_header = []
        for col in header:
            if col in sample_meta_map:
                new_header.append(sample_meta_map[col])
            else:
                new_header.append(col)

        writer.writerow(new_header)

        # Copy all data rows unchanged
        for row in reader:
            writer.writerow(row)

    print(f"✅ Annotated CSV saved to: {out_path}")
    print(f"   Original columns: {len([h for h in header if h != 'Metabolite'])}")
    print(f"   Annotated columns: {len([h for h in new_header if h != 'Metabolite'])}")


def main():
    parser = argparse.ArgumentParser(
        description="Annotate mwTab CSV sample columns with metadata from mwTab JSON."
    )
    parser.add_argument(
        "--study-id",
        required=False,
        help="Metabolomics Workbench study ID (e.g., ST004396). Used to locate CSV & mwTab JSON if paths not provided.",
    )
    parser.add_argument(
        "--csv-path",
        required=False,
        help="Path to existing mwTab CSV (e.g., data/mwtab/ST004396.csv).",
    )
    parser.add_argument(
        "--mwtab-json",
        required=False,
        help="Optional path to a mwTab JSON file. If not provided, fetch_mw_mwtab(study_id) is used.",
    )
    parser.add_argument(
        "--output-csv",
        required=False,
        help="Path to write annotated CSV. Defaults to <input>_annotated.csv.",
    )
    args = parser.parse_args()

    if not args.study_id and not args.csv_path:
        parser.error("Either --study-id or --csv-path must be provided.")

    csv_path = resolve_csv_path(args.study_id, args.csv_path)
    if csv_path is None:
        return 1

    study_id = args.study_id
    if not study_id and args.csv_path:
        # Try to extract study ID from CSV filename
        csv_name = Path(args.csv_path).stem
        if csv_name.upper().startswith("ST"):
            study_id = csv_name.upper()

    print(f"Loading mwTab JSON for annotation...")
    mwtab_data = load_mwtab_json(study_id, args.mwtab_json)
    if mwtab_data is None:
        return 1

    print(f"Extracting sample metadata...")
    sample_meta_map = build_sample_metadata_map(mwtab_data)

    if not sample_meta_map:
        print(
            "⚠️  No sample metadata found in mwTab JSON. CSV headers will not be annotated."
        )
    else:
        print(f"Found metadata for {len(sample_meta_map)} samples:")
        for sample_id, annotation in list(sample_meta_map.items())[:5]:
            print(f"  {sample_id} -> {annotation}")
        if len(sample_meta_map) > 5:
            print(f"  ... and {len(sample_meta_map) - 5} more")

    print(f"\nAnnotating CSV: {csv_path}")
    annotate_csv_headers(csv_path, sample_meta_map, args.output_csv)

    return 0


if __name__ == "__main__":
    sys.exit(main())
