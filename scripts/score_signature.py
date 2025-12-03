#!/usr/bin/env python3
"""
CLI script to score a lipid signature against a dataset.

Loads a signature definition (TSV) and a dataset (CSV), performs species matching,
and computes a direction-aware, weight-aware signature score.

Usage:
    python scripts/score_signature.py \
      --signature-tsv data/signatures/ALS_CSF_Core_6Ceramides.tsv \
      --dataset-csv data/mwtab/ST004396.csv \
      --output-report data/reports/ST004396_ALS_signature_score.txt
"""

import argparse
import csv
import json
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.signatures import (load_signature_from_tsv, match_species,
                                     score_signature)
from scripts.scan_sphingolipids_in_mwtab_csv import scan_sphingolipids


def load_dataset_species(csv_path: Path) -> set[str]:
    """
    Load species names from a dataset CSV.

    Uses the sphingolipid scanner to identify sphingolipid species,
    or falls back to reading all metabolite names.

    Args:
        csv_path: Path to dataset CSV file

    Returns:
        Set of species names found in the dataset
    """
    # Try to use sphingolipid scanner to get sphingolipids
    try:
        _, species_counts, _ = scan_sphingolipids(csv_path, refmet_map={})

        # Collect all species from all classes
        all_species = set()
        for class_counter in species_counts.values():
            all_species.update(class_counter.keys())

        if all_species:
            return all_species
    except Exception as e:
        print(f"Warning: Could not use sphingolipid scanner: {e}")
        print("Falling back to reading all metabolite names from CSV...")

    # Fallback: read all metabolite names from CSV
    species = set()

    with csv_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.reader(f)
        try:
            header = next(reader)
        except StopIteration:
            return species

        header_lower = [h.lower() for h in header]

        # Find metabolite name column
        name_col_candidates = [
            "metabolite_name",
            "metabolite",
            "refmet_name",
            "compound",
            "species",
        ]

        name_idx = None
        for c in name_col_candidates:
            if c in header_lower:
                name_idx = header_lower.index(c)
                break

        if name_idx is None:
            name_idx = 0  # fallback to first column

        for row in reader:
            if row and len(row) > name_idx:
                species_name = row[name_idx].strip()
                if species_name:
                    species.add(species_name)

    return species


def main():
    parser = argparse.ArgumentParser(
        description="Score a lipid signature against a dataset CSV."
    )
    parser.add_argument(
        "--signature-tsv",
        required=True,
        help="Path to signature TSV file (columns: species, direction, weight).",
    )
    parser.add_argument(
        "--dataset-csv",
        required=True,
        help="Path to dataset CSV file (from convert_mwtab_to_csv.py or similar).",
    )
    parser.add_argument(
        "--output-report",
        required=False,
        help="Optional path to save a detailed report (TSV or JSON).",
    )
    parser.add_argument(
        "--refmet-map",
        required=False,
        help="Optional path to RefMet mapping TSV for species matching.",
    )
    parser.add_argument(
        "--format",
        choices=["text", "tsv", "json"],
        default="text",
        help="Output format for report (default: text).",
    )
    args = parser.parse_args()

    # Load signature
    signature_path = Path(args.signature_tsv)
    if not signature_path.is_file():
        print(f"Error: Signature file not found: {signature_path}")
        return 1

    print(f"Loading signature from {signature_path}...")
    try:
        signature = load_signature_from_tsv(signature_path)
        print(
            f"Loaded signature '{signature.name}' with {len(signature.components)} components"
        )
    except Exception as e:
        print(f"Error loading signature: {e}")
        return 1

    # Load dataset
    dataset_path = Path(args.dataset_csv)
    if not dataset_path.is_file():
        print(f"Error: Dataset file not found: {dataset_path}")
        return 1

    print(f"Loading dataset from {dataset_path}...")
    dataset_species = load_dataset_species(dataset_path)
    print(f"Found {len(dataset_species)} species in dataset")

    # Load RefMet map if provided
    refmet_map = {}
    if args.refmet_map:
        refmet_map_path = Path(args.refmet_map)
        if refmet_map_path.is_file():
            print(f"Loading RefMet map from {refmet_map_path}...")
            # Simple RefMet map loader (can be enhanced)
            with refmet_map_path.open("r", encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    if "refmet_id" in row and "refmet_name" in row:
                        refmet_map[row["refmet_id"].lower()] = row["refmet_name"]
                        refmet_map[row["refmet_name"].lower()] = row["refmet_name"]
            print(f"Loaded {len(refmet_map)} RefMet mappings")
        else:
            print(f"Warning: RefMet map file not found: {refmet_map_path}")

    # Score signature
    print("\nScoring signature...")
    result = score_signature(
        signature=signature,
        dataset_species=dataset_species,
        dataset_directions=None,  # TODO: Support direction inference from multi-condition data
        refmet_map=refmet_map if refmet_map else None,
    )

    # Print summary
    print(f"\n{'='*60}")
    print(
        f"Signature Score: {result.total_score:.3f} (0.0 = no match, 1.0 = perfect match)"
    )
    print(f"{'='*60}")
    print(
        f"\nMatched: {len(result.matched_species)}/{len(signature.components)} components"
    )
    print(f"Missing: {len(result.missing_species)} components")
    print(f"Conflicts: {len(result.conflicting_species)} components")

    if result.missing_species:
        print(f"\nMissing species:")
        for species in result.missing_species:
            print(f"  - {species}")

    if result.conflicting_species:
        print(f"\nConflicting species (direction mismatch):")
        for species in result.conflicting_species:
            print(f"  - {species}")

    print(f"\nComponent details:")
    print(
        f"{'Species':<30} {'Matched':<30} {'Match':<12} {'Direction':<12} {'Weight':<8}"
    )
    print("-" * 92)
    for comp in result.component_matches:
        matched = comp.matched_dataset_species or "N/A"
        match_type = comp.match_type
        direction = comp.direction_match
        weight = comp.weight

        print(
            f"{comp.signature_species[:28]:<30} {matched[:28]:<30} {match_type:<12} {direction:<12} {weight:<8.2f}"
        )

    # Save report if requested
    if args.output_report:
        report_path = Path(args.output_report)
        report_path.parent.mkdir(parents=True, exist_ok=True)

        if args.format == "json":
            # JSON report
            report_data = {
                "signature_name": signature.name,
                "dataset_path": str(dataset_path),
                "total_score": result.total_score,
                "matched_count": len(result.matched_species),
                "missing_count": len(result.missing_species),
                "conflict_count": len(result.conflicting_species),
                "missing_species": result.missing_species,
                "conflicting_species": result.conflicting_species,
                "matched_species": result.matched_species,
                "component_matches": [
                    {
                        "signature_species": comp.signature_species,
                        "matched_dataset_species": comp.matched_dataset_species,
                        "match_type": comp.match_type,
                        "direction_match": comp.direction_match,
                        "weight": comp.weight,
                    }
                    for comp in result.component_matches
                ],
            }

            with report_path.open("w", encoding="utf-8") as f:
                json.dump(report_data, f, indent=2)

        elif args.format == "tsv":
            # TSV report
            with report_path.open("w", encoding="utf-8", newline="") as f:
                writer = csv.writer(f, delimiter="\t")
                writer.writerow(
                    [
                        "signature_species",
                        "matched_dataset_species",
                        "match_type",
                        "direction_match",
                        "weight",
                    ]
                )
                for comp in result.component_matches:
                    writer.writerow(
                        [
                            comp.signature_species,
                            comp.matched_dataset_species or "",
                            comp.match_type,
                            comp.direction_match,
                            comp.weight,
                        ]
                    )

        else:
            # Text report
            with report_path.open("w", encoding="utf-8") as f:
                f.write(f"Signature Scoring Report\n")
                f.write(f"{'='*60}\n\n")
                f.write(f"Signature: {signature.name}\n")
                f.write(f"Dataset: {dataset_path}\n")
                f.write(f"Total Score: {result.total_score:.3f}\n\n")
                f.write(f"Summary:\n")
                f.write(
                    f"  Matched: {len(result.matched_species)}/{len(signature.components)}\n"
                )
                f.write(f"  Missing: {len(result.missing_species)}\n")
                f.write(f"  Conflicts: {len(result.conflicting_species)}\n\n")

                if result.missing_species:
                    f.write(f"Missing Species:\n")
                    for species in result.missing_species:
                        f.write(f"  - {species}\n")
                    f.write("\n")

                if result.conflicting_species:
                    f.write(f"Conflicting Species:\n")
                    for species in result.conflicting_species:
                        f.write(f"  - {species}\n")
                    f.write("\n")

                f.write(f"Component Details:\n")
                f.write(
                    f"{'Species':<30} {'Matched':<30} {'Match':<12} {'Direction':<12} {'Weight':<8}\n"
                )
                f.write("-" * 92 + "\n")
                for comp in result.component_matches:
                    matched = comp.matched_dataset_species or "N/A"
                    f.write(
                        f"{comp.signature_species[:28]:<30} {matched[:28]:<30} {comp.match_type:<12} {comp.direction_match:<12} {comp.weight:<8.2f}\n"
                    )

        print(f"\n✅ Report saved to {report_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
