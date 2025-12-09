#!/usr/bin/env python3
"""
Bulk ingest lipid signatures from a directory of TSV/CSV files.

This script scans a configured directory for signature files and ingests
them all into Notion databases:
- Lipid Signatures database
- Lipid Signature Components database
- Lipid Species database

Usage:
    python scripts/bulk_ingest_signatures.py
    python scripts/bulk_ingest_signatures.py --signatures-dir /path/to/signatures
"""

import argparse
import sys
from pathlib import Path
from typing import Any, Dict, List, Set

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.config import get_config
from amprenta_rag.ingestion.signature_ingestion import \
    ingest_signature_from_file
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def find_signature_files(signatures_dir: Path) -> List[Path]:
    """
    Find all signature TSV/CSV files in a directory.

    Args:
        signatures_dir: Directory to scan

    Returns:
        List of Path objects for signature files
    """
    signature_files: List[Path] = []

    if not signatures_dir.exists():
        logger.warning(
            "[INGEST][SIGNATURES][BULK] Directory does not exist: %s",
            signatures_dir,
        )
        return signature_files

    if not signatures_dir.is_dir():
        logger.warning(
            "[INGEST][SIGNATURES][BULK] Path is not a directory: %s",
            signatures_dir,
        )
        return signature_files

    # Find all TSV and CSV files
    for pattern in ["*.tsv", "*.csv"]:
        signature_files.extend(signatures_dir.glob(pattern))

    # Filter out hidden files and directories
    signature_files = [
        f for f in signature_files if f.is_file() and not f.name.startswith(".")
    ]

    # Sort for consistent processing order
    signature_files.sort()

    return signature_files


def extract_signature_metadata_from_filename(file_path: Path) -> Dict[str, Any]:
    """
    Extract optional metadata from filename patterns.

    Attempts to infer signature type, version, etc. from filename.

    Args:
        file_path: Path to signature file

    Returns:
        Dictionary with extracted metadata (if any)
    """
    metadata: Dict[str, Any] = {}
    filename = file_path.stem.lower()

    # Try to extract version from filename (e.g., "signature_v1.tsv" -> "1")
    import re

    version_match = re.search(r"[vV](\d+(?:\.\d+)?)", filename)
    if version_match:
        metadata["version"] = version_match.group(1)

    # Try to infer signature type from filename
    if "consortium" in filename or "blsa" in filename or "adni" in filename:
        metadata["signature_type"] = "Consortium"
    elif "literature" in filename or "paper" in filename:
        metadata["signature_type"] = "Literature-derived"
    elif "dataset" in filename or "mw" in filename:
        metadata["signature_type"] = "Open Dataset"
    else:
        metadata["signature_type"] = "Literature-derived"  # Default

    # Try to extract disease context from filename
    disease_keywords = ["als", "ad", "alzheimer", "parkinson", "pd", "fxs", "fragile"]
    disease_context: List[str] = []
    for keyword in disease_keywords:
        if keyword in filename:
            if keyword == "als":
                disease_context.append("ALS")
            elif keyword in ["ad", "alzheimer"]:
                disease_context.append("Alzheimer's disease")
            elif keyword in ["pd", "parkinson"]:
                disease_context.append("Parkinson's disease")
            elif keyword in ["fxs", "fragile"]:
                disease_context.append("Fragile X Syndrome")

    if disease_context:
        metadata["disease_context"] = disease_context

    # Try to extract matrix from filename
    matrix_keywords = ["csf", "plasma", "serum", "tissue", "cell"]
    matrix: List[str] = []
    for keyword in matrix_keywords:
        if keyword in filename:
            if keyword == "csf":
                matrix.append("CSF")
            elif keyword == "plasma":
                matrix.append("Plasma")
            elif keyword == "serum":
                matrix.append("Serum")
            elif keyword == "tissue":
                matrix.append("Tissue")
            elif keyword == "cell":
                matrix.append("Cell")

    if matrix:
        metadata["matrix"] = matrix

    return metadata


def bulk_ingest_signatures(signatures_dir: Path) -> Dict[str, Any]:
    """
    Bulk ingest all signature files from a directory.

    Args:
        signatures_dir: Directory containing signature TSV/CSV files

    Returns:
        Dictionary with summary statistics:
        - files_processed: Number of files processed
        - files_failed: Number of files that failed
        - signatures_created: Number of new signatures created
        - signatures_updated: Number of existing signatures updated
        - components_created: Total components created/updated
        - species_created: Total lipid species created/linked
        - warnings: List of all warnings
        - errors: List of all errors
    """
    logger.info(
        "[INGEST][SIGNATURES][BULK] Starting bulk ingestion from directory: %s",
        signatures_dir,
    )

    # Find all signature files
    signature_files = find_signature_files(signatures_dir)

    if not signature_files:
        logger.warning(
            "[INGEST][SIGNATURES][BULK] No signature files found in directory: %s",
            signatures_dir,
        )
        return {
            "files_processed": 0,
            "files_failed": 0,
            "signatures_created": 0,
            "signatures_updated": 0,
            "components_created": 0,
            "species_created": 0,
            "warnings": [],
            "errors": [],
        }

    logger.info(
        "[INGEST][SIGNATURES][BULK] Found %d signature file(s) to process",
        len(signature_files),
    )

    # Statistics
    files_processed = 0
    files_failed = 0
    signatures_created = 0
    signatures_updated = 0
    total_components = 0
    total_species = 0
    all_warnings: List[str] = []
    all_errors: List[str] = []
    processed_signature_ids: Set[str] = set()

    # Process each signature file
    for sig_file in signature_files:
        logger.info(
            "[INGEST][SIGNATURES][BULK] Processing signature file: %s",
            sig_file.name,
        )

        try:
            # Extract metadata from filename (optional hints)
            file_metadata = extract_signature_metadata_from_filename(sig_file)

            # Ingest signature
            result = ingest_signature_from_file(
                signature_path=sig_file,
                signature_type=file_metadata.get(
                    "signature_type", "Literature-derived"
                ),
                data_ownership="Public",
                version=file_metadata.get("version"),
                description=None,
                disease_context=file_metadata.get("disease_context"),
                matrix=file_metadata.get("matrix"),
            )

            if result["signature_page_id"]:
                # Track if this is a new signature or existing
                if result["signature_page_id"] in processed_signature_ids:
                    # This signature was already processed (shouldn't happen with current logic)
                    signatures_updated += 1
                else:
                    # Check if signature was created or found existing
                    # We can't easily distinguish, so we'll count as processed
                    signatures_created += 1
                    processed_signature_ids.add(result["signature_page_id"])

                total_components += result["component_count"]
                total_species += result["species_count"]

                files_processed += 1

                logger.info(
                    "[INGEST][SIGNATURES][BULK] Successfully processed %s: "
                    "%d components, %d species",
                    sig_file.name,
                    result["component_count"],
                    result["species_count"],
                )

                # Collect warnings
                if result.get("warnings"):
                    for warning in result["warnings"]:
                        all_warnings.append(f"{sig_file.name}: {warning}")
            else:
                # Signature page creation failed
                error_msg = f"{sig_file.name}: Failed to create signature page"
                logger.error(f"[INGEST][SIGNATURES][BULK] {error_msg}")
                all_errors.append(error_msg)
                files_failed += 1

        except Exception as e:
            error_msg = f"{sig_file.name}: {str(e)}"
            logger.error(
                "[INGEST][SIGNATURES][BULK] Error processing %s: %r",
                sig_file.name,
                e,
            )
            all_errors.append(error_msg)
            files_failed += 1

    logger.info(
        "[INGEST][SIGNATURES][BULK] Bulk ingestion complete: "
        "%d files processed, %d failed",
        files_processed,
        files_failed,
    )

    return {
        "files_processed": files_processed,
        "files_failed": files_failed,
        "signatures_created": signatures_created,
        "signatures_updated": signatures_updated,
        "components_created": total_components,
        "species_created": total_species,
        "warnings": all_warnings,
        "errors": all_errors,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Bulk ingest lipid signatures from a directory of TSV/CSV files."
    )
    parser.add_argument(
        "--signatures-dir",
        type=Path,
        default=None,
        help="Directory containing signature files (overrides SIGNATURES_DIR from config)",
    )

    args = parser.parse_args()

    # Determine signatures directory
    if args.signatures_dir:
        signatures_dir = args.signatures_dir.resolve()
    else:
        cfg = get_config()
        signatures_dir_str = cfg.pipeline.signatures_dir

        if not signatures_dir_str:
            print(
                "Error: SIGNATURES_DIR is not configured. "
                "Set it in your .env file or use --signatures-dir option.",
                file=sys.stderr,
            )
            print(
                "\nExample:",
                file=sys.stderr,
            )
            print("  SIGNATURES_DIR=data/signatures", file=sys.stderr)
            print("\nOr use --signatures-dir:", file=sys.stderr)
            print(
                "  python scripts/bulk_ingest_signatures.py --signatures-dir /path/to/signatures",
                file=sys.stderr,
            )
            return 1

        signatures_dir = Path(signatures_dir_str).resolve()

    if not signatures_dir.exists():
        print(
            f"Error: Signatures directory does not exist: {signatures_dir}",
            file=sys.stderr,
        )
        return 1

    if not signatures_dir.is_dir():
        print(f"Error: Path is not a directory: {signatures_dir}", file=sys.stderr)
        return 1

    print(f"\n{'=' * 80}")
    print("Bulk Signature Ingestion")
    print(f"{'=' * 80}")
    print(f"\nDirectory: {signatures_dir}")
    print()

    # Run bulk ingestion
    try:
        summary = bulk_ingest_signatures(signatures_dir)

        # Print summary
        print(f"\n{'=' * 80}")
        print("Summary")
        print(f"{'=' * 80}")
        print(f"\nFiles processed:      {summary['files_processed']}")
        print(f"Files failed:          {summary['files_failed']}")
        print(f"Signatures created:    {summary['signatures_created']}")
        print(f"Signatures updated:    {summary['signatures_updated']}")
        print(f"Components created:    {summary['components_created']}")
        print(f"Lipid species created: {summary['species_created']}")

        if summary["warnings"]:
            print(f"\n⚠️  Warnings ({len(summary['warnings'])}):")
            for warning in summary["warnings"][:10]:  # Show first 10
                print(f"  - {warning}")
            if len(summary["warnings"]) > 10:
                print(f"  ... and {len(summary['warnings']) - 10} more warnings")

        if summary["errors"]:
            print(f"\n❌ Errors ({len(summary['errors'])}):")
            for error in summary["errors"]:
                print(f"  - {error}")

        print(f"\n{'=' * 80}\n")

        return 0 if summary["files_failed"] == 0 else 1

    except Exception as e:
        print(f"\n❌ Fatal error during bulk ingestion: {e}", file=sys.stderr)
        logger.exception("Fatal error in bulk ingestion")
        return 1


if __name__ == "__main__":
    sys.exit(main())
