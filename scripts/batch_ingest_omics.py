#!/usr/bin/env python3
"""
Batch ingestion framework for multi-omics datasets.

Automatically detects omics type from files and ingests them in batch.

Usage:
    python scripts/batch_ingest_omics.py --directory <path>
    python scripts/batch_ingest_omics.py --directory <path> --omics-type lipidomics
    python scripts/batch_ingest_omics.py --file <path> --file <path> ...
    python scripts/batch_ingest_omics.py --directory <path> --parallel --max-workers 4
"""

import argparse
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.ingestion.lipidomics_ingestion import ingest_lipidomics_file
from amprenta_rag.ingestion.metabolomics_ingestion import ingest_metabolomics_file
from amprenta_rag.ingestion.proteomics_ingestion import ingest_proteomics_file
from amprenta_rag.ingestion.transcriptomics_ingestion import ingest_transcriptomics_file
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Supported file extensions
SUPPORTED_EXTENSIONS = {".csv", ".tsv", ".txt"}

# Omics type detection patterns
OMICS_PATTERNS = {
    "lipidomics": {
        "filename_keywords": ["lipid", "cer", "sphing", "sm", "hexcer", "laccer"],
        "column_keywords": ["lipid", "species", "cer", "sphing", "sm"],
        "content_patterns": ["cer(", "sm(", "pc(", "pe(", "ceramide", "sphingomyelin"],
    },
    "metabolomics": {
        "filename_keywords": ["metabol", "metab", "polar"],
        "column_keywords": ["metabolite", "compound", "metab"],
        "content_patterns": ["glutamate", "citrate", "glucose", "lactate"],
    },
    "proteomics": {
        "filename_keywords": ["protein", "proteo", "lfq"],
        "column_keywords": ["protein", "gene", "uniprot", "accession"],
        "content_patterns": ["P0", "P1", "Q", "ENS", "_HUMAN", "_MOUSE"],
    },
    "transcriptomics": {
        "filename_keywords": ["rna", "transcript", "dge", "deseq", "edge"],
        "column_keywords": ["gene", "geneid", "ensembl", "symbol"],
        "content_patterns": ["ENSG", "ENSMUS", "log2fc", "padj", "pvalue"],
    },
}


def detect_omics_type(file_path: str, force_type: Optional[str] = None) -> Optional[str]:
    """
    Auto-detect omics type from file path and content.

    Args:
        file_path: Path to the file
        force_type: If provided, skip detection and use this type

    Returns:
        Detected omics type (lipidomics, metabolomics, proteomics, transcriptomics) or None
    """
    if force_type:
        if force_type.lower() in OMICS_PATTERNS:
            return force_type.lower()
        logger.warning(
            "[BATCH][DETECT] Unknown omics type '%s', attempting auto-detection",
            force_type,
        )

    file_path_lower = str(file_path).lower()
    filename = Path(file_path).stem.lower()

    # Score each omics type based on filename and path
    scores: Dict[str, int] = defaultdict(int)

    for omics_type, patterns in OMICS_PATTERNS.items():
        # Check filename keywords
        for keyword in patterns["filename_keywords"]:
            if keyword in filename or keyword in file_path_lower:
                scores[omics_type] += 2

    # Try to read file headers for column-based detection
    try:
        import pandas as pd

        # Read first few rows to detect delimiter and get headers
        with open(file_path, "r", encoding="utf-8") as f:
            first_line = f.readline()
            delimiter = "\t" if "\t" in first_line else ","
            f.seek(0)

        df = pd.read_csv(file_path, delimiter=delimiter, encoding="utf-8", nrows=5)
        columns_lower = [str(col).lower() for col in df.columns]

        for omics_type, patterns in OMICS_PATTERNS.items():
            # Check column keywords
            for keyword in patterns["column_keywords"]:
                if any(keyword in col for col in columns_lower):
                    scores[omics_type] += 3

            # Check content patterns (sample first few rows)
            for col in df.columns:
                sample_values = df[col].astype(str).str.lower().head(10).tolist()
                for pattern in patterns["content_patterns"]:
                    if any(pattern in val for val in sample_values):
                        scores[omics_type] += 1

    except Exception as e:
        logger.debug(
            "[BATCH][DETECT] Could not read file %s for content-based detection: %r",
            file_path,
            e,
        )

    # Return highest scoring omics type
    if scores:
        best_type = max(scores.items(), key=lambda x: x[1])
        if best_type[1] > 0:
            logger.info(
                "[BATCH][DETECT] Detected omics type '%s' for file %s (score: %d)",
                best_type[0],
                file_path,
                best_type[1],
            )
            return best_type[0]

    logger.warning(
        "[BATCH][DETECT] Could not auto-detect omics type for file %s",
        file_path,
    )
    return None


def collect_files(
    directory: Optional[str] = None, files: Optional[List[str]] = None
) -> List[Path]:
    """
    Collect all supported files from directory or file list.

    Args:
        directory: Directory to scan
        files: Explicit list of file paths

    Returns:
        List of file paths
    """
    file_paths: List[Path] = []

    if directory:
        dir_path = Path(directory)
        if not dir_path.exists():
            raise ValueError(f"Directory does not exist: {directory}")

        if not dir_path.is_dir():
            raise ValueError(f"Path is not a directory: {directory}")

        # Recursively find all CSV/TSV files
        for ext in SUPPORTED_EXTENSIONS:
            file_paths.extend(dir_path.rglob(f"*{ext}"))

        logger.info(
            "[BATCH] Found %d files in directory %s",
            len(file_paths),
            directory,
        )

    if files:
        for file_path in files:
            path = Path(file_path)
            if not path.exists():
                logger.warning("[BATCH] File does not exist: %s", file_path)
                continue

            if path.suffix.lower() not in SUPPORTED_EXTENSIONS:
                logger.warning(
                    "[BATCH] Unsupported file extension: %s (file: %s)",
                    path.suffix,
                    file_path,
                )
                continue

            file_paths.append(path)

    # Remove duplicates
    file_paths = list(set(file_paths))

    logger.info("[BATCH] Total files to process: %d", len(file_paths))
    return file_paths


def ingest_single_file(
    file_path: Path,
    omics_type: Optional[str] = None,
    create_page: bool = True,
    program_ids: Optional[List[str]] = None,
    experiment_ids: Optional[List[str]] = None,
) -> Tuple[str, bool, Optional[str], Optional[str]]:
    """
    Ingest a single file.

    Args:
        file_path: Path to the file
        omics_type: Omics type (if None, will auto-detect)
        create_page: Whether to create a new Notion page
        program_ids: Optional list of program IDs to link
        experiment_ids: Optional list of experiment IDs to link

    Returns:
        Tuple of (file_path, success, page_id, error_message)
    """
    file_str = str(file_path)

    try:
        # Auto-detect omics type if not provided
        if not omics_type:
            omics_type = detect_omics_type(file_str)
            if not omics_type:
                return (file_str, False, None, "Could not auto-detect omics type")

        # Call appropriate ingestion function
        if omics_type == "lipidomics":
            page_id = ingest_lipidomics_file(
                file_path=file_str,
                notion_page_id=None,
                create_page=create_page,
                program_ids=program_ids,
                experiment_ids=experiment_ids,
            )
        elif omics_type == "metabolomics":
            page_id = ingest_metabolomics_file(
                file_path=file_str,
                notion_page_id=None,
                create_page=create_page,
                program_ids=program_ids,
                experiment_ids=experiment_ids,
            )
        elif omics_type == "proteomics":
            page_id = ingest_proteomics_file(
                file_path=file_str,
                notion_page_id=None,
                create_page=create_page,
                program_ids=program_ids,
                experiment_ids=experiment_ids,
            )
        elif omics_type == "transcriptomics":
            page_id = ingest_transcriptomics_file(
                file_path=file_str,
                notion_page_id=None,
                create_page=create_page,
                program_ids=program_ids,
                experiment_ids=experiment_ids,
            )
        else:
            return (file_str, False, None, f"Unknown omics type: {omics_type}")

        logger.info(
            "[BATCH] Successfully ingested %s (%s) -> page_id: %s",
            file_str,
            omics_type,
            page_id,
        )
        return (file_str, True, page_id, None)

    except Exception as e:
        error_msg = str(e)
        logger.error(
            "[BATCH] Failed to ingest %s (%s): %r",
            file_str,
            omics_type or "unknown",
            e,
        )
        return (file_str, False, None, error_msg)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Batch ingest multi-omics datasets with automatic type detection."
    )
    parser.add_argument(
        "--directory",
        help="Directory to scan for CSV/TSV files (recursive)",
    )
    parser.add_argument(
        "--file",
        dest="files",
        action="append",
        help="Explicit file path(s) to ingest (can be specified multiple times)",
    )
    parser.add_argument(
        "--omics-type",
        choices=["lipidomics", "metabolomics", "proteomics", "transcriptomics"],
        help="Force omics type (skip auto-detection)",
    )
    parser.add_argument(
        "--no-create-page",
        action="store_true",
        help="Do not create new Notion pages (requires --dataset-page-id, not supported in batch)",
    )
    parser.add_argument(
        "--program-id",
        dest="program_ids",
        action="append",
        default=[],
        help="Program page ID to link to all datasets (can be specified multiple times)",
    )
    parser.add_argument(
        "--experiment-id",
        dest="experiment_ids",
        action="append",
        default=[],
        help="Experiment page ID to link to all datasets (can be specified multiple times)",
    )
    parser.add_argument(
        "--parallel",
        action="store_true",
        help="Process files in parallel (use with caution for Notion API rate limits)",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=2,
        help="Maximum number of parallel workers (default: 2, recommended for Notion API)",
    )

    args = parser.parse_args()

    # Validate arguments
    if not args.directory and not args.files:
        parser.error("Either --directory or --file must be provided")

    if args.no_create_page:
        logger.warning(
            "[BATCH] --no-create-page is not supported in batch mode. "
            "All files will create new pages."
        )

    # Collect files
    try:
        file_paths = collect_files(directory=args.directory, files=args.files)
    except Exception as e:
        logger.error("[BATCH] Error collecting files: %r", e)
        print(f"\n❌ Error: {e}", file=sys.stderr)
        sys.exit(1)

    if not file_paths:
        logger.warning("[BATCH] No files found to process")
        print("\n⚠️  No files found to process")
        sys.exit(0)

    # Prepare program/experiment IDs
    program_ids = args.program_ids if args.program_ids else None
    experiment_ids = args.experiment_ids if args.experiment_ids else None

    # Process files
    start_time = time.time()
    results: List[Tuple[str, bool, Optional[str], Optional[str]]] = []

    if args.parallel:
        logger.info(
            "[BATCH] Processing %d files in parallel (max_workers=%d)",
            len(file_paths),
            args.max_workers,
        )

        with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
            futures = {
                executor.submit(
                    ingest_single_file,
                    file_path,
                    args.omics_type,
                    create_page=True,
                    program_ids=program_ids,
                    experiment_ids=experiment_ids,
                ): file_path
                for file_path in file_paths
            }

            for future in as_completed(futures):
                result = future.result()
                results.append(result)

                # Print progress
                completed = len(results)
                total = len(file_paths)
                status = "✅" if result[1] else "❌"
                print(
                    f"\r[{completed}/{total}] {status} {Path(result[0]).name}",
                    end="",
                    flush=True,
                )

    else:
        logger.info("[BATCH] Processing %d files sequentially", len(file_paths))

        for idx, file_path in enumerate(file_paths, 1):
            result = ingest_single_file(
                file_path,
                args.omics_type,
                create_page=True,
                program_ids=program_ids,
                experiment_ids=experiment_ids,
            )
            results.append(result)

            # Print progress
            status = "✅" if result[1] else "❌"
            print(
                f"\r[{idx}/{len(file_paths)}] {status} {file_path.name}",
                end="",
                flush=True,
            )

    print()  # New line after progress

    # Aggregate results
    elapsed_time = time.time() - start_time
    successful = [r for r in results if r[1]]
    failed = [r for r in results if not r[1]]

    # Group by omics type
    by_omics: Dict[str, List[Tuple[str, bool, Optional[str], Optional[str]]]] = (
        defaultdict(list)
    )
    for result in results:
        file_path = result[0]
        detected_type = detect_omics_type(file_path)
        omics_key = detected_type or "unknown"
        by_omics[omics_key].append(result)

    # Print summary
    print("\n" + "=" * 80)
    print("📊 BATCH INGESTION SUMMARY")
    print("=" * 80)
    print(f"\n⏱️  Total time: {elapsed_time:.2f} seconds")
    print(f"📁 Total files: {len(file_paths)}")
    print(f"✅ Successful: {len(successful)}")
    print(f"❌ Failed: {len(failed)}")

    print("\n📋 Results by Omics Type:")
    for omics_type, omics_results in sorted(by_omics.items()):
        omics_successful = [r for r in omics_results if r[1]]
        omics_failed = [r for r in omics_results if not r[1]]
        print(
            f"  {omics_type.upper()}: {len(omics_successful)}/{len(omics_results)} successful"
        )

    if successful:
        print("\n✅ Successfully Ingested:")
        for result in successful:
            print(f"  • {Path(result[0]).name} -> {result[2]}")

    if failed:
        print("\n❌ Failed Files:")
        for result in failed:
            print(f"  • {Path(result[0]).name}: {result[3]}")

    print("\n" + "=" * 80)

    # Exit with error code if any failures
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()

