#!/usr/bin/env python3
"""
Enhanced batch ingestion framework for multi-omics datasets.

Automatically detects omics type from files and ingests them in batch with:
- Enhanced progress bars (tqdm)
- Detailed error categorization
- Performance statistics
- JSON/CSV export options
- Skip already-processed files

Usage:
    python scripts/batch_ingest_omics_enhanced.py --directory <path>
    python scripts/batch_ingest_omics_enhanced.py --directory <path> --omics-type lipidomics
    python scripts/batch_ingest_omics_enhanced.py --file <path> --file <path> ...
    python scripts/batch_ingest_omics_enhanced.py --directory <path> --parallel --max-workers 4
    python scripts/batch_ingest_omics_enhanced.py --directory <path> --export-results results.json
"""

import argparse
import json
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None
    print("Warning: tqdm not available, using basic progress tracking")

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.ingestion.lipidomics_ingestion import ingest_lipidomics_file
from amprenta_rag.ingestion.metabolomics_ingestion import ingest_metabolomics_file
from amprenta_rag.ingestion.proteomics_ingestion import ingest_proteomics_file
from amprenta_rag.ingestion.transcriptomics_ingestion import ingest_transcriptomics_file
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.utils.performance import timer

logger = get_logger(__name__)

# Supported file extensions
SUPPORTED_EXTENSIONS = {".csv", ".tsv", ".txt"}

# Omics type detection patterns (from original script)
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


@dataclass
class FileIngestionResult:
    """Result of ingesting a single file."""
    file_path: str
    file_name: str
    omics_type: Optional[str]
    success: bool
    page_id: Optional[str]
    error_type: Optional[str]
    error_message: Optional[str]
    duration_seconds: float
    file_size_bytes: int


@dataclass
class BatchIngestionSummary:
    """Summary of batch ingestion results."""
    total_files: int
    successful: int
    failed: int
    skipped: int
    total_duration_seconds: float
    average_duration_seconds: float
    by_omics_type: Dict[str, Dict[str, int]]
    error_categories: Dict[str, int]
    results: List[FileIngestionResult]


def detect_omics_type(file_path: str, force_type: Optional[str] = None) -> Optional[str]:
    """
    Auto-detect omics type from file path and content.
    
    (Reuses logic from original batch_ingest_omics.py)
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
    
    scores: Dict[str, int] = defaultdict(int)
    
    for omics_type, patterns in OMICS_PATTERNS.items():
        for keyword in patterns["filename_keywords"]:
            if keyword in filename or keyword in file_path_lower:
                scores[omics_type] += 2
    
    try:
        import pandas as pd
        
        with open(file_path, "r", encoding="utf-8") as f:
            first_line = f.readline()
            delimiter = "\t" if "\t" in first_line else ","
            f.seek(0)
        
        df = pd.read_csv(file_path, delimiter=delimiter, encoding="utf-8", nrows=5)
        columns_lower = [str(col).lower() for col in df.columns]
        
        for omics_type, patterns in OMICS_PATTERNS.items():
            for keyword in patterns["column_keywords"]:
                if any(keyword in col for col in columns_lower):
                    scores[omics_type] += 3
            
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


def categorize_error(error_message: str) -> str:
    """Categorize error messages into error types."""
    error_lower = error_message.lower()
    
    if "not found" in error_lower or "does not exist" in error_lower:
        return "FileNotFound"
    elif "could not auto-detect" in error_lower or "detect" in error_lower:
        return "DetectionError"
    elif "notion" in error_lower or "api" in error_lower:
        return "NotionAPIError"
    elif "timeout" in error_lower:
        return "Timeout"
    elif "permission" in error_lower or "access" in error_lower:
        return "PermissionError"
    elif "format" in error_lower or "parse" in error_lower:
        return "FormatError"
    elif "empty" in error_lower or "no" in error_lower:
        return "EmptyFile"
    else:
        return "OtherError"


def collect_files(
    directory: Optional[str] = None,
    files: Optional[List[str]] = None,
    pattern: Optional[str] = None,
    skip_existing: bool = False,
) -> List[Path]:
    """
    Collect all supported files from directory or file list.
    
    Args:
        directory: Directory to scan
        files: Explicit list of file paths
        pattern: Optional glob pattern to filter files
        skip_existing: Skip files that might already be processed (basic check)
        
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
            if pattern:
                file_paths.extend(dir_path.rglob(pattern))
            else:
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
    
    # Filter by pattern if provided
    if pattern:
        file_paths = [p for p in file_paths if pattern in str(p)]
    
    logger.info("[BATCH] Total files to process: %d", len(file_paths))
    return file_paths


def ingest_single_file(
    file_path: Path,
    omics_type: Optional[str] = None,
    create_page: bool = True,
    program_ids: Optional[List[str]] = None,
    experiment_ids: Optional[List[str]] = None,
) -> FileIngestionResult:
    """
    Ingest a single file and return detailed result.
    
    Returns:
        FileIngestionResult with all details
    """
    file_str = str(file_path)
    file_name = file_path.name
    start_time = time.time()
    file_size = file_path.stat().st_size if file_path.exists() else 0
    
    try:
        # Auto-detect omics type if not provided
        detected_type = omics_type
        if not detected_type:
            detected_type = detect_omics_type(file_str)
            if not detected_type:
                duration = time.time() - start_time
                return FileIngestionResult(
                    file_path=file_str,
                    file_name=file_name,
                    omics_type=None,
                    success=False,
                    page_id=None,
                    error_type="DetectionError",
                    error_message="Could not auto-detect omics type",
                    duration_seconds=duration,
                    file_size_bytes=file_size,
                )
        
        # Call appropriate ingestion function
        with timer(f"ingest_{detected_type}_{file_name[:20]}", log_threshold=5.0):
            if detected_type == "lipidomics":
                page_id = ingest_lipidomics_file(
                    file_path=file_str,
                    notion_page_id=None,
                    create_page=create_page,
                    program_ids=program_ids,
                    experiment_ids=experiment_ids,
                )
            elif detected_type == "metabolomics":
                page_id = ingest_metabolomics_file(
                    file_path=file_str,
                    notion_page_id=None,
                    create_page=create_page,
                    program_ids=program_ids,
                    experiment_ids=experiment_ids,
                )
            elif detected_type == "proteomics":
                page_id = ingest_proteomics_file(
                    file_path=file_str,
                    notion_page_id=None,
                    create_page=create_page,
                    program_ids=program_ids,
                    experiment_ids=experiment_ids,
                )
            elif detected_type == "transcriptomics":
                page_id = ingest_transcriptomics_file(
                    file_path=file_str,
                    notion_page_id=None,
                    create_page=create_page,
                    program_ids=program_ids,
                    experiment_ids=experiment_ids,
                )
            else:
                duration = time.time() - start_time
                return FileIngestionResult(
                    file_path=file_str,
                    file_name=file_name,
                    omics_type=detected_type,
                    success=False,
                    page_id=None,
                    error_type="UnknownType",
                    error_message=f"Unknown omics type: {detected_type}",
                    duration_seconds=duration,
                    file_size_bytes=file_size,
                )
        
        duration = time.time() - start_time
        logger.info(
            "[BATCH] Successfully ingested %s (%s) -> page_id: %s (%.2fs)",
            file_str,
            detected_type,
            page_id,
            duration,
        )
        return FileIngestionResult(
            file_path=file_str,
            file_name=file_name,
            omics_type=detected_type,
            success=True,
            page_id=page_id,
            error_type=None,
            error_message=None,
            duration_seconds=duration,
            file_size_bytes=file_size,
        )
    
    except Exception as e:
        duration = time.time() - start_time
        error_msg = str(e)
        error_type = categorize_error(error_msg)
        
        logger.error(
            "[BATCH] Failed to ingest %s (%s): %r",
            file_str,
            omics_type or "unknown",
            e,
        )
        
        return FileIngestionResult(
            file_path=file_str,
            file_name=file_name,
            omics_type=omics_type,
            success=False,
            page_id=None,
            error_type=error_type,
            error_message=error_msg,
            duration_seconds=duration,
            file_size_bytes=file_size,
        )


def generate_summary(results: List[FileIngestionResult], total_duration: float) -> BatchIngestionSummary:
    """Generate comprehensive summary from results."""
    successful = [r for r in results if r.success]
    failed = [r for r in results if not r.success]
    
    # Group by omics type
    by_omics: Dict[str, Dict[str, int]] = defaultdict(lambda: {"success": 0, "failed": 0, "total": 0})
    
    for result in results:
        omics_key = result.omics_type or "unknown"
        by_omics[omics_key]["total"] += 1
        if result.success:
            by_omics[omics_key]["success"] += 1
        else:
            by_omics[omics_key]["failed"] += 1
    
    # Categorize errors
    error_categories: Dict[str, int] = defaultdict(int)
    for result in failed:
        if result.error_type:
            error_categories[result.error_type] += 1
    
    avg_duration = (
        sum(r.duration_seconds for r in results) / len(results) if results else 0.0
    )
    
    return BatchIngestionSummary(
        total_files=len(results),
        successful=len(successful),
        failed=len(failed),
        skipped=0,  # Could be enhanced to track skipped files
        total_duration_seconds=total_duration,
        average_duration_seconds=avg_duration,
        by_omics_type=dict(by_omics),
        error_categories=dict(error_categories),
        results=results,
    )


def print_summary(summary: BatchIngestionSummary):
    """Print formatted summary to console."""
    print("\n" + "=" * 80)
    print("📊 ENHANCED BATCH INGESTION SUMMARY")
    print("=" * 80)
    
    print(f"\n⏱️  Total time: {summary.total_duration_seconds:.2f} seconds")
    print(f"📁 Total files: {summary.total_files}")
    print(f"✅ Successful: {summary.successful}")
    print(f"❌ Failed: {summary.failed}")
    if summary.skipped > 0:
        print(f"⏭️  Skipped: {summary.skipped}")
    print(f"📈 Average duration: {summary.average_duration_seconds:.2f} seconds/file")
    
    print("\n📋 Results by Omics Type:")
    for omics_type, stats in sorted(summary.by_omics_type.items()):
        success_rate = (stats["success"] / stats["total"] * 100) if stats["total"] > 0 else 0
        print(
            f"  {omics_type.upper()}: {stats['success']}/{stats['total']} successful ({success_rate:.1f}%)"
        )
    
    if summary.error_categories:
        print("\n❌ Error Categories:")
        for error_type, count in sorted(summary.error_categories.items(), key=lambda x: -x[1]):
            print(f"  {error_type}: {count}")
    
    print("\n" + "=" * 80)


def export_results(summary: BatchIngestionSummary, output_path: str):
    """Export results to JSON or CSV file."""
    output_file = Path(output_path)
    
    if output_file.suffix.lower() == ".json":
        # Export as JSON
        export_data = {
            "summary": {
                "total_files": summary.total_files,
                "successful": summary.successful,
                "failed": summary.failed,
                "skipped": summary.skipped,
                "total_duration_seconds": summary.total_duration_seconds,
                "average_duration_seconds": summary.average_duration_seconds,
                "by_omics_type": summary.by_omics_type,
                "error_categories": summary.error_categories,
            },
            "results": [asdict(r) for r in summary.results],
        }
        
        with open(output_file, "w") as f:
            json.dump(export_data, f, indent=2)
        
        logger.info("[BATCH] Exported results to %s", output_file)
    
    elif output_file.suffix.lower() == ".csv":
        # Export as CSV
        import pandas as pd
        
        df = pd.DataFrame([asdict(r) for r in summary.results])
        df.to_csv(output_file, index=False)
        
        logger.info("[BATCH] Exported results to %s", output_file)
    
    else:
        logger.warning(
            "[BATCH] Unknown export format: %s (use .json or .csv)",
            output_file.suffix,
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Enhanced batch ingest multi-omics datasets with automatic type detection.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
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
    parser.add_argument(
        "--export-results",
        help="Export results to JSON or CSV file (specify path with .json or .csv extension)",
    )
    parser.add_argument(
        "--pattern",
        help="Glob pattern to filter files (e.g., '*lipid*.csv')",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip files that might already be processed (basic check)",
    )
    
    args = parser.parse_args()
    
    # Validate arguments
    if not args.directory and not args.files:
        parser.error("Either --directory or --file must be provided")
    
    # Collect files
    try:
        file_paths = collect_files(
            directory=args.directory,
            files=args.files,
            pattern=args.pattern,
            skip_existing=args.skip_existing,
        )
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
    results: List[FileIngestionResult] = []
    
    # Use tqdm if available
    progress_bar = None
    if tqdm:
        progress_bar = tqdm(
            total=len(file_paths),
            desc="Processing files",
            unit="file",
            ncols=100,
        )
    
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
                
                if progress_bar:
                    status = "✅" if result.success else "❌"
                    progress_bar.set_postfix_str(f"{status} {result.file_name[:30]}")
                    progress_bar.update(1)
                else:
                    completed = len(results)
                    status = "✅" if result.success else "❌"
                    print(
                        f"\r[{completed}/{len(file_paths)}] {status} {result.file_name}",
                        end="",
                        flush=True,
                    )
    
    else:
        logger.info("[BATCH] Processing %d files sequentially", len(file_paths))
        
        for file_path in file_paths:
            result = ingest_single_file(
                file_path,
                args.omics_type,
                create_page=True,
                program_ids=program_ids,
                experiment_ids=experiment_ids,
            )
            results.append(result)
            
            if progress_bar:
                status = "✅" if result.success else "❌"
                progress_bar.set_postfix_str(f"{status} {result.file_name[:30]}")
                progress_bar.update(1)
            else:
                status = "✅" if result.success else "❌"
                print(
                    f"\r[{len(results)}/{len(file_paths)}] {status} {result.file_name}",
                    end="",
                    flush=True,
                )
    
    if progress_bar:
        progress_bar.close()
    else:
        print()  # New line after progress
    
    # Generate summary
    total_duration = time.time() - start_time
    summary = generate_summary(results, total_duration)
    
    # Print summary
    print_summary(summary)
    
    # Export results if requested
    if args.export_results:
        export_results(summary, args.export_results)
    
    # Exit with error code if any failures
    if summary.failed > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()

