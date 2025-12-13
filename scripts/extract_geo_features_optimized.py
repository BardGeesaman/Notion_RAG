#!/usr/bin/env python3
"""
Optimized feature extraction for GEO - streams and processes efficiently.

Downloads only what's needed and processes incrementally to handle large files.
"""

import argparse
import gzip
import re
import sys
from pathlib import Path
from typing import Set
from uuid import UUID

import requests

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import Dataset
from amprenta_rag.ingestion.features.postgres_linking import (
    batch_link_features_to_dataset_in_postgres,
)
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.models.domain import FeatureType
from amprenta_rag.ingestion.transcriptomics.normalization import normalize_gene_identifier

logger = get_logger(__name__)


def extract_genes_from_tsv_stream(
    tsv_url: str,
    max_rows: int = 50000,
    download_dir: Path = None,
) -> Set[str]:
    """
    Stream and extract genes from TSV file efficiently.
    
    Only processes the first column and limits rows to avoid full file download.
    
    Args:
        tsv_url: URL to TSV file
        max_rows: Maximum rows to process (default: 50k - usually enough for all genes)
        download_dir: Directory for temporary files
        
    Returns:
        Set of normalized gene identifiers
    """
    gene_set: Set[str] = set()
    
    # Convert FTP to HTTP
    http_url = tsv_url.replace("ftp://", "https://")
    filename = tsv_url.split("/")[-1]
    
    logger.info("[GEO-FEATURES] Streaming TSV file: %s", filename)
    logger.info("[GEO-FEATURES] Processing first %d rows (should capture all unique genes)", max_rows)
    
    try:
        # Stream download with progress
        response = requests.get(http_url, timeout=300, stream=True)
        response.raise_for_status()
        
        # Check file size
        total_size = int(response.headers.get('Content-Length', 0))
        if total_size > 0:
            size_mb = total_size / (1024 * 1024)
            logger.info("[GEO-FEATURES] File size: %.2f MB", size_mb)
        
        # Stream and decompress incrementally
        decompressor = gzip.decompressobj()
        buffer = b""
        lines_read = 0
        bytes_read = 0
        
        logger.info("[GEO-FEATURES] Streaming and parsing...")
        
        for chunk in response.iter_content(chunk_size=8192):
            bytes_read += len(chunk)
            buffer += chunk
            
            # Decompress available data
            try:
                decompressed = decompressor.decompress(buffer)
                buffer = decompressor.unconsumed_tail
                
                if decompressed:
                    # Process decompressed text
                    text = decompressed.decode('utf-8', errors='ignore')
                    lines = text.split('\n')
                    
                    # Process lines
                    for line in lines:
                        if lines_read >= max_rows:
                            break
                        
                        if not line.strip():
                            continue
                        
                        # First line is header - skip or use to identify gene column
                        if lines_read == 0:
                            lines_read += 1
                            continue
                        
                        # Extract first column (gene ID)
                        parts = line.split('\t')
                        if parts:
                            gene_id = parts[0].strip().strip('"')
                            
                            if gene_id and not gene_id.startswith('!'):
                                normalized = normalize_gene_identifier(gene_id)
                                if normalized:
                                    gene_set.add(normalized)
                        
                        lines_read += 1
                        
                        # Progress update every 1000 rows
                        if lines_read % 1000 == 0:
                            logger.debug("[GEO-FEATURES] Processed %d rows, found %d unique genes", 
                                       lines_read, len(gene_set))
                    
                    if lines_read >= max_rows:
                        break
                        
            except Exception:
                # If decompression fails, we might not have enough data yet
                continue
        
        # Process remaining buffer
        if buffer and lines_read < max_rows:
            try:
                remaining = decompressor.flush()
                if remaining:
                    text = remaining.decode('utf-8', errors='ignore')
                    for line in text.split('\n'):
                        if lines_read >= max_rows:
                            break
                        if line.strip():
                            parts = line.split('\t')
                            if parts:
                                gene_id = parts[0].strip().strip('"')
                                if gene_id and not gene_id.startswith('!'):
                                    normalized = normalize_gene_identifier(gene_id)
                                    if normalized:
                                        gene_set.add(normalized)
                            lines_read += 1
            except:
                pass
        
        logger.info(
            "[GEO-FEATURES] Processed %d rows, extracted %d unique genes",
            lines_read,
            len(gene_set),
        )
        
        return gene_set
        
    except Exception as e:
        logger.error("[GEO-FEATURES] Error streaming TSV file: %r", e)
        return gene_set


def extract_genes_from_supplementary_files_optimized(
    study_id: str,
    download_dir: Path = None,
) -> Set[str]:
    """
    Optimized extraction from GEO supplementary TSV files.
    
    Uses streaming to avoid downloading entire large files.
    
    Args:
        study_id: GEO Series ID
        download_dir: Directory for temporary files (not used in streaming mode)
        
    Returns:
        Set of normalized gene identifiers
    """
    gene_set: Set[str] = set()
    
    # Get Series Matrix to find supplementary file URLs
    series_num = study_id.replace("GSE", "")
    series_prefix = int(series_num) // 1000
    series_dir = f"GSE{series_prefix}nnn"
    
    matrix_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{series_dir}/{study_id}/matrix/{study_id}_series_matrix.txt.gz"
    
    logger.info("[GEO-FEATURES] Finding supplementary TSV files for %s", study_id)
    
    try:
        # Download Series Matrix header to get supplementary file URLs
        response = requests.get(matrix_url, timeout=30, stream=True)
        response.raise_for_status()
        
        content = b""
        for chunk in response.iter_content(chunk_size=8192):
            content += chunk
            if len(content) > 300000:  # Enough to find URLs
                break
        
        decompressed = gzip.decompress(content)
        text = decompressed.decode("utf-8", errors="ignore")
        
        # Find supplementary file URLs
        supp_file_urls = []
        for line in text.split("\n"):
            if "!Series_supplementary_file" in line:
                urls = re.findall(r'ftp://[^\s"]+\.tsv\.gz', line)
                supp_file_urls.extend(urls)
        
        if not supp_file_urls:
            logger.warning("[GEO-FEATURES] No TSV files found in supplementary files")
            return gene_set
        
        # Prefer raw counts files
        raw_counts_urls = [url for url in supp_file_urls if "raw" in url.lower() and "count" in url.lower()]
        if not raw_counts_urls:
            tsv_urls = [url for url in supp_file_urls if url.endswith(".tsv.gz")]
            if tsv_urls:
                raw_counts_urls = [tsv_urls[0]]
        
        if not raw_counts_urls:
            logger.warning("[GEO-FEATURES] No suitable TSV files found")
            return gene_set
        
        # Stream and extract from first TSV file
        tsv_url = raw_counts_urls[0]
        logger.info("[GEO-FEATURES] Using file: %s", tsv_url.split("/")[-1])
        
        gene_set = extract_genes_from_tsv_stream(tsv_url, max_rows=50000)
        
        return gene_set
        
    except Exception as e:
        logger.error("[GEO-FEATURES] Error finding/processing supplementary files: %r", e)
        return gene_set


def extract_features_from_geo_study_optimized(
    study_id: str,
    dataset_id: UUID,
    download_dir: Path = None,
) -> int:
    """
    Optimized gene feature extraction from GEO study.
    
    Uses streaming and incremental processing to handle large files efficiently.
    
    Args:
        study_id: GEO Series ID
        dataset_id: Postgres dataset UUID
        download_dir: Directory for temporary files
        
    Returns:
        Number of features linked
    """
    logger.info("[GEO-FEATURES] Starting optimized feature extraction for %s", study_id)
    
    # Try supplementary files first (RNA-seq), then Series Matrix (microarray)
    gene_set = extract_genes_from_supplementary_files_optimized(study_id)
    
    if not gene_set:
        logger.info("[GEO-FEATURES] No genes from supplementary files, trying Series Matrix...")
        # Fall back to Series Matrix if needed (smaller files, usually works)
        from scripts.extract_geo_features import extract_genes_from_geo_matrix, download_geo_series_matrix
        
        if download_dir is None:
            download_dir = Path("/tmp")
        download_dir.mkdir(parents=True, exist_ok=True)
        
        matrix_path = download_dir / f"{study_id}_series_matrix.txt.gz"
        if not matrix_path.exists():
            download_geo_series_matrix(study_id, matrix_path)
        
        if matrix_path.exists():
            gene_set = extract_genes_from_geo_matrix(matrix_path)
    
    if not gene_set:
        logger.warning("[GEO-FEATURES] No genes extracted from any source")
        return 0
    
    # Link genes to dataset
    logger.info("[GEO-FEATURES] Linking %d genes to dataset %s", len(gene_set), dataset_id)
    
    features_to_link = [(gene, FeatureType.GENE) for gene in gene_set]
    
    with db_session() as db:
        results = batch_link_features_to_dataset_in_postgres(
            features=features_to_link,
            dataset_id=dataset_id,
            db=db,
        )
        
        linked_count = sum(1 for v in results.values() if v)
        logger.info(
            "[GEO-FEATURES] Successfully linked %d/%d genes to dataset %s",
            linked_count,
            len(gene_set),
            dataset_id,
        )
        
        return linked_count


def main():
    parser = argparse.ArgumentParser(
        description="Extract gene features from GEO expression data (optimized)"
    )
    parser.add_argument(
        "--study-id",
        required=True,
        help="GEO Series ID (e.g., GSE275841)",
    )
    parser.add_argument(
        "--dataset-id",
        required=True,
        help="Postgres dataset UUID to link features to",
    )
    parser.add_argument(
        "--download-dir",
        type=Path,
        default=Path("/tmp"),
        help="Directory for temporary files (default: /tmp)",
    )
    
    args = parser.parse_args()
    
    try:
        dataset_id = UUID(args.dataset_id)
    except ValueError:
        logger.error("[GEO-FEATURES] Invalid dataset ID: %s", args.dataset_id)
        sys.exit(1)
    
    # Verify dataset exists
    with db_session() as db:
        dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
        if not dataset:
            logger.error("[GEO-FEATURES] Dataset %s not found", dataset_id)
            sys.exit(1)
        logger.info("[GEO-FEATURES] Found dataset: %s", dataset.name if hasattr(dataset, 'name') else dataset_id)
    
    # Extract and link features
    linked_count = extract_features_from_geo_study_optimized(
        study_id=args.study_id,
        dataset_id=dataset_id,
        download_dir=args.download_dir,
    )
    
    if linked_count > 0:
        logger.info("✅ Successfully linked %d gene features to dataset", linked_count)
        sys.exit(0)
    else:
        logger.error("❌ No features were linked")
        sys.exit(1)


if __name__ == "__main__":
    main()

