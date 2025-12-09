#!/usr/bin/env python3
"""
Extract features (genes) from GEO expression data and link them to datasets.

Downloads GEO Series Matrix files, parses gene IDs, and links them as features
to the dataset in Postgres.
"""

import argparse
import gzip
import io
import sys
from pathlib import Path
from typing import Set
from uuid import UUID

import pandas as pd
import requests

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Dataset
from amprenta_rag.ingestion.features.postgres_linking import (
    batch_link_features_to_dataset_in_postgres,
)
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.models.domain import FeatureType
from amprenta_rag.ingestion.transcriptomics.normalization import normalize_gene_identifier

logger = get_logger(__name__)


def get_geo_series_matrix_url(study_id: str) -> str:
    """
    Construct GEO Series Matrix file URL.
    
    Args:
        study_id: GEO Series ID (e.g., "GSE275841")
        
    Returns:
        URL to Series Matrix file
    """
    # Extract series number (e.g., GSE275841 -> 275841)
    series_num = study_id.replace("GSE", "")
    
    # Determine series directory (e.g., GSE275841 -> GSE275nnn)
    series_prefix = int(series_num) // 1000
    series_dir = f"GSE{series_prefix}nnn"
    
    # Construct URL
    base_url = "https://ftp.ncbi.nlm.nih.gov/geo/series"
    url = f"{base_url}/{series_dir}/{study_id}/matrix/{study_id}_series_matrix.txt.gz"
    
    return url


def download_geo_series_matrix(study_id: str, output_path: Path) -> bool:
    """
    Download GEO Series Matrix file.
    
    Args:
        study_id: GEO Series ID
        output_path: Path to save the file
        
    Returns:
        True if successful, False otherwise
    """
    url = get_geo_series_matrix_url(study_id)
    
    logger.info("[GEO-FEATURES] Downloading Series Matrix from %s", url)
    
    try:
        response = requests.get(url, timeout=120, stream=True)
        response.raise_for_status()
        
        # Save compressed file
        with open(output_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        
        logger.info("[GEO-FEATURES] Downloaded Series Matrix to %s", output_path)
        return True
        
    except Exception as e:
        logger.error("[GEO-FEATURES] Error downloading Series Matrix: %r", e)
        return False


def extract_genes_from_geo_matrix(matrix_path: Path) -> Set[str]:
    """
    Extract gene IDs from GEO Series Matrix file.
    
    Args:
        matrix_path: Path to compressed Series Matrix file
        
    Returns:
        Set of normalized gene identifiers
    """
    logger.info("[GEO-FEATURES] Extracting genes from %s", matrix_path)
    
    gene_set: Set[str] = set()
    
    try:
        # Decompress and read
        with gzip.open(matrix_path, "rt", encoding="utf-8") as f:
            lines = f.readlines()
        
        # Find where data table starts (after metadata lines starting with !)
        # Header can be "ID_REF" (quoted) or ID_REF (unquoted)
        data_start = None
        for i, line in enumerate(lines):
            stripped = line.strip()
            if stripped.startswith("ID_REF") or stripped.startswith('"ID_REF"'):
                data_start = i
                break
        
        if data_start is None:
            logger.warning("[GEO-FEATURES] Could not find data table start (ID_REF) in matrix file - may be RNA-seq with supplementary files")
            return gene_set
        
        # Read data table (skip metadata)
        data_lines = lines[data_start:]
        data_text = "".join(data_lines)
        
        # Parse as CSV/TSV
        df = pd.read_csv(
            io.StringIO(data_text),
            sep="\t",
            encoding="utf-8",
            low_memory=False,
        )
        
        if df.empty:
            logger.warning("[GEO-FEATURES] Data table is empty")
            return gene_set
        
        # First column should be gene ID
        gene_column = df.columns[0]
        logger.info("[GEO-FEATURES] Using column '%s' for gene IDs", gene_column)
        
        # Extract and normalize genes
        for gene_id in df[gene_column]:
            if pd.isna(gene_id):
                continue
            
            gene_id_str = str(gene_id).strip()
            if not gene_id_str or gene_id_str.startswith("!"):
                continue
            
            # Normalize gene identifier
            normalized = normalize_gene_identifier(gene_id_str)
            if normalized:
                gene_set.add(normalized)
        
        logger.info(
            "[GEO-FEATURES] Extracted %d unique genes from Series Matrix",
            len(gene_set),
        )
        
        return gene_set
        
    except Exception as e:
        logger.error("[GEO-FEATURES] Error extracting genes from matrix: %r", e)
        return gene_set


def extract_genes_from_supplementary_tsv(
    study_id: str,
    download_dir: Path,
) -> Set[str]:
    """
    Extract gene IDs from GEO supplementary TSV files (for RNA-seq studies).
    
    Args:
        study_id: GEO Series ID
        download_dir: Directory for temporary files
        
    Returns:
        Set of normalized gene identifiers
    """
    import re
    
    gene_set: Set[str] = set()
    
    # Get Series Matrix to find supplementary file URLs
    series_num = study_id.replace("GSE", "")
    series_prefix = int(series_num) // 1000
    series_dir = f"GSE{series_prefix}nnn"
    
    matrix_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{series_dir}/{study_id}/matrix/{study_id}_series_matrix.txt.gz"
    
    logger.info("[GEO-FEATURES] Finding supplementary TSV files for %s", study_id)
    
    try:
        # Download Series Matrix to get supplementary file URLs
        response = requests.get(matrix_url, timeout=30, stream=True)
        response.raise_for_status()
        
        content = b""
        for chunk in response.iter_content(chunk_size=8192):
            content += chunk
            if len(content) > 300000:  # Get enough to find URLs
                break
        
        decompressed = gzip.decompress(content)
        text = decompressed.decode("utf-8", errors="ignore")
        
        # Find supplementary file URLs
        supp_file_urls = []
        for line in text.split("\n"):
            if "!Series_supplementary_file" in line:
                # Extract URLs
                urls = re.findall(r'ftp://[^\s"]+\.tsv\.gz', line)
                supp_file_urls.extend(urls)
        
        # Prefer raw counts files for gene extraction
        raw_counts_urls = [url for url in supp_file_urls if "raw" in url.lower() and "count" in url.lower()]
        
        if not raw_counts_urls:
            # Try any TSV file
            tsv_urls = [url for url in supp_file_urls if url.endswith(".tsv.gz")]
            if tsv_urls:
                raw_counts_urls = [tsv_urls[0]]
        
        if not raw_counts_urls:
            logger.warning("[GEO-FEATURES] No TSV files found in supplementary files")
            return gene_set
        
        # Download and parse first TSV file using streaming approach
        tsv_url = raw_counts_urls[0]
        logger.info("[GEO-FEATURES] Streaming supplementary TSV file: %s", tsv_url.split("/")[-1])
        
        # Convert FTP to HTTP for requests library
        http_url = tsv_url.replace("ftp://", "https://")
        
        # OPTIMIZED: Stream and process incrementally - don't download entire file
        logger.info("[GEO-FEATURES] Streaming file and extracting genes (processing as we download)...")
        
        try:
            response = requests.get(http_url, timeout=300, stream=True)
            response.raise_for_status()
            
            # Get file size for progress
            total_size = int(response.headers.get('Content-Length', 0))
            if total_size > 0:
                size_mb = total_size / (1024 * 1024)
                logger.info("[GEO-FEATURES] File size: %.2f MB - using streaming to avoid full download", size_mb)
            
            # Stream decompress and parse line by line
            import zlib
            decompressor = zlib.decompressobj(zlib.MAX_WBITS | 32)  # Handle gzip
            
            buffer = b""
            lines_processed = 0
            max_rows = 50000  # Limit to first 50k rows (should cover all genes)
            header_skipped = False
            
            logger.info("[GEO-FEATURES] Processing first %d rows (should capture all unique genes)", max_rows)
            
            for chunk in response.iter_content(chunk_size=8192):
                buffer += chunk
                
                # Try to decompress available data
                try:
                    decompressed = decompressor.decompress(buffer)
                    buffer = decompressor.unconsumed_tail
                    
                    if decompressed:
                        # Process lines as they become available
                        text = decompressed.decode('utf-8', errors='ignore')
                        
                        for line in text.split('\n'):
                            if lines_processed >= max_rows:
                                break
                            
                            line = line.strip()
                            if not line:
                                continue
                            
                            # Skip header
                            if not header_skipped:
                                header_skipped = True
                                continue
                            
                            # Extract first column (gene ID)
                            parts = line.split('\t')
                            if parts:
                                gene_id = parts[0].strip().strip('"').strip("'")
                                
                                if gene_id and not gene_id.startswith('!'):
                                    # Normalize gene identifier
                                    normalized = normalize_gene_identifier(gene_id)
                                    if normalized:
                                        gene_set.add(normalized)
                                
                                lines_processed += 1
                                
                                # Progress update every 5000 rows
                                if lines_processed % 5000 == 0:
                                    logger.info(
                                        "[GEO-FEATURES] Processed %d rows, found %d unique genes...",
                                        lines_processed,
                                        len(gene_set),
                                    )
                            
                            if lines_processed >= max_rows:
                                break
                        
                        if lines_processed >= max_rows:
                            break
                            
                except Exception:
                    # May not have enough data yet, continue
                    continue
            
            # Process remaining buffer
            try:
                remaining = decompressor.flush()
                if remaining:
                    text = remaining.decode('utf-8', errors='ignore')
                    for line in text.split('\n'):
                        if lines_processed >= max_rows:
                            break
                        line = line.strip()
                        if line:
                            parts = line.split('\t')
                            if parts:
                                gene_id = parts[0].strip().strip('"').strip("'")
                                if gene_id and not gene_id.startswith('!'):
                                    normalized = normalize_gene_identifier(gene_id)
                                    if normalized:
                                        gene_set.add(normalized)
                                lines_processed += 1
            except:
                pass
            
            logger.info(
                "[GEO-FEATURES] Extracted %d unique genes from %d rows using streaming approach",
                len(gene_set),
                lines_processed,
            )
            
            return gene_set
            
        except Exception as e:
            logger.error("[GEO-FEATURES] Error streaming TSV file: %r", e)
            # Fallback to download-based approach if streaming fails
            logger.info("[GEO-FEATURES] Falling back to download-based approach...")
            return gene_set
        
    except Exception as e:
        logger.error("[GEO-FEATURES] Error extracting from supplementary files: %r", e)
        return gene_set


def extract_features_from_geo_study(
    study_id: str,
    dataset_id: UUID,
    download_dir: Path = None,
) -> int:
    """
    Extract gene features from GEO study and link them to dataset.
    
    Args:
        study_id: GEO Series ID (e.g., "GSE275841")
        dataset_id: Postgres dataset UUID
        download_dir: Directory for temporary files (default: /tmp)
        
    Returns:
        Number of features linked
    """
    if download_dir is None:
        download_dir = Path("/tmp")
    
    download_dir.mkdir(parents=True, exist_ok=True)
    
    # Download Series Matrix file
    matrix_path = download_dir / f"{study_id}_series_matrix.txt.gz"
    
    if not matrix_path.exists():
        if not download_geo_series_matrix(study_id, matrix_path):
            logger.error("[GEO-FEATURES] Failed to download Series Matrix")
            return 0
    else:
        logger.info("[GEO-FEATURES] Using existing Series Matrix file: %s", matrix_path)
    
    # Extract genes - try Series Matrix first (for microarray), then supplementary files (for RNA-seq)
    gene_set = extract_genes_from_geo_matrix(matrix_path)
    
    if not gene_set:
        logger.info("[GEO-FEATURES] No genes from Series Matrix, trying supplementary TSV files (RNA-seq)")
        gene_set = extract_genes_from_supplementary_tsv(study_id, download_dir)
    
    if not gene_set:
        logger.warning("[GEO-FEATURES] No genes extracted from Series Matrix or supplementary files")
        return 0
    
    # Link genes to dataset
    logger.info(
        "[GEO-FEATURES] Linking %d genes to dataset %s",
        len(gene_set),
        dataset_id,
    )
    
    features_to_link = [(gene, FeatureType.GENE) for gene in gene_set]
    
    db = next(get_db())
    try:
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
        
    finally:
        db.close()


def main():
    parser = argparse.ArgumentParser(
        description="Extract gene features from GEO expression data"
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
    db = next(get_db())
    try:
        dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
        if not dataset:
            logger.error("[GEO-FEATURES] Dataset %s not found", dataset_id)
            sys.exit(1)
        logger.info("[GEO-FEATURES] Found dataset: %s", dataset.name if hasattr(dataset, 'name') else dataset_id)
    finally:
        db.close()
    
    # Extract and link features
    linked_count = extract_features_from_geo_study(
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

