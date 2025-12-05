"""
Unified feature extraction for all repository types.

Extracts features (genes, proteins, metabolites, lipids) from repository datasets
and links them to Postgres datasets.

Supports:
- GEO: Genes from expression matrices (using GEOparse library)
- PRIDE: Proteins from identification files
- MetaboLights: Metabolites from ISA-Tab data files
- MW: Metabolites from REST API /data endpoint
"""

from __future__ import annotations

import io
import os
import re
import time
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from uuid import UUID

import GEOparse
import pandas as pd
import requests
from Bio import Entrez

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Dataset
from amprenta_rag.ingestion.features.postgres_linking import (
    batch_link_features_to_dataset_in_postgres,
)
from amprenta_rag.ingestion.transcriptomics.normalization import normalize_gene_identifier
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.models.domain import FeatureType

logger = get_logger(__name__)


def extract_geo_features_with_geoparse(
    study_id: str,
    download_dir: Path = None,
) -> Set[str]:
    """
    Extract gene features from GEO study using GEOparse library.
    
    This is the recommended approach - cleaner and more robust than custom parsing.
    
    Args:
        study_id: GEO Series ID (e.g., "GSE12251")
        download_dir: Directory for caching downloaded files
        
    Returns:
        Set of normalized gene identifiers
    """
    if download_dir is None:
        download_dir = Path("/tmp/geo_cache")
    
    download_dir.mkdir(parents=True, exist_ok=True)
    
    gene_set: Set[str] = set()
    
    logger.info("[FEATURE-EXTRACT][GEO] Extracting genes from %s using GEOparse", study_id)
    
    try:
        # Configure Entrez (for GEOparse's internal use)
        from amprenta_rag.config import get_config
        cfg = get_config()
        
        email = getattr(cfg, "ncbi_email", None) or os.getenv("NCBI_EMAIL", "")
        api_key = getattr(cfg, "geo_api_key", None) or os.getenv("GEO_API_KEY", "")
        
        if email:
            Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        
        # Download and parse using GEOparse
        # destdir caches files locally so they aren't re-downloaded
        gse = GEOparse.get_GEO(geo=study_id, destdir=str(download_dir), silent=True)
        
        logger.info("[FEATURE-EXTRACT][GEO] Successfully parsed %s: %s", study_id, gse.name)
        
        # Extract genes from expression matrix
        # GEOparse provides pivot_samples('VALUE') which gives expression matrix
        # Rows = probes/genes, Columns = samples
        try:
            expression_matrix = gse.pivot_samples('VALUE')
            
            if expression_matrix is not None and not expression_matrix.empty:
                logger.info(
                    "[FEATURE-EXTRACT][GEO] Expression matrix: %d probes/genes x %d samples",
                    len(expression_matrix),
                    len(expression_matrix.columns),
                )
                
                # Extract gene IDs from row index
                for gene_id in expression_matrix.index:
                    if pd.isna(gene_id):
                        continue
                    
                    gene_id_str = str(gene_id).strip()
                    if not gene_id_str:
                        continue
                    
                    # Normalize gene identifier
                    normalized = normalize_gene_identifier(gene_id_str)
                    if normalized:
                        gene_set.add(normalized)
                
                logger.info(
                    "[FEATURE-EXTRACT][GEO] Extracted %d unique genes from expression matrix",
                    len(gene_set),
                )
            else:
                logger.warning("[FEATURE-EXTRACT][GEO] Expression matrix is empty for %s", study_id)
                
                # Try to get genes from platform annotation if available
                if hasattr(gse, 'gpls') and gse.gpls:
                    for gpl_id, gpl in gse.gpls.items():
                        if hasattr(gpl, 'table') and gpl.table is not None:
                            # Platform table might have gene annotations
                            logger.info("[FEATURE-EXTRACT][GEO] Checking platform %s for gene annotations", gpl_id)
                            # Could extract genes from platform table if needed
                
        except Exception as e:
            logger.warning("[FEATURE-EXTRACT][GEO] Could not extract expression matrix: %r", e)
            
            # Fallback: Try to extract from platform data
            if hasattr(gse, 'gpls') and gse.gpls:
                logger.info("[FEATURE-EXTRACT][GEO] Attempting fallback extraction from platform data")
                for gpl_id, gpl in gse.gpls.items():
                    if hasattr(gpl, 'table') and gpl.table is not None:
                        # Try to find gene column in platform table
                        table = gpl.table
                        gene_columns = ['GENE', 'Gene Symbol', 'GENE_SYMBOL', 'GENE_NAME', 'ID']
                        
                        for col in gene_columns:
                            if col in table.columns:
                                for gene_id in table[col]:
                                    if pd.isna(gene_id):
                                        continue
                                    gene_id_str = str(gene_id).strip()
                                    if gene_id_str:
                                        normalized = normalize_gene_identifier(gene_id_str)
                                        if normalized:
                                            gene_set.add(normalized)
                                break
                        
                        if gene_set:
                            logger.info(
                                "[FEATURE-EXTRACT][GEO] Extracted %d genes from platform annotation",
                                len(gene_set),
                            )
                            break
        
        return gene_set
        
    except Exception as e:
        logger.error("[FEATURE-EXTRACT][GEO] Error extracting genes with GEOparse: %r", e)
        return gene_set


def extract_pride_proteins_from_data_files(
    study_id: str,
    download_dir: Path = None,
) -> Set[str]:
    """
    Extract protein features from PRIDE data files using optimized priority-based approach.
    
    Priority order:
    1. *.mzTab files (Community standard)
    2. protein_groups.txt (MaxQuant standard)
    3. *.xls/*.xlsx (Supplied spreadsheets)
    4. Other TSV/CSV result files (fallback)
    
    Uses pandas for efficient data extraction.
    
    Args:
        study_id: PRIDE project ID (e.g., "PXD012345")
        download_dir: Directory for temporary files (optional, for caching)
        
    Returns:
        Set of normalized protein identifiers (accessions)
    """
    if download_dir is None:
        download_dir = Path("/tmp")
    
    download_dir.mkdir(parents=True, exist_ok=True)
    
    protein_set: Set[str] = set()
    
    # Ensure study_id starts with PXD
    if not study_id.startswith("PXD"):
        study_id = f"PXD{study_id}"
    
    logger.info("[FEATURE-EXTRACT][PRIDE] Finding best protein data file for %s", study_id)
    
    try:
        # Import here to avoid circular imports
        from amprenta_rag.ingestion.repositories.pride import PRIDERepository
        
        # Get repository instance
        repo = PRIDERepository()
        
        # Fetch ALL available data files (no filtering yet)
        data_files = repo.fetch_study_data_files(study_id=study_id)
        
        if not data_files:
            logger.warning("[FEATURE-EXTRACT][PRIDE] No data files found for %s", study_id)
            return protein_set
        
        logger.info("[FEATURE-EXTRACT][PRIDE] Found %d files, applying priority selection...", len(data_files))
        
        # Priority-based file selection
        target_file = None
        
        # Priority 1: Look for *.mzTab files (Community standard)
        mztab_files = [
            f for f in data_files
            if f.filename.lower().endswith(".mztab")
        ]
        if mztab_files:
            target_file = mztab_files[0]
            logger.info("[FEATURE-EXTRACT][PRIDE] Priority 1: Selected mzTab file: %s", target_file.filename)
        
        # Priority 2: Look for protein_groups.txt (MaxQuant standard)
        if not target_file:
            maxquant_files = [
                f for f in data_files
                if "protein_groups.txt" in f.filename.lower()
            ]
            if maxquant_files:
                target_file = maxquant_files[0]
                logger.info("[FEATURE-EXTRACT][PRIDE] Priority 2: Selected MaxQuant file: %s", target_file.filename)
        
        # Priority 3: Look for Excel files (*.xls, *.xlsx)
        if not target_file:
            excel_files = [
                f for f in data_files
                if f.filename.lower().endswith((".xls", ".xlsx"))
            ]
            if excel_files:
                target_file = excel_files[0]
                logger.info("[FEATURE-EXTRACT][PRIDE] Priority 3: Selected Excel file: %s", target_file.filename)
        
        # Fallback: Look for TSV/CSV result files
        if not target_file:
            result_files = [
                f for f in data_files
                if (f.filename.lower().endswith((".tsv", ".csv", ".txt")) and
                    (f.file_type in ["TSV", "CSV", "TXT", "SEARCH"] or
                     any(kw in f.filename.lower() for kw in ["protein", "result", "peptide", "identification"])))
            ]
            if result_files:
                target_file = result_files[0]
                logger.info("[FEATURE-EXTRACT][PRIDE] Fallback: Selected result file: %s", target_file.filename)
        
        if not target_file:
            logger.warning("[FEATURE-EXTRACT][PRIDE] No suitable result file found for %s", study_id)
            return protein_set
        
        # Convert FTP URL to HTTP for requests library
        download_url = target_file.download_url
        if download_url.startswith("ftp://"):
            download_url = download_url.replace("ftp://ftp.pride.ebi.ac.uk", "https://ftp.pride.ebi.ac.uk")
            logger.info("[FEATURE-EXTRACT][PRIDE] Converted FTP to HTTP")
        
        # Download and parse based on file type
        logger.info("[FEATURE-EXTRACT][PRIDE] Downloading and parsing %s...", target_file.filename)
        
        time.sleep(1.0)  # Rate limiting
        from amprenta_rag.ingestion.repositories import REPOSITORY_USER_AGENT
        headers = {"User-Agent": REPOSITORY_USER_AGENT}
        response = requests.get(download_url, headers=headers, timeout=300)
        response.raise_for_status()
        
        df = None
        
        # Parse based on file type using pandas
        if target_file.filename.lower().endswith('.mztab'):
            # Parse mzTab format
            # mzTab format: Header lines start with 'PRH', data lines start with 'PRT'
            content = response.text
            lines = content.split('\n')
            
            # Find PRH (Protein Header) line
            header_line = None
            for line in lines:
                if line.startswith('PRH'):
                    header_line = line
                    break
            
            if header_line:
                # Extract column names (skip 'PRH\t' prefix)
                cols = header_line.replace('PRH\t', '').split('\t')
                
                # Extract PRT (Protein) data lines
                data_lines = [line.replace('PRT\t', '') for line in lines if line.startswith('PRT')]
                
                if data_lines:
                    df = pd.read_csv(
                        io.StringIO('\n'.join(data_lines)),
                        sep='\t',
                        names=cols,
                        low_memory=False,
                    )
            else:
                logger.warning("[FEATURE-EXTRACT][PRIDE] No PRH header found in mzTab file")
                
        elif 'protein_groups.txt' in target_file.filename.lower():
            # Parse MaxQuant protein_groups.txt (tab-separated)
            df = pd.read_csv(
                io.StringIO(response.text),
                sep='\t',
                low_memory=False,
            )
            
        elif target_file.filename.lower().endswith(('.xls', '.xlsx')):
            # Parse Excel file
            # Note: openpyxl needed for .xlsx files - will fail gracefully if not installed
            try:
                engine = 'openpyxl' if target_file.filename.lower().endswith('.xlsx') else None
                df = pd.read_excel(
                    io.BytesIO(response.content),
                    engine=engine,
                )
            except ImportError:
                logger.warning(
                    "[FEATURE-EXTRACT][PRIDE] openpyxl not installed, cannot parse Excel file. "
                    "Install with: pip install openpyxl"
                )
                return protein_set
            
        else:
            # Fallback: Parse as TSV/CSV
            separator = '\t' if target_file.filename.lower().endswith('.tsv') else ','
            df = pd.read_csv(
                io.StringIO(response.text),
                sep=separator,
                low_memory=False,
            )
        
        if df is None or df.empty:
            logger.warning("[FEATURE-EXTRACT][PRIDE] Could not parse file or file is empty")
            return protein_set
        
        logger.info("[FEATURE-EXTRACT][PRIDE] Parsed DataFrame: %d rows x %d columns", len(df), len(df.columns))
        
        # Find protein accession column
        accession_col = None
        for col in df.columns:
            col_lower = str(col).lower().strip()
            if any(keyword in col_lower for keyword in [
                "accession", "accessions", "protein", "protein id",
                "proteinid", "uniprot", "uniprot id"
            ]):
                accession_col = col
                break
        
        if not accession_col:
            # Default to first column
            accession_col = df.columns[0]
            logger.warning("[FEATURE-EXTRACT][PRIDE] No accession column found, using first column: %s", accession_col)
        
        # Extract and normalize protein accessions
        from amprenta_rag.ingestion.proteomics.normalization import normalize_protein_identifier
        
        for protein_id in df[accession_col]:
            if pd.isna(protein_id):
                continue
            
            protein_id_str = str(protein_id).strip().strip('"').strip("'")
            
            if protein_id_str and not protein_id_str.startswith('#'):
                # Normalize protein identifier
                normalized = normalize_protein_identifier(protein_id_str)
                if normalized:
                    protein_set.add(normalized)
        
        logger.info(
            "[FEATURE-EXTRACT][PRIDE] Extracted %d unique protein accessions from %s",
            len(protein_set),
            target_file.filename,
        )
        
        return protein_set
        
    except Exception as e:
        logger.error("[FEATURE-EXTRACT][PRIDE] Error extracting proteins: %r", e)
        import traceback
        logger.debug("[FEATURE-EXTRACT][PRIDE] Traceback: %s", traceback.format_exc())
        return protein_set


def extract_metabolights_metabolites_from_isa_tab(
    study_id: str,
    download_dir: Path = None,
) -> Set[str]:
    """
    Extract metabolite features from MetaboLights ISA-Tab files.
    
    Uses the study details to get HTTP URL and finds metabolite data files (m_*.tsv).
    Parses ISA-Tab format to extract metabolite identifiers.
    
    Args:
        study_id: MetaboLights study ID (e.g., "MTBLS1")
        download_dir: Directory for temporary files
        
    Returns:
        Set of normalized metabolite identifiers
    """
    if download_dir is None:
        download_dir = Path("/tmp")
    
    download_dir.mkdir(parents=True, exist_ok=True)
    
    metabolite_set: Set[str] = set()
    
    # Ensure study_id format
    study_id = study_id.upper()
    if not study_id.startswith("MTBLS"):
        logger.error("[FEATURE-EXTRACT][METABOLIGHTS] Invalid study ID format: %s", study_id)
        return metabolite_set
    
    logger.info("[FEATURE-EXTRACT][METABOLIGHTS] Extracting metabolites from %s", study_id)
    
    try:
        # Import here to avoid circular imports
        from amprenta_rag.ingestion.repositories.metabolights import MetaboLightsRepository
        
        # Get repository instance
        repo = MetaboLightsRepository()
        
        # Fetch study metadata to get HTTP URL
        metadata = repo.fetch_study_metadata(study_id)
        if not metadata or not metadata.raw_metadata:
            logger.warning("[FEATURE-EXTRACT][METABOLIGHTS] Could not fetch metadata for %s", study_id)
            return metabolite_set
        
        # Get HTTP URL from study details
        mtbls_study = metadata.raw_metadata.get("mtblsStudy", {})
        http_url = mtbls_study.get("studyHttpUrl", "")
        
        if not http_url:
            logger.warning("[FEATURE-EXTRACT][METABOLIGHTS] No HTTP URL found for %s", study_id)
            return metabolite_set
        
        # Convert FTP to HTTPS
        if http_url.startswith("ftp://"):
            https_url = http_url.replace("ftp://", "https://")
        elif http_url.startswith("http://ftp."):
            https_url = http_url.replace("http://ftp.", "https://ftp.")
        else:
            https_url = http_url
        
        logger.info("[FEATURE-EXTRACT][METABOLIGHTS] Study directory: %s", https_url)
        
        # Find MAF (Metabolite Assignment File) files (m_*_maf.tsv)
        # MAF files are the core metabolite data files containing metabolite identification
        # in metadata columns (left side) followed by sample abundance columns (right side)
        metabolite_file_url = None
        
        # Check investigation file to find actual file names
        inv_url = f"{https_url}/i_Investigation.txt"
        try:
            from amprenta_rag.ingestion.repositories import REPOSITORY_USER_AGENT
            headers = {"User-Agent": REPOSITORY_USER_AGENT}
            inv_response = requests.get(inv_url, headers=headers, timeout=30)
            if inv_response.status_code == 200:
                inv_content = inv_response.text
                
                # Look for MAF files first (m_*_maf.tsv) - these are the metabolite data files
                maf_files = re.findall(r'(m_[^\s\t\n]*maf[^\s\t\n]*\.(?:tsv|txt))', inv_content, re.IGNORECASE)
                if maf_files:
                    metabolite_file_url = f"{https_url}/{maf_files[0]}"
                    logger.info("[FEATURE-EXTRACT][METABOLIGHTS] Found MAF file from investigation: %s", maf_files[0])
                
                # Fall back to other m_ files (metabolite data files)
                if not metabolite_file_url:
                    m_files = re.findall(r'(m_[^\s\t\n]+\.(?:tsv|txt))', inv_content)
                    if m_files:
                        metabolite_file_url = f"{https_url}/{m_files[0]}"
                        logger.info("[FEATURE-EXTRACT][METABOLIGHTS] Found metabolite data file from investigation: %s", m_files[0])
                
                # Last resort: assay files (a_*.txt) - these contain sample metadata, not metabolite data
                if not metabolite_file_url:
                    a_files = re.findall(r'(a_[^\s\t\n]+metabolite[^\s\t\n]+\.txt)', inv_content, re.IGNORECASE)
                    if a_files:
                        metabolite_file_url = f"{https_url}/{a_files[0]}"
                        logger.warning("[FEATURE-EXTRACT][METABOLIGHTS] Using assay file (may contain sample metadata, not metabolite data): %s", a_files[0])
        except Exception as e:
            logger.debug("[FEATURE-EXTRACT][METABOLIGHTS] Could not read investigation file: %r", e)
        
        # If not found, try common patterns
        if not metabolite_file_url:
            metabolite_file_patterns = [
                f"m_{study_id}_LC-MS_positive_maf.tsv",
                f"m_{study_id}_LC-MS_negative_maf.tsv",
                f"m_{study_id}_GC-MS_maf.tsv",
                f"m_{study_id}_NMR_maf.tsv",
                f"m_{study_id}_metabolite_profiling_mass_spectrometry.tsv",
                f"m_{study_id}_metabolite_profiling_NMR_spectroscopy.tsv",
            ]
            
            for pattern in metabolite_file_patterns:
                test_url = f"{https_url}/{pattern}"
                test_response = requests.head(test_url, timeout=10)
                if test_response.status_code == 200:
                    metabolite_file_url = test_url
                    logger.info("[FEATURE-EXTRACT][METABOLIGHTS] Found metabolite file: %s", pattern)
                    break
        
        if not metabolite_file_url:
            logger.warning("[FEATURE-EXTRACT][METABOLIGHTS] No metabolite data file found for %s", study_id)
            return metabolite_set
        
        # Download and parse metabolite file
        logger.info("[FEATURE-EXTRACT][METABOLIGHTS] Downloading metabolite file: %s", metabolite_file_url)
        
        from amprenta_rag.ingestion.repositories import REPOSITORY_USER_AGENT
        headers = {"User-Agent": REPOSITORY_USER_AGENT}
        response = requests.get(metabolite_file_url, headers=headers, timeout=120, stream=True)
        response.raise_for_status()
        
        # Parse ISA-Tab metabolite file
        # First column is usually metabolite identifier/name
        is_gzipped = metabolite_file_url.endswith('.gz') or metabolite_file_url.endswith('.gzip')
        
        import zlib
        if is_gzipped:
            decompressor = zlib.decompressobj(zlib.MAX_WBITS | 32)
        
        buffer = b""
        lines_processed = 0
        max_rows = 50000
        header_skipped = False
        metabolite_column_idx = None
        
        logger.info("[FEATURE-EXTRACT][METABOLIGHTS] Streaming file and extracting metabolites...")
        
        for chunk in response.iter_content(chunk_size=8192):
            buffer += chunk
            
            if is_gzipped:
                try:
                    decompressed = decompressor.decompress(buffer)
                    buffer = decompressor.unconsumed_tail
                except:
                    continue
            else:
                decompressed = buffer
                buffer = b""
            
            if decompressed:
                text = decompressed.decode('utf-8', errors='ignore')
                
                for line in text.split('\n'):
                    if lines_processed >= max_rows:
                        break
                    
                    line = line.strip()
                    if not line:
                        continue
                    
                    # ISA-Tab uses tab separator
                    parts = line.split('\t')
                    
                    if not header_skipped:
                        header_skipped = True
                        # For MAF files, look for metabolite identification column
                        # Priority: metabolite_identification > database_identifier > chemical_formula
                        # Sample columns come after metadata columns, so we stop when we see sample names
                        for i, col in enumerate(parts):
                            col_lower = col.lower().strip()
                            
                            # Check if this is a sample column (stops metadata section)
                            if col_lower in ["sample name", "sample_name"] or col_lower.startswith("sample_"):
                                # We've reached sample columns, stop looking
                                break
                            
                            # Look for metabolite identification column (highest priority)
                            if any(keyword in col_lower for keyword in [
                                "metabolite_identification", "metabolite identification",
                                "metabolite identifier", "metabolite_identifier"
                            ]):
                                metabolite_column_idx = i
                                logger.info("[FEATURE-EXTRACT][METABOLIGHTS] Using column '%s' for metabolite identification", col)
                                break
                        
                        # If not found, try database_identifier
                        if metabolite_column_idx is None:
                            for i, col in enumerate(parts):
                                col_lower = col.lower().strip()
                                if col_lower in ["database_identifier", "database identifier"]:
                                    metabolite_column_idx = i
                                    logger.info("[FEATURE-EXTRACT][METABOLIGHTS] Using column '%s' for metabolite identification", col)
                                    break
                        
                        # If still not found, try chemical_formula or compound name
                        if metabolite_column_idx is None:
                            for i, col in enumerate(parts):
                                col_lower = col.lower().strip()
                                if any(keyword in col_lower for keyword in [
                                    "chemical_formula", "compound", "metabolite_name",
                                    "name", "identifier"
                                ]):
                                    metabolite_column_idx = i
                                    logger.info("[FEATURE-EXTRACT][METABOLIGHTS] Using column '%s' for metabolite identification", col)
                                    break
                        
                        # Default to first column if not found
                        if metabolite_column_idx is None:
                            metabolite_column_idx = 0
                            logger.warning("[FEATURE-EXTRACT][METABOLIGHTS] No metabolite column found, using first column")
                        continue
                    
                    if len(parts) > (metabolite_column_idx or 0):
                        metabolite_id = parts[metabolite_column_idx or 0].strip().strip('"').strip("'")
                        
                        if metabolite_id and not metabolite_id.startswith('#') and metabolite_id != "NA":
                            # Skip sample names or other non-metabolite identifiers
                            # Sample names often look like "Sample_01" or are empty
                            if not metabolite_id.lower().startswith("sample") and len(metabolite_id) > 1:
                                # Normalize metabolite identifier
                                from amprenta_rag.ingestion.metabolomics.normalization import normalize_metabolite_name
                                normalized = normalize_metabolite_name(metabolite_id)
                                if normalized and len(normalized) > 1:
                                    metabolite_set.add(normalized)
                        
                        lines_processed += 1
                        
                        if lines_processed % 5000 == 0:
                            logger.info(
                                "[FEATURE-EXTRACT][METABOLIGHTS] Processed %d rows, found %d unique metabolites...",
                                lines_processed,
                                len(metabolite_set),
                            )
                    
                    if lines_processed >= max_rows:
                        break
                
                if lines_processed >= max_rows:
                    break
        
        logger.info(
            "[FEATURE-EXTRACT][METABOLIGHTS] Extracted %d unique metabolites from %d rows",
            len(metabolite_set),
            lines_processed,
        )
        
        return metabolite_set
        
    except Exception as e:
        logger.error("[FEATURE-EXTRACT][METABOLIGHTS] Error extracting metabolites: %r", e)
        return metabolite_set


def extract_mw_metabolites_from_data_endpoint(
    study_id: str,
) -> Set[str]:
    """
    Extract metabolite features from Metabolomics Workbench using the cleaner /data endpoint.
    
    This is much more robust than MetaboLights because:
    - Returns clean JSON directly (no file parsing)
    - No brittle API endpoints
    - More stable (NIH-hosted)
    
    Args:
        study_id: MW study ID (e.g., "ST000001")
        
    Returns:
        Set of normalized metabolite identifiers
    """
    metabolite_set: Set[str] = set()
    
    # Ensure study_id format
    study_id = study_id.upper()
    if not study_id.startswith("ST"):
        logger.error("[FEATURE-EXTRACT][MW] Invalid study ID format: %s", study_id)
        return metabolite_set
    
    logger.info("[FEATURE-EXTRACT][MW] Extracting metabolites from %s using /data endpoint", study_id)
    
    try:
        # Use the cleaner REST API endpoint
        base_url = "https://www.metabolomicsworkbench.org/rest"
        url = f"{base_url}/study/study_id/{study_id}/data"
        
        logger.info("[FEATURE-EXTRACT][MW] Fetching data from: %s", url)
        
        from amprenta_rag.ingestion.repositories import REPOSITORY_USER_AGENT
        headers = {"User-Agent": REPOSITORY_USER_AGENT}
        response = requests.get(url, headers=headers, timeout=30)
        
        # CASE 1: Success
        if response.status_code == 200:
            data = response.json()
            
            # MW API returns a dictionary where keys are row numbers (e.g., '1', '2', '3')
            # Each value is a dict with metabolite information including:
            # - 'metabolite_name': The metabolite name
            # - 'refmet_name': Reference metabolite name
            # - 'metabolite_id': Metabolite ID
            # - 'DATA': Quantification values
            
            from amprenta_rag.ingestion.metabolomics.normalization import normalize_metabolite_name
            
            if isinstance(data, dict):
                # Iterate through all rows in the data dictionary
                for row_key, metabolite_data in data.items():
                    if isinstance(metabolite_data, dict):
                        # Extract metabolite name from various possible fields
                        metabolite_id = (
                            metabolite_data.get('refmet_name') or  # Reference metabolite name (preferred)
                            metabolite_data.get('metabolite_name') or  # Metabolite name
                            metabolite_data.get('metabolite_identification') or
                            metabolite_data.get('database_identifier') or
                            metabolite_data.get('name') or
                            metabolite_data.get('chemical_name')
                        )
                        
                        if metabolite_id:
                            # Normalize metabolite identifier
                            normalized = normalize_metabolite_name(str(metabolite_id))
                            if normalized and len(normalized) > 1:
                                metabolite_set.add(normalized)
            
            elif isinstance(data, list):
                # Alternative structure: list of metabolite dictionaries
                for item in data:
                    if isinstance(item, dict):
                        metabolite_id = (
                            item.get('refmet_name') or
                            item.get('metabolite_name') or
                            item.get('metabolite_identification') or
                            item.get('database_identifier') or
                            item.get('name')
                        )
                        if metabolite_id:
                            from amprenta_rag.ingestion.metabolomics.normalization import normalize_metabolite_name
                            normalized = normalize_metabolite_name(str(metabolite_id))
                            if normalized and len(normalized) > 1:
                                metabolite_set.add(normalized)
        
        # CASE 2: Known error codes
        elif response.status_code in [403, 404]:
            logger.warning("[FEATURE-EXTRACT][MW] Study %s not found or is private (403/404)", study_id)
        
        # CASE 3: Server error
        elif response.status_code >= 500:
            logger.warning("[FEATURE-EXTRACT][MW] Server error (500) for study %s - skipping", study_id)
        
        else:
            logger.warning("[FEATURE-EXTRACT][MW] Unexpected status %d for study %s", response.status_code, study_id)
        
        logger.info(
            "[FEATURE-EXTRACT][MW] Extracted %d unique metabolites from %s",
            len(metabolite_set),
            study_id,
        )
        
        return metabolite_set
        
    except requests.exceptions.RequestException as e:
        logger.error("[FEATURE-EXTRACT][MW] Connection error for study %s: %r", study_id, e)
        return metabolite_set
    except Exception as e:
        logger.error("[FEATURE-EXTRACT][MW] Error extracting metabolites from %s: %r", study_id, e)
        return metabolite_set


def extract_features_from_repository_dataset(
    dataset_id: UUID,
    repository: str,
    study_id: Optional[str] = None,
    download_dir: Path = None,
) -> int:
    """
    Extract features from a repository dataset and link them to Postgres.
    
    Supports:
    - GEO: Genes from expression data (using GEOparse)
    - PRIDE: Proteins from TSV/CSV protein tables
    - MetaboLights: Metabolites from ISA-Tab data files
    - MW: Metabolites from REST API /data endpoint
    
    Args:
        dataset_id: Postgres dataset UUID
        repository: Repository name (GEO, PRIDE, MetaboLights, MW)
        study_id: Optional repository-specific study ID
        download_dir: Directory for temporary files
        
    Returns:
        Number of features linked
    """
    db = next(get_db())
    
    try:
        # Get dataset
        dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
        if not dataset:
            logger.error("[FEATURE-EXTRACT] Dataset %s not found", dataset_id)
            return 0
        
        # Get study_id from dataset if not provided
        if not study_id:
            # Try to extract from dataset metadata or external_ids
            if hasattr(dataset, "external_ids") and dataset.external_ids:
                if isinstance(dataset.external_ids, dict):
                    study_id = dataset.external_ids.get("study_id") or dataset.external_ids.get("accession")
        
        if not study_id:
            logger.warning(
                "[FEATURE-EXTRACT] No study_id provided and couldn't extract from dataset %s",
                dataset_id,
            )
            return 0
        
        logger.info(
            "[FEATURE-EXTRACT] Extracting features from %s study %s for dataset %s",
            repository,
            study_id,
            dataset_id,
        )
        
        # Extract features based on repository type
        features: List[Tuple[str, FeatureType]] = []
        
        if repository.upper() == "GEO":
            # Use GEOparse for cleaner extraction
            gene_set = extract_geo_features_with_geoparse(
                study_id=study_id,
                download_dir=download_dir,
            )
            features = [(gene, FeatureType.GENE) for gene in gene_set]
            
        elif repository.upper() == "PRIDE":
            protein_set = extract_pride_proteins_from_data_files(
                study_id=study_id,
                download_dir=download_dir,
            )
            features = [(protein, FeatureType.PROTEIN) for protein in protein_set]
            
        elif repository.upper() == "METABOLIGHTS":
            metabolite_set = extract_metabolights_metabolites_from_isa_tab(
                study_id=study_id,
                download_dir=download_dir,
            )
            features = [(metabolite, FeatureType.METABOLITE) for metabolite in metabolite_set]
            
        elif repository.upper() in ["MW", "METABOLOMICS WORKBENCH"]:
            metabolite_set = extract_mw_metabolites_from_data_endpoint(
                study_id=study_id,
            )
            features = [(metabolite, FeatureType.METABOLITE) for metabolite in metabolite_set]
        
        else:
            logger.warning("[FEATURE-EXTRACT] Unknown repository: %s", repository)
            return 0
        
        if not features:
            logger.warning("[FEATURE-EXTRACT] No features extracted from %s study %s", repository, study_id)
            return 0
        
        # Link features to dataset
        logger.info(
            "[FEATURE-EXTRACT] Linking %d features to dataset %s",
            len(features),
            dataset_id,
        )
        
        results = batch_link_features_to_dataset_in_postgres(
            features=features,
            dataset_id=dataset_id,
            db=db,
        )
        
        linked_count = sum(1 for v in results.values() if v)
        logger.info(
            "[FEATURE-EXTRACT] Successfully linked %d/%d features to dataset %s",
            linked_count,
            len(features),
            dataset_id,
        )
        
        return linked_count
        
    except Exception as e:
        logger.error(
            "[FEATURE-EXTRACT] Error extracting features for dataset %s: %r",
            dataset_id,
            e,
        )
        return 0
    finally:
        db.close()
