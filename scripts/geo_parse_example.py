#!/usr/bin/env python3
"""
Example script for searching and extracting data from NCBI GEO using GEOparse.

This demonstrates the recommended approach using GEOparse library.
"""

import os
import sys
import time
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import GEOparse
import pandas as pd
from Bio import Entrez

from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def configure_entrez():
    """Configure Entrez with email and API key from config."""
    cfg = get_config()
    
    # Get email and API key from config
    email = getattr(cfg, "ncbi_email", None) or os.getenv("NCBI_EMAIL", "")
    api_key = getattr(cfg, "geo_api_key", None) or os.getenv("GEO_API_KEY", "")
    
    if email:
        Entrez.email = email
        logger.info("[GEO] Configured Entrez.email = %s", email)
    else:
        logger.warning("[GEO] No email configured. Set NCBI_EMAIL in config or .env")
    
    if api_key:
        Entrez.api_key = api_key
        logger.info("[GEO] Configured Entrez.api_key")
    
    # Rate limiting
    return 0.1 if api_key else 0.34


def search_geo_studies(keyword: str, max_results: int = 5) -> list[str]:
    """
    Step 1: Find GSE IDs using Biopython Entrez.
    
    Args:
        keyword: Search term (e.g., "breast cancer")
        max_results: Maximum number of results
        
    Returns:
        List of GSE IDs (e.g., ["GSE1001", "GSE1234"])
    """
    logger.info("[GEO] Searching GEO for: %s (max_results=%d)", keyword, max_results)
    
    try:
        # Search for Series (GSE) specifically
        term = f"{keyword} AND GSE[Entry Type]"
        
        handle = Entrez.esearch(db="gds", term=term, retmax=max_results)
        results = Entrez.read(handle)
        handle.close()
        
        # ID list from Entrez is internal IDs, we need to fetch the GSE accession
        id_list = results['IdList']
        if not id_list:
            logger.warning("[GEO] No studies found for keyword: %s", keyword)
            return []
        
        # Fetch Summary to get the actual GSE accession (e.g., GSE1001)
        # Note: This step is needed because 'IdList' returns internal numbers like '20001001'
        time.sleep(0.1)  # Rate limit
        
        handle = Entrez.esummary(db="gds", id=",".join(id_list))
        summaries = Entrez.read(handle)
        handle.close()
        
        gse_ids = []
        for summary in summaries:
            # The accession is usually in the 'Accession' field e.g., 'GSE12345'
            accession = summary.get('Accession', '')
            if accession.startswith('GSE'):
                gse_ids.append(accession)
        
        logger.info("[GEO] Found %d studies: %s", len(gse_ids), gse_ids[:5])
        return gse_ids
        
    except Exception as e:
        logger.error("[GEO] Search error: %r", e)
        return []


def extract_geo_data(gse_id: str, cache_dir: str = "./geo_cache") -> tuple:
    """
    Step 2: Download and Extract Data using GEOparse.
    
    Args:
        gse_id: GEO Series ID (e.g., "GSE1001")
        cache_dir: Directory to cache downloaded files
        
    Returns:
        Tuple of (clinical_data DataFrame, expression_matrix DataFrame) or (None, None) on error
    """
    logger.info("[GEO] Downloading and parsing %s...", gse_id)
    
    try:
        # Create cache directory
        os.makedirs(cache_dir, exist_ok=True)
        
        # Download and parse using GEOparse
        # destdir saves the file locally so you don't download it twice
        gse = GEOparse.get_GEO(geo=gse_id, destdir=cache_dir, silent=True)
        
        # 1. Extract Clinical Data (Sample Metadata)
        # Rows = Samples, Columns = Attributes (Age, Treatment, etc.)
        clinical_data = gse.phenotype_data
        
        # 2. Extract Expression Matrix
        # Rows = Probes/Genes, Columns = Samples
        # 'VALUE' is the standard column for normalized data
        expression_matrix = gse.pivot_samples('VALUE')
        
        logger.info("[GEO] Successfully loaded %s", gse_id)
        logger.info("[GEO]   Samples: %d", len(clinical_data))
        logger.info("[GEO]   Probes/Genes: %d", len(expression_matrix))
        
        return clinical_data, expression_matrix
        
    except Exception as e:
        logger.error("[GEO] Extraction error for %s: %r", gse_id, e)
        return None, None


def extract_genes_from_geo_parse(gse_id: str, cache_dir: str = "./geo_cache") -> set[str]:
    """
    Extract gene identifiers from GEO study using GEOparse.
    
    This is a cleaner approach than custom parsing.
    
    Args:
        gse_id: GEO Series ID
        cache_dir: Cache directory for downloaded files
        
    Returns:
        Set of normalized gene identifiers
    """
    from amprenta_rag.ingestion.transcriptomics.normalization import normalize_gene_identifier
    
    gene_set: set[str] = set()
    
    try:
        os.makedirs(cache_dir, exist_ok=True)
        
        # Download and parse using GEOparse
        gse = GEOparse.get_GEO(geo=gse_id, destdir=cache_dir, silent=True)
        
        # Get expression matrix
        expression_matrix = gse.pivot_samples('VALUE')
        
        if expression_matrix is not None and not expression_matrix.empty:
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
                "[GEO] Extracted %d unique genes from %s using GEOparse",
                len(gene_set),
                gse_id,
            )
        
    except Exception as e:
        logger.error("[GEO] Error extracting genes from %s with GEOparse: %r", gse_id, e)
    
    return gene_set


def main():
    """Example usage."""
    # Configure Entrez
    rate_delay = configure_entrez()
    
    # Create cache dir
    cache_dir = "./geo_cache"
    os.makedirs(cache_dir, exist_ok=True)
    
    # Example 1: Search for studies
    keyword = "breast cancer"
    studies = search_geo_studies(keyword, max_results=3)
    
    print(f"\n✅ Found {len(studies)} studies: {studies}")
    
    # Example 2: Extract data from first study
    if studies:
        target_study = studies[0]
        print(f"\n📊 Extracting data from {target_study}...")
        
        meta, data = extract_geo_data(target_study, cache_dir=cache_dir)
        
        if data is not None:
            print("\n--- Preview: Clinical Data ---")
            print(meta.iloc[:, :5].head())  # Print first 5 cols only
            
            print("\n--- Preview: Expression Matrix ---")
            print(data.head())
            
            # Example 3: Extract genes using GEOparse
            print(f"\n🧬 Extracting gene identifiers from {target_study}...")
            genes = extract_genes_from_geo_parse(target_study, cache_dir=cache_dir)
            print(f"✅ Extracted {len(genes)} unique genes")
            if genes:
                print(f"Sample genes: {list(genes)[:10]}")


if __name__ == "__main__":
    main()

