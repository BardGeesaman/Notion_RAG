"""
ID Mapping services for pathway analysis.

Provides functions to map feature identifiers (genes, proteins, metabolites)
to pathway database IDs (KEGG, Reactome, UniProt).
"""

from __future__ import annotations

import re
import time
from typing import Dict, List, Optional, Set

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Create a session with retry strategy
_session = requests.Session()
retry_strategy = Retry(
    total=3,
    backoff_factor=1,
    status_forcelist=[429, 500, 502, 503, 504],
)
adapter = HTTPAdapter(max_retries=retry_strategy)
_session.mount("http://", adapter)
_session.mount("https://", adapter)


# Cache for ID mappings to avoid repeated API calls
_id_mapping_cache: Dict[tuple, Optional[str]] = {}


def _get_cached_or_fetch(cache_key: tuple, fetch_func) -> Optional[str]:
    """Get from cache or fetch using provided function."""
    if cache_key in _id_mapping_cache:
        return _id_mapping_cache[cache_key]
    
    result = fetch_func()
    _id_mapping_cache[cache_key] = result
    return result


# ============================================================================
# UniProt Mapping
# ============================================================================

def map_protein_to_uniprot(protein_id: str) -> Optional[str]:
    """
    Map a protein identifier to UniProt ID.
    
    Handles various input formats:
    - Gene symbols (e.g., "TP53")
    - UniProt IDs (already in correct format)
    - FASTA headers (e.g., "sp|P04637|TP53_HUMAN")
    - Ensembl protein IDs
    
    Args:
        protein_id: Protein identifier
        
    Returns:
        UniProt ID (e.g., "P04637") or None if not found
    """
    if not protein_id:
        return None
    
    # Clean the input
    protein_id = protein_id.strip().upper()
    
    # Check if already a UniProt ID (format: P12345 or P12345-1)
    uniprot_pattern = r'^[OPQ][0-9][A-Z0-9]{3}[0-9](-[0-9]+)?$'
    if re.match(uniprot_pattern, protein_id):
        return protein_id.split('-')[0]  # Return base ID without isoform
    
    # Extract from FASTA header format: sp|P04637|TP53_HUMAN
    fasta_match = re.search(r'\|([OPQ][0-9][A-Z0-9]{3}[0-9])\|', protein_id)
    if fasta_match:
        return fasta_match.group(1)
    
    # Try mapping via UniProt API
    cache_key = ("uniprot", protein_id)
    return _get_cached_or_fetch(cache_key, lambda: _fetch_uniprot_id(protein_id))


def _fetch_uniprot_id(identifier: str) -> Optional[str]:
    """
    Fetch UniProt ID from UniProt REST API.
    
    Uses the UniProt search API to find matching entries.
    """
    try:
        # UniProt REST API: search by gene name or identifier
        url = "https://rest.uniprot.org/uniprotkb/search"
        # Fix query format - use proper field syntax
        query_parts = []
        # Try gene name first (most common)
        query_parts.append(f'gene:"{identifier}"')
        # Try accession if it looks like a UniProt ID
        if re.match(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]', identifier):
            query_parts.append(f'accession:{identifier}')
        # Try as identifier
        query_parts.append(f'id:{identifier}')
        
        params = {
            "query": " OR ".join(query_parts),
            "format": "json",
            "size": 1,
            "fields": "accession",
        }
        
        resp = _session.get(url, params=params, timeout=10)
        resp.raise_for_status()
        
        data = resp.json()
        results = data.get("results", [])
        
        if results:
            uniprot_id = results[0].get("primaryAccession")
            logger.debug(
                "[ANALYSIS][ID-MAP] Mapped %s → UniProt %s",
                identifier,
                uniprot_id,
            )
            return uniprot_id
        
        logger.debug(
            "[ANALYSIS][ID-MAP] No UniProt ID found for %s",
            identifier,
        )
        return None
        
    except Exception as e:
        logger.warning(
            "[ANALYSIS][ID-MAP] Error fetching UniProt ID for %s: %r",
            identifier,
            e,
        )
        return None


# ============================================================================
# KEGG Mapping
# ============================================================================

def map_gene_to_kegg(
    gene_symbol: str,
    organism: str = "hsa",  # hsa = Homo sapiens
) -> Optional[str]:
    """
    Map a gene symbol to KEGG gene ID.
    
    Args:
        gene_symbol: Gene symbol (e.g., "TP53")
        organism: KEGG organism code (default: "hsa" for human)
        
    Returns:
        KEGG gene ID (e.g., "hsa:7157") or None if not found
    """
    if not gene_symbol:
        return None
    
    gene_symbol = gene_symbol.strip().upper()
    
    # Check cache
    cache_key = ("kegg_gene", organism, gene_symbol)
    return _get_cached_or_fetch(
        cache_key,
        lambda: _fetch_kegg_gene_id(gene_symbol, organism),
    )


def _fetch_kegg_gene_id(gene_symbol: str, organism: str) -> Optional[str]:
    """
    Fetch KEGG gene ID from KEGG REST API.
    
    Uses the KEGG conv API to convert gene symbols to KEGG IDs.
    """
    try:
        # KEGG conv API: ncbi-geneid:7157 → kegg
        # First, try to get NCBI gene ID via MyGene.info
        ncbi_id = _get_ncbi_gene_id(gene_symbol)
        
        if ncbi_id:
            # Convert NCBI ID to KEGG
            url = f"https://rest.kegg.jp/conv/{organism}/ncbi-geneid:{ncbi_id}"
            resp = _session.get(url, timeout=10)
            
            if resp.status_code == 200 and resp.text.strip():
                # Response format: ncbi-geneid:7157	hsa:7157 (note: columns can be in either order)
                lines = resp.text.strip().split("\n")
                if lines:
                    parts = lines[0].split("\t")
                    # Find the KEGG ID (format: hsa:####)
                    for part in parts:
                        if ":" in part and part.startswith(organism):
                            kegg_id = part
                            logger.debug(
                                "[ANALYSIS][ID-MAP] Mapped %s → KEGG %s",
                                gene_symbol,
                                kegg_id,
                            )
                            return kegg_id
        
        # Fallback: try direct KEGG search
        url = f"https://rest.kegg.jp/find/{organism}/gene"
        params = {"query": gene_symbol}
        resp = _session.get(url, params=params, timeout=10)
        
        if resp.status_code == 200 and resp.text.strip():
            # Response format: hsa:7157	TP53; tumor protein p53
            lines = resp.text.strip().split("\n")
            for line in lines:
                if gene_symbol.upper() in line.upper():
                    kegg_id = line.split("\t")[0]
                    logger.debug(
                        "[ANALYSIS][ID-MAP] Mapped %s → KEGG %s (direct)",
                        gene_symbol,
                        kegg_id,
                    )
                    return kegg_id
        
        logger.debug(
            "[ANALYSIS][ID-MAP] No KEGG gene ID found for %s",
            gene_symbol,
        )
        return None
        
    except Exception as e:
        logger.warning(
            "[ANALYSIS][ID-MAP] Error fetching KEGG gene ID for %s: %r",
            gene_symbol,
            e,
        )
        return None


def map_metabolite_to_kegg(metabolite_name: str) -> Optional[str]:
    """
    Map a metabolite name to KEGG compound ID.
    
    Args:
        metabolite_name: Metabolite name (e.g., "Glucose", "ATP")
        
    Returns:
        KEGG compound ID (e.g., "C00031") or None if not found
    """
    if not metabolite_name:
        return None
    
    # Clean metabolite name
    metabolite_name = metabolite_name.strip()
    
    # Check cache
    cache_key = ("kegg_compound", metabolite_name)
    return _get_cached_or_fetch(
        cache_key,
        lambda: _fetch_kegg_compound_id(metabolite_name),
    )


def _fetch_kegg_compound_id(metabolite_name: str) -> Optional[str]:
    """
    Fetch KEGG compound ID from KEGG REST API.
    """
    try:
        # KEGG find API: search compounds by name
        # Note: KEGG REST API doesn't support query params in the same way
        # We need to use the URL directly
        url = f"https://rest.kegg.jp/find/compound/{metabolite_name}"
        resp = _session.get(url, timeout=10)
        
        if resp.status_code == 200 and resp.text.strip():
            # Response format: C00031	D-Glucose
            lines = resp.text.strip().split("\n")
            for line in lines:
                parts = line.split("\t")
                if len(parts) >= 2:
                    kegg_id = parts[0]
                    name = parts[1]
                    # Check if metabolite name matches (case-insensitive)
                    if metabolite_name.upper() in name.upper() or name.upper() in metabolite_name.upper():
                        logger.debug(
                            "[ANALYSIS][ID-MAP] Mapped %s → KEGG %s",
                            metabolite_name,
                            kegg_id,
                        )
                        return kegg_id
        
        logger.debug(
            "[ANALYSIS][ID-MAP] No KEGG compound ID found for %s",
            metabolite_name,
        )
        return None
        
    except Exception as e:
        logger.warning(
            "[ANALYSIS][ID-MAP] Error fetching KEGG compound ID for %s: %r",
            metabolite_name,
            e,
        )
        return None


def map_protein_to_kegg(protein_id: str, organism: str = "hsa") -> Optional[str]:
    """
    Map a protein identifier to KEGG gene ID (via UniProt).
    
    Args:
        protein_id: Protein identifier
        organism: KEGG organism code
        
    Returns:
        KEGG gene ID or None if not found
    """
    # First, get UniProt ID
    uniprot_id = map_protein_to_uniprot(protein_id)
    if not uniprot_id:
        return None
    
    # Then, convert UniProt to KEGG
    try:
        url = f"https://rest.kegg.jp/conv/{organism}/uniprot:{uniprot_id}"
        resp = _session.get(url, timeout=10)
        
        if resp.status_code == 200 and resp.text.strip():
            # Response format: hsa:7157	uniprot:P04637
            lines = resp.text.strip().split("\n")
            if lines:
                kegg_id = lines[0].split("\t")[0]
                logger.debug(
                    "[ANALYSIS][ID-MAP] Mapped %s (UniProt %s) → KEGG %s",
                    protein_id,
                    uniprot_id,
                    kegg_id,
                )
                return kegg_id
        
        return None
        
    except Exception as e:
        logger.warning(
            "[ANALYSIS][ID-MAP] Error converting UniProt %s to KEGG: %r",
            uniprot_id,
            e,
        )
        return None


# ============================================================================
# Reactome Mapping
# ============================================================================

def map_gene_to_reactome(gene_symbol: str) -> Optional[str]:
    """
    Map a gene symbol to Reactome stable identifier.
    
    Args:
        gene_symbol: Gene symbol
        
    Returns:
        Reactome stable identifier or None if not found
    """
    if not gene_symbol:
        return None
    
    gene_symbol = gene_symbol.strip().upper()
    
    # Check cache
    cache_key = ("reactome_gene", gene_symbol)
    return _get_cached_or_fetch(
        cache_key,
        lambda: _fetch_reactome_gene_id(gene_symbol),
    )


def _fetch_reactome_gene_id(gene_symbol: str) -> Optional[str]:
    """
    Fetch Reactome stable identifier from Reactome API.
    """
    try:
        # Reactome identifier mapping API
        url = "https://reactome.org/ContentService/data/query/ids/map"
        params = {
            "species": "9606",  # Homo sapiens
            "ids": gene_symbol,
        }
        
        resp = _session.post(url, json=params, timeout=10)
        resp.raise_for_status()
        
        data = resp.json()
        if data and len(data) > 0:
            stable_id = data[0].get("stableIdentifier", {}).get("identifier")
            if stable_id:
                logger.debug(
                    "[ANALYSIS][ID-MAP] Mapped %s → Reactome %s",
                    gene_symbol,
                    stable_id,
                )
                return stable_id
        
        logger.debug(
            "[ANALYSIS][ID-MAP] No Reactome ID found for %s",
            gene_symbol,
        )
        return None
        
    except Exception as e:
        logger.warning(
            "[ANALYSIS][ID-MAP] Error fetching Reactome ID for %s: %r",
            gene_symbol,
            e,
        )
        return None


def map_protein_to_reactome(protein_id: str) -> Optional[str]:
    """
    Map a protein identifier to Reactome stable identifier (via UniProt).
    """
    # First, get UniProt ID
    uniprot_id = map_protein_to_uniprot(protein_id)
    if not uniprot_id:
        return None
    
    # Reactome uses UniProt IDs directly in many cases
    # Try to find Reactome entity for this UniProt ID
    try:
        url = f"https://reactome.org/ContentService/data/query/map/uniprot/{uniprot_id}"
        resp = _session.get(url, timeout=10)
        
        if resp.status_code == 200:
            data = resp.json()
            if data and len(data) > 0:
                stable_id = data[0].get("stableIdentifier", {}).get("identifier")
                if stable_id:
                    logger.debug(
                        "[ANALYSIS][ID-MAP] Mapped %s (UniProt %s) → Reactome %s",
                        protein_id,
                        uniprot_id,
                        stable_id,
                    )
                    return stable_id
        
        return None
        
    except Exception as e:
        logger.warning(
            "[ANALYSIS][ID-MAP] Error mapping UniProt %s to Reactome: %r",
            uniprot_id,
            e,
        )
        return None


# ============================================================================
# Helper Functions
# ============================================================================

def _get_ncbi_gene_id(gene_symbol: str) -> Optional[str]:
    """
    Get NCBI gene ID from MyGene.info API.
    
    This is used as an intermediate step for KEGG mapping.
    """
    try:
        url = "https://mygene.info/v3/query"
        params = {
            "q": gene_symbol,
            "species": "human",
            "fields": "entrezgene",
            "size": 1,
        }
        
        resp = _session.get(url, params=params, timeout=10)
        resp.raise_for_status()
        
        data = resp.json()
        hits = data.get("hits", [])
        
        if hits:
            ncbi_id = hits[0].get("entrezgene")
            if ncbi_id:
                return str(ncbi_id)
        
        return None
        
    except Exception as e:
        logger.debug(
            "[ANALYSIS][ID-MAP] Error fetching NCBI gene ID for %s: %r",
            gene_symbol,
            e,
        )
        return None


def batch_map_features_to_pathway_ids(
    features_by_type: Dict[str, Set[str]],
    pathway_source: str = "KEGG",
    organism: str = "hsa",
) -> Dict[str, Dict[str, Optional[str]]]:
    """
    Batch map features to pathway database IDs.
    
    Args:
        features_by_type: Dictionary mapping feature types to sets of feature names
        pathway_source: "KEGG" or "Reactome"
        organism: KEGG organism code (for KEGG only)
        
    Returns:
        Dictionary mapping feature_type -> feature_name -> pathway_id
    """
    results: Dict[str, Dict[str, Optional[str]]] = {}
    
    for feature_type, features in features_by_type.items():
        results[feature_type] = {}
        
        for feature in features:
            if pathway_source == "KEGG":
                if feature_type == "gene":
                    pathway_id = map_gene_to_kegg(feature, organism=organism)
                elif feature_type == "protein":
                    pathway_id = map_protein_to_kegg(feature, organism=organism)
                elif feature_type == "metabolite":
                    pathway_id = map_metabolite_to_kegg(feature)
                else:
                    pathway_id = None
            elif pathway_source == "Reactome":
                if feature_type == "gene":
                    pathway_id = map_gene_to_reactome(feature)
                elif feature_type == "protein":
                    pathway_id = map_protein_to_reactome(feature)
                else:
                    pathway_id = None
            else:
                pathway_id = None
            
            results[feature_type][feature] = pathway_id
            
            # Rate limiting: be nice to APIs
            time.sleep(0.1)
    
    return results

