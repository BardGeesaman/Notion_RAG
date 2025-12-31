"""
ID Mapping Services for Pathway Analysis.

Provides mapping functions to convert feature identifiers to pathway database IDs:
- Genes → KEGG IDs, Reactome symbols
- Proteins → UniProt IDs, KEGG IDs, Reactome IDs
- Metabolites → KEGG Compound IDs

Uses external services: UniProt REST API, KEGG API, Reactome API, MyGene.info

KEGG API Rate Limits (Important):
- No official rate limit published by KEGG
- Recommended: Maximum 3 requests/second
- Current implementation: 0.1s sleep between requests (line 413)
- Bulk downloads require paid KEGG license ($2,000+/year)
- Academic use only without subscription
- This implementation uses on-demand caching to comply with terms
"""

from __future__ import annotations

import re
import time
from typing import Dict, Optional, Set, Tuple

import requests

from amprenta_rag.logging_utils import get_logger
from amprenta_rag.services.id_mapping_service import (
    get_mapping,
    save_mapping,
)

logger = get_logger(__name__)

# Cache for ID mappings to avoid repeated API calls
_id_mapping_cache: Dict[Tuple[str, str, str], Optional[str]] = {}


def map_protein_to_uniprot(protein_id: str) -> Optional[str]:
    """
    Map protein identifier to UniProt ID.

    Handles various protein ID formats:
    - UniProt IDs (P12345) - returned as-is
    - Gene symbols (TP53) - mapped via UniProt API
    - Ensembl protein IDs (ENSP00000...) - mapped via UniProt API

    Args:
        protein_id: Protein identifier

    Returns:
        UniProt ID or None if mapping fails
    """
    # 1. Check in-memory cache first (fast path)
    cache_key = ("uniprot", protein_id, "")
    if cache_key in _id_mapping_cache:
        return _id_mapping_cache[cache_key]

    # 2. Check database
    db_result = get_mapping("protein", protein_id, "uniprot")
    if db_result:
        _id_mapping_cache[cache_key] = db_result
        return db_result

    # 3. Check if already a UniProt ID (format: [OPQ][0-9][A-Z0-9]{3}[0-9])
    uniprot_pattern = r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"
    if re.match(uniprot_pattern, protein_id):
        _id_mapping_cache[cache_key] = protein_id
        return protein_id

    # 4. Fall back to API (existing logic)
    try:
        # Try UniProt mapping service
        # https://www.uniprot.org/help/id_mapping
        url = "https://rest.uniprot.org/idmapping/run"

        # Try as gene name first, then as other ID types
        for from_db in ["Gene_Name", "Ensembl_Protein", "RefSeq_Protein"]:
            try:
                data = {
                    "from": from_db,
                    "to": "UniProtKB",
                    "ids": protein_id,
                }

                resp = requests.post(url, data=data, timeout=10)

                if resp.status_code == 200:
                    job_id = resp.json().get("jobId")
                    if job_id:
                        # Poll for results
                        results_url = f"https://rest.uniprot.org/idmapping/results/{job_id}"

                        for _ in range(10):  # Max 10 attempts
                            time.sleep(0.5)  # Wait for job completion
                            results_resp = requests.get(results_url, timeout=10)

                            if results_resp.status_code == 200:
                                results = results_resp.json()
                                if results.get("results"):
                                    uniprot_id = results["results"][0].get("to", {}).get("primaryAccession")
                                    if uniprot_id:
                                        logger.debug(
                                            "[ANALYSIS][ID-MAP] Mapped protein %s -> UniProt %s",
                                            protein_id,
                                            uniprot_id,
                                        )
                                        _id_mapping_cache[cache_key] = uniprot_id
                                        # 5. Cache successful API result in database
                                        save_mapping("protein", protein_id, "uniprot", uniprot_id, ttl_days=90)
                                        return uniprot_id
                                break

            except Exception as e:
                logger.debug(
                    "[ANALYSIS][ID-MAP] UniProt mapping attempt failed for %s from %s: %r",
                    protein_id,
                    from_db,
                    e,
                )
                continue

        # If all attempts fail, return None
        logger.debug(
            "[ANALYSIS][ID-MAP] Could not map protein %s to UniProt",
            protein_id,
        )
        _id_mapping_cache[cache_key] = None
        return None

    except Exception as e:
        logger.warning(
            "[ANALYSIS][ID-MAP] Error mapping protein %s to UniProt: %r",
            protein_id,
            e,
        )
        _id_mapping_cache[cache_key] = None
        return None


def map_gene_to_kegg(gene_symbol: str, organism: str = "hsa") -> Optional[str]:
    """
    Map gene symbol to KEGG gene ID.

    Args:
        gene_symbol: Gene symbol (e.g., "TP53")
        organism: KEGG organism code (default: "hsa" for human)

    Returns:
        KEGG gene ID (e.g., "hsa:7157") or None if mapping fails
    """
    # 1. Check in-memory cache first (fast path)
    cache_key = ("kegg_gene", gene_symbol, organism)
    if cache_key in _id_mapping_cache:
        return _id_mapping_cache[cache_key]

    # 2. Check database
    organism_key = "human" if organism == "hsa" else organism
    db_result = get_mapping("gene", gene_symbol, "kegg_gene", organism_key)
    if db_result:
        _id_mapping_cache[cache_key] = db_result
        return db_result

    # 3. Fall back to API (existing logic)
    try:
        # Use KEGG find API to search by gene symbol
        # https://www.kegg.jp/kegg/rest/keggapi.html
        url = f"https://rest.kegg.jp/find/{organism}/{gene_symbol}"

        resp = requests.get(url, timeout=10)

        if resp.status_code == 200 and resp.text.strip():
            # Parse response (format: kegg_id\tdescription)
            lines = resp.text.strip().split("\n")
            if lines:
                first_line = lines[0]
                kegg_id = first_line.split("\t")[0]
                logger.debug(
                    "[ANALYSIS][ID-MAP] Mapped gene %s -> KEGG %s",
                    gene_symbol,
                    kegg_id,
                )
                _id_mapping_cache[cache_key] = kegg_id
                # Cache successful API result in database
                save_mapping("gene", gene_symbol, "kegg_gene", kegg_id, organism_key, ttl_days=90)
                return kegg_id

        logger.debug(
            "[ANALYSIS][ID-MAP] Could not map gene %s to KEGG",
            gene_symbol,
        )
        _id_mapping_cache[cache_key] = None
        return None

    except Exception as e:
        logger.warning(
            "[ANALYSIS][ID-MAP] Error mapping gene %s to KEGG: %r",
            gene_symbol,
            e,
        )
        _id_mapping_cache[cache_key] = None
        return None


def map_protein_to_kegg(protein_id: str, organism: str = "hsa") -> Optional[str]:
    """
    Map protein identifier to KEGG gene ID.

    Strategy:
    1. Map protein to UniProt ID
    2. Map UniProt ID to KEGG gene ID using KEGG conv API

    Args:
        protein_id: Protein identifier
        organism: KEGG organism code (default: "hsa" for human)

    Returns:
        KEGG gene ID or None if mapping fails
    """
    # 1. Check in-memory cache first (fast path)
    cache_key = ("kegg_protein", protein_id, organism)
    if cache_key in _id_mapping_cache:
        return _id_mapping_cache[cache_key]

    # 2. Check database
    organism_key = "human" if organism == "hsa" else organism
    db_result = get_mapping("protein", protein_id, "kegg_gene", organism_key)
    if db_result:
        _id_mapping_cache[cache_key] = db_result
        return db_result

    # 3. Fall back to API (existing logic)
    try:
        # First get UniProt ID
        uniprot_id = map_protein_to_uniprot(protein_id)
        if not uniprot_id:
            _id_mapping_cache[cache_key] = None
            return None

        # Use KEGG conv API to map UniProt -> KEGG
        # https://rest.kegg.jp/conv/genes/uniprot:{uniprot_id}
        url = f"https://rest.kegg.jp/conv/{organism}/uniprot:{uniprot_id}"

        resp = requests.get(url, timeout=10)

        if resp.status_code == 200 and resp.text.strip():
            # Parse response (format: uniprot_id\tkegg_id)
            lines = resp.text.strip().split("\n")
            if lines:
                first_line = lines[0]
                parts = first_line.split("\t")
                if len(parts) >= 2:
                    kegg_id = parts[1]
                    logger.debug(
                        "[ANALYSIS][ID-MAP] Mapped protein %s -> KEGG %s (via UniProt %s)",
                        protein_id,
                        kegg_id,
                        uniprot_id,
                    )
                    _id_mapping_cache[cache_key] = kegg_id
                    # Cache successful API result in database
                    save_mapping("protein", protein_id, "kegg_gene", kegg_id, organism_key, ttl_days=90)
                    return kegg_id

        logger.debug(
            "[ANALYSIS][ID-MAP] Could not map protein %s to KEGG",
            protein_id,
        )
        _id_mapping_cache[cache_key] = None
        return None

    except Exception as e:
        logger.warning(
            "[ANALYSIS][ID-MAP] Error mapping protein %s to KEGG: %r",
            protein_id,
            e,
        )
        _id_mapping_cache[cache_key] = None
        return None


def map_metabolite_to_kegg(metabolite_name: str) -> Optional[str]:
    """
    Map metabolite name to KEGG Compound ID.

    Args:
        metabolite_name: Metabolite name (e.g., "glucose", "ATP")

    Returns:
        KEGG Compound ID (e.g., "cpd:C00031") or None if mapping fails
    """
    # 1. Check in-memory cache first (fast path)
    cache_key = ("kegg_metabolite", metabolite_name, "")
    if cache_key in _id_mapping_cache:
        return _id_mapping_cache[cache_key]

    # 2. Check database
    db_result = get_mapping("metabolite", metabolite_name, "kegg_compound")
    if db_result:
        _id_mapping_cache[cache_key] = db_result
        return db_result

    # 3. Fall back to API (existing logic)
    try:
        # Use KEGG find API to search compound database
        url = f"https://rest.kegg.jp/find/compound/{metabolite_name}"

        resp = requests.get(url, timeout=10)

        if resp.status_code == 200 and resp.text.strip():
            # Parse response (format: cpd:C00001\tname; synonyms)
            lines = resp.text.strip().split("\n")

            # Look for exact or best match
            for line in lines:
                parts = line.split("\t")
                if len(parts) >= 2:
                    compound_id = parts[0]
                    names = parts[1].lower()

                    # Check if metabolite name is in compound names
                    if metabolite_name.lower() in names:
                        logger.debug(
                            "[ANALYSIS][ID-MAP] Mapped metabolite %s -> KEGG %s",
                            metabolite_name,
                            compound_id,
                        )
                        _id_mapping_cache[cache_key] = compound_id
                        # Cache successful API result in database
                        save_mapping("metabolite", metabolite_name, "kegg_compound", compound_id, ttl_days=90)
                        return compound_id

            # If no exact match, return first result
            if lines:
                first_line = lines[0]
                compound_id = first_line.split("\t")[0]
                logger.debug(
                    "[ANALYSIS][ID-MAP] Mapped metabolite %s -> KEGG %s (fuzzy match)",
                    metabolite_name,
                    compound_id,
                )
                _id_mapping_cache[cache_key] = compound_id
                # Cache successful API result in database
                save_mapping("metabolite", metabolite_name, "kegg_compound", compound_id, ttl_days=90)
                return compound_id

        logger.debug(
            "[ANALYSIS][ID-MAP] Could not map metabolite %s to KEGG",
            metabolite_name,
        )
        _id_mapping_cache[cache_key] = None
        return None

    except Exception as e:
        logger.warning(
            "[ANALYSIS][ID-MAP] Error mapping metabolite %s to KEGG: %r",
            metabolite_name,
            e,
        )
        _id_mapping_cache[cache_key] = None
        return None


def map_gene_to_reactome(gene_symbol: str) -> Optional[str]:
    """
    Map gene symbol to Reactome-compatible identifier.

    Reactome uses gene symbols directly for human genes.

    Args:
        gene_symbol: Gene symbol (e.g., "TP53")

    Returns:
        Gene symbol (Reactome uses symbols directly) or None
    """
    # Reactome uses gene symbols directly, so return as-is if valid
    # Valid gene symbols are typically alphanumeric with underscores/hyphens
    if re.match(r"^[A-Z0-9_-]+$", gene_symbol, re.IGNORECASE):
        return gene_symbol

    logger.debug(
        "[ANALYSIS][ID-MAP] Invalid gene symbol format for Reactome: %s",
        gene_symbol,
    )
    return None


def map_protein_to_reactome(protein_id: str) -> Optional[str]:
    """
    Map protein identifier to Reactome-compatible identifier.

    Reactome uses UniProt IDs for proteins.

    Args:
        protein_id: Protein identifier

    Returns:
        UniProt ID or None if mapping fails
    """
    # Reactome uses UniProt IDs, so map to UniProt
    return map_protein_to_uniprot(protein_id)


def batch_map_features_to_pathway_ids(
    features: Set[str],
    feature_type: str,
    pathway_source: str = "KEGG",
    organism: str = "hsa",
) -> Dict[str, Optional[str]]:
    """
    Batch map features to pathway database IDs.

    Args:
        features: Set of feature identifiers
        feature_type: One of "gene", "protein", "metabolite", "lipid"
        pathway_source: "KEGG" or "Reactome"
        organism: KEGG organism code (default: "hsa" for human)

    Returns:
        Dictionary mapping feature_id -> pathway_db_id
    """
    logger.info(
        "[ANALYSIS][ID-MAP] Batch mapping %d %s features to %s IDs",
        len(features),
        feature_type,
        pathway_source,
    )

    mapping: Dict[str, Optional[str]] = {}

    for feature in features:
        if pathway_source == "KEGG":
            if feature_type == "gene":
                mapped_id = map_gene_to_kegg(feature, organism)
            elif feature_type == "protein":
                mapped_id = map_protein_to_kegg(feature, organism)
            elif feature_type == "metabolite":
                mapped_id = map_metabolite_to_kegg(feature)
            else:
                logger.warning(
                    "[ANALYSIS][ID-MAP] KEGG mapping not supported for feature_type '%s'",
                    feature_type,
                )
                mapped_id = None

        elif pathway_source == "Reactome":
            if feature_type == "gene":
                mapped_id = map_gene_to_reactome(feature)
            elif feature_type == "protein":
                mapped_id = map_protein_to_reactome(feature)
            else:
                logger.warning(
                    "[ANALYSIS][ID-MAP] Reactome mapping not supported for feature_type '%s'",
                    feature_type,
                )
                mapped_id = None

        else:
            logger.warning(
                "[ANALYSIS][ID-MAP] Unknown pathway source '%s'",
                pathway_source,
            )
            mapped_id = None

        mapping[feature] = mapped_id

        # Rate limiting to respect API limits
        time.sleep(0.1)

    successful_mappings = sum(1 for v in mapping.values() if v is not None)
    logger.info(
        "[ANALYSIS][ID-MAP] Successfully mapped %d/%d features to %s IDs",
        successful_mappings,
        len(features),
        pathway_source,
    )

    return mapping


def clear_id_mapping_cache():
    """Clear the ID mapping cache."""
    global _id_mapping_cache
    _id_mapping_cache.clear()
    logger.info("[ANALYSIS][ID-MAP] Cleared ID mapping cache")


def get_cache_stats() -> Dict[str, int]:
    """
    Get statistics about the ID mapping cache.

    Returns:
        Dictionary with cache statistics
    """
    total = len(_id_mapping_cache)
    successful = sum(1 for v in _id_mapping_cache.values() if v is not None)
    failed = total - successful

    return {
        "total_cached": total,
        "successful_mappings": successful,
        "failed_mappings": failed,
    }
