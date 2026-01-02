"""STRING database API client for protein-protein interaction data."""

from __future__ import annotations

import logging
import time
from typing import Any, Optional

import requests
from pydantic import BaseModel
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type

logger = logging.getLogger(__name__)

STRING_API_BASE = "https://string-db.org/api"
RATE_LIMIT_DELAY = 1.0  # 1 request per second
DEFAULT_SPECIES = 9606  # Human


class StringInteraction(BaseModel):
    """Protein-protein interaction from STRING."""
    protein_a: str
    protein_b: str
    score: int  # Combined score 0-1000
    gene_a: Optional[str] = None
    gene_b: Optional[str] = None
    # Evidence channels
    nscore: int = 0  # Neighborhood
    fscore: int = 0  # Fusion
    pscore: int = 0  # Phylogenetic
    ascore: int = 0  # Coexpression
    escore: int = 0  # Experimental
    dscore: int = 0  # Database
    tscore: int = 0  # Textmining


class StringProtein(BaseModel):
    """Protein metadata from STRING."""
    string_id: str
    preferred_name: str  # Gene symbol
    annotation: Optional[str] = None


@retry(
    wait=wait_exponential(multiplier=2, min=4, max=60),
    stop=stop_after_attempt(3),
    retry=retry_if_exception_type((requests.HTTPError, requests.ConnectionError)),
)
def _string_request(endpoint: str, params: dict, timeout: int = 30) -> Any:
    """Make STRING API request with retry logic."""
    time.sleep(RATE_LIMIT_DELAY)
    url = f"{STRING_API_BASE}/{endpoint}"
    response = requests.get(url, params=params, timeout=timeout)
    response.raise_for_status()
    return response.json()


def resolve_gene_symbols(
    symbols: list[str],
    species: int = DEFAULT_SPECIES,
) -> dict[str, str]:
    """
    Resolve gene symbols to STRING protein IDs.
    
    Args:
        symbols: List of gene symbols (e.g., ["TP53", "BRCA1"])
        species: NCBI taxonomy ID (9606=human, 10090=mouse)
    
    Returns:
        Dict mapping gene symbol to STRING ID
    """
    if not symbols:
        return {}
    
    identifiers = "\r".join(symbols)
    params = {
        "identifiers": identifiers,
        "species": species,
        "limit": 1,  # Best match per query
        "echo_query": 1,
    }
    
    try:
        data = _string_request("json/get_string_ids", params)
        result = {}
        for item in data:
            query = item.get("queryItem", "").upper()
            string_id = item.get("stringId")
            if query and string_id:
                result[query] = string_id
        return result
    except Exception as e:
        logger.warning(f"Failed to resolve gene symbols: {e}")
        return {}


def get_interactions(
    proteins: list[str],
    species: int = DEFAULT_SPECIES,
    min_score: int = 400,
    network_type: str = "functional",
) -> list[StringInteraction]:
    """
    Fetch protein-protein interactions from STRING.
    
    Args:
        proteins: List of gene symbols or STRING IDs
        species: NCBI taxonomy ID
        min_score: Minimum combined score (0-1000, default 400 = medium confidence)
        network_type: "functional" or "physical"
    
    Returns:
        List of StringInteraction objects
    """
    if not proteins:
        return []
    
    identifiers = "\r".join(proteins)
    params = {
        "identifiers": identifiers,
        "species": species,
        "required_score": min_score,
        "network_type": network_type,
    }
    
    try:
        data = _string_request("json/network", params, timeout=60)
        interactions = []
        for item in data:
            interactions.append(StringInteraction(
                protein_a=item.get("stringId_A", ""),
                protein_b=item.get("stringId_B", ""),
                gene_a=item.get("preferredName_A"),
                gene_b=item.get("preferredName_B"),
                score=int(item.get("score", 0) * 1000),
                nscore=int(item.get("nscore", 0) * 1000),
                fscore=int(item.get("fscore", 0) * 1000),
                pscore=int(item.get("pscore", 0) * 1000),
                ascore=int(item.get("ascore", 0) * 1000),
                escore=int(item.get("escore", 0) * 1000),
                dscore=int(item.get("dscore", 0) * 1000),
                tscore=int(item.get("tscore", 0) * 1000),
            ))
        logger.info(f"STRING returned {len(interactions)} interactions for {len(proteins)} proteins")
        return interactions
    except Exception as e:
        logger.error(f"STRING network request failed: {e}")
        return []


def get_interaction_partners(
    protein: str,
    species: int = DEFAULT_SPECIES,
    limit: int = 50,
    min_score: int = 400,
) -> list[StringInteraction]:
    """
    Get interaction partners for a single protein.
    
    Args:
        protein: Gene symbol or STRING ID
        species: NCBI taxonomy ID
        limit: Maximum number of partners
        min_score: Minimum combined score
    
    Returns:
        List of interactions involving the query protein
    """
    params = {
        "identifiers": protein,
        "species": species,
        "limit": limit,
        "required_score": min_score,
    }
    
    try:
        data = _string_request("json/interaction_partners", params, timeout=60)
        interactions = []
        for item in data:
            interactions.append(StringInteraction(
                protein_a=item.get("stringId_A", ""),
                protein_b=item.get("stringId_B", ""),
                gene_a=item.get("preferredName_A"),
                gene_b=item.get("preferredName_B"),
                score=int(item.get("score", 0) * 1000),
            ))
        return interactions
    except Exception as e:
        logger.error(f"STRING partners request failed: {e}")
        return []


def get_enrichment(
    proteins: list[str],
    species: int = DEFAULT_SPECIES,
) -> list[dict[str, Any]]:
    """
    Get functional enrichment for a protein set.
    
    Args:
        proteins: List of gene symbols
        species: NCBI taxonomy ID
    
    Returns:
        List of enriched terms (GO, KEGG, etc.)
    """
    if not proteins:
        return []
    
    identifiers = "\r".join(proteins)
    params = {
        "identifiers": identifiers,
        "species": species,
    }
    
    try:
        data = _string_request("json/enrichment", params, timeout=60)
        return data if isinstance(data, list) else []
    except Exception as e:
        logger.warning(f"STRING enrichment request failed: {e}")
        return []


def interactions_to_cytoscape(
    interactions: list[StringInteraction],
) -> tuple[list[dict], list[dict]]:
    """
    Convert STRING interactions to Cytoscape.js format.
    
    Returns:
        (nodes, edges) tuple for Cytoscape rendering
    """
    nodes_set = set()
    nodes = []
    edges = []
    
    for i in interactions:
        gene_a = i.gene_a or i.protein_a
        gene_b = i.gene_b or i.protein_b
        
        if gene_a not in nodes_set:
            nodes_set.add(gene_a)
            nodes.append({
                "data": {
                    "id": gene_a,
                    "label": gene_a,
                    "string_id": i.protein_a,
                    "node_type": "protein",
                }
            })
        
        if gene_b not in nodes_set:
            nodes_set.add(gene_b)
            nodes.append({
                "data": {
                    "id": gene_b,
                    "label": gene_b,
                    "string_id": i.protein_b,
                    "node_type": "protein",
                }
            })
        
        edges.append({
            "data": {
                "id": f"{gene_a}->{gene_b}",
                "source": gene_a,
                "target": gene_b,
                "score": i.score,
                "confidence": i.score / 1000.0,
            }
        })
    
    return nodes, edges


__all__ = [
    "StringInteraction",
    "StringProtein",
    "resolve_gene_symbols",
    "get_interactions",
    "get_interaction_partners",
    "get_enrichment",
    "interactions_to_cytoscape",
    "DEFAULT_SPECIES",
]
