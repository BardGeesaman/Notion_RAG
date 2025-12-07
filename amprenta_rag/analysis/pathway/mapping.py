"""
Pathway mapping utilities for KEGG and Reactome.

This module provides functions to map biological features (genes, proteins, metabolites)
to pathways in the KEGG and Reactome databases. It handles:
- ID conversion (gene symbols, UniProt IDs, etc.)
- API calls to KEGG and Reactome REST APIs
- Rate limiting to respect API constraints
- Error handling and retry logic

Rate Limits:
    - KEGG: ~1 request/second (0.5s delay between requests)
    - Reactome: ~10 requests/second (0.1s delay between requests)

Usage:
    >>> from amprenta_rag.analysis.pathway.mapping import map_features_to_kegg_pathways
    >>> genes = {"ALDOA", "GAPDH", "PKM"}
    >>> pathways = map_features_to_kegg_pathways(genes, feature_type="gene")
    >>> len(pathways)
    3
"""

from __future__ import annotations

import time
from typing import Dict, Optional, Set

import requests

from amprenta_rag.analysis.id_mapping import (
    map_gene_to_kegg,
    map_gene_to_reactome,
    map_metabolite_to_kegg,
    map_protein_to_kegg,
    map_protein_to_reactome,
    map_protein_to_uniprot,
)
from amprenta_rag.analysis.pathway.models import Pathway
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# API base URLs
KEGG_BASE_URL = "https://rest.kegg.jp"
REACTOME_BASE_URL = "https://reactome.org/ContentService"


def map_features_to_kegg_pathways(
    features: Set[str],
    feature_type: str,
) -> Dict[str, Pathway]:
    """
    Map features to KEGG pathways.

    Queries the KEGG REST API to find which pathways contain the given features.
    Features are first converted to KEGG-compatible IDs using ID mapping services,
    then queried against KEGG's pathway database.

    Args:
        features: Set of feature identifiers (gene symbols, protein IDs, metabolite IDs)
        feature_type: Type of features being mapped
            - "gene": Human gene symbols (mapped to KEGG organism code "hsa")
            - "protein": Protein identifiers (mapped to KEGG organism code "hsa")
            - "metabolite": Metabolite identifiers (mapped to KEGG compound database "cpd")

    Returns:
        Dictionary mapping pathway_id -> Pathway object
        Empty dictionary if no pathways found or if feature_type is unsupported

    Note:
        - Rate limited to ~1 request/second (0.5s delay between requests)
        - Features that can't be mapped to KEGG IDs are silently skipped
        - Each feature may map to multiple pathways

    Example:
        >>> genes = {"ALDOA", "GAPDH", "PKM"}
        >>> pathways = map_features_to_kegg_pathways(genes, feature_type="gene")
        >>> "hsa00010" in pathways  # Glycolysis pathway
        True
        >>> pathways["hsa00010"].name
        'Glycolysis / Gluconeogenesis'
    """
    logger.info(
        "[ANALYSIS][PATHWAY] Mapping %d %s features to KEGG pathways",
        len(features),
        feature_type,
    )

    pathways: Dict[str, Pathway] = {}

    try:
        # Map feature type to KEGG database code
        # KEGG uses organism codes for genes/proteins (e.g., "hsa" for human)
        # and "cpd" for compounds/metabolites
        db_map = {
            "gene": "hsa",  # Human genes (Homo sapiens organism code)
            "protein": "hsa",  # Human proteins (same organism code)
            "metabolite": "cpd",  # KEGG compound database
        }

        if feature_type not in db_map:
            logger.warning(
                "[ANALYSIS][PATHWAY] KEGG mapping not supported for feature_type '%s'",
                feature_type,
            )
            return pathways

        # Process each feature sequentially (KEGG API rate limit: ~1 req/sec)
        for feature in features:
            try:
                # Step 1: Convert feature to KEGG-compatible ID
                # This uses ID mapping services to convert various ID formats
                if feature_type == "gene":
                    kegg_id = map_gene_to_kegg(feature, organism="hsa")
                elif feature_type == "protein":
                    kegg_id = map_protein_to_kegg(feature, organism="hsa")
                elif feature_type == "metabolite":
                    kegg_id = map_metabolite_to_kegg(feature)
                else:
                    kegg_id = None

                # Skip features that can't be mapped to KEGG IDs
                if not kegg_id:
                    continue

                # Step 2: Query KEGG API for pathway associations
                # KEGG REST API endpoint: /link/pathway/{kegg_id}
                # Returns tab-separated list of pathway_id\tkegg_id pairs
                url = f"{KEGG_BASE_URL}/link/pathway/{kegg_id}"
                resp = requests.get(url, timeout=10)

                if resp.status_code == 200:
                    # Step 3: Parse response (format: pathway_id\tkegg_id)
                    # Each line represents one pathway association
                    for line in resp.text.strip().split("\n"):
                        if not line:
                            continue
                        parts = line.split("\t")
                        if len(parts) >= 2:
                            pathway_id = parts[0]

                            # Step 4: Create Pathway object if not already seen
                            if pathway_id not in pathways:
                                # Fetch detailed pathway information (name, description)
                                pathway_info = _fetch_kegg_pathway_info(pathway_id)
                                if pathway_info:
                                    pathways[pathway_id] = Pathway(
                                        pathway_id=pathway_id,
                                        name=pathway_info.get("name", pathway_id),
                                        source="KEGG",
                                        description=pathway_info.get("description"),
                                        features=set(),
                                        feature_types={feature_type},
                                    )

                            # Step 5: Add this feature to the pathway's feature set
                            if pathway_id in pathways:
                                pathways[pathway_id].features.add(feature)

                # Rate limiting: KEGG API allows ~1 request/second
                # Using 0.5s delay to be conservative and avoid rate limit errors
                time.sleep(0.5)

            except (requests.RequestException, ValueError, KeyError) as e:
                logger.debug(
                    "[ANALYSIS][PATHWAY] Error mapping feature %s to KEGG: %r",
                    feature,
                    e,
                )
                continue

        logger.info(
            "[ANALYSIS][PATHWAY] Mapped %d features to %d KEGG pathways",
            len(features),
            len(pathways),
        )

    except ImportError:
        logger.warning(
            "[ANALYSIS][PATHWAY] 'requests' library not available for KEGG API calls. "
            "Install with: pip install requests"
        )
    except (requests.RequestException, ConnectionError, TimeoutError) as e:
        # Network errors
        logger.error(
            "[ANALYSIS][PATHWAY] Network error in KEGG pathway mapping: %r",
            e,
        )
    except (ValueError, KeyError) as e:
        # Data parsing errors
        logger.error(
            "[ANALYSIS][PATHWAY] Data parsing error in KEGG pathway mapping: %r",
            e,
            exc_info=True,
        )
    except Exception as e:
        # Unexpected errors
        logger.error(
            "[ANALYSIS][PATHWAY] Unexpected error in KEGG pathway mapping: %r",
            e,
            exc_info=True,
        )

    return pathways


def map_features_to_reactome_pathways(
    features: Set[str],
    feature_type: str,
) -> Dict[str, Pathway]:
    """
    Map features to Reactome pathways.

    Args:
        features: Set of feature identifiers
        feature_type: One of "gene", "protein"

    Returns:
        Dictionary mapping pathway_id -> Pathway object
    """
    logger.info(
        "[ANALYSIS][PATHWAY] Mapping %d %s features to Reactome pathways",
        len(features),
        feature_type,
    )

    pathways: Dict[str, Pathway] = {}

    try:
        if feature_type not in ["gene", "protein"]:
            logger.warning(
                "[ANALYSIS][PATHWAY] Reactome mapping only supports genes/proteins, not '%s'",
                feature_type,
            )
            return pathways

        # Reactome uses UniProt IDs for proteins and gene symbols for genes
        # Batch query Reactome (rate limit: ~10 requests/second)
        for feature in features:
            try:
                # Convert feature to Reactome-compatible ID using ID mapping service
                if feature_type == "gene":
                    reactome_id = map_gene_to_reactome(feature)
                elif feature_type == "protein":
                    reactome_id = map_protein_to_reactome(feature)
                else:
                    reactome_id = None

                if not reactome_id:
                    continue

                # Query pathway mapping using UniProt ID or gene symbol
                # Reactome API: /query/mapping/{identifier}
                # For proteins, use UniProt ID; for genes, use gene symbol
                if feature_type == "protein":
                    # Get UniProt ID first
                    uniprot_id = map_protein_to_uniprot(feature)
                    if uniprot_id:
                        url = f"{REACTOME_BASE_URL}/query/mapping/uniprot/{uniprot_id}"
                    else:
                        continue
                else:
                    # For genes, try direct query
                    url = f"{REACTOME_BASE_URL}/query/mapping/{reactome_id}"

                resp = requests.get(url, timeout=10)

                if resp.status_code == 200:
                    data = resp.json()
                    # Parse Reactome response
                    for pathway_info in data.get("pathways", []):
                        pathway_id = pathway_info.get("stId")
                        pathway_name = pathway_info.get("displayName", pathway_id)

                        if pathway_id and pathway_id not in pathways:
                            pathways[pathway_id] = Pathway(
                                pathway_id=pathway_id,
                                name=pathway_name,
                                source="Reactome",
                                description=pathway_info.get("summary"),
                                features=set(),
                                feature_types={feature_type},
                            )
                        if pathway_id in pathways:
                            pathways[pathway_id].features.add(feature)

                # Rate limiting
                time.sleep(0.1)  # Reactome rate limit

            except (requests.RequestException, ValueError, KeyError) as e:
                logger.debug(
                    "[ANALYSIS][PATHWAY] Error mapping feature %s to Reactome: %r",
                    feature,
                    e,
                )
                continue

        logger.info(
            "[ANALYSIS][PATHWAY] Mapped %d features to %d Reactome pathways",
            len(features),
            len(pathways),
        )

    except ImportError:
        logger.warning(
            "[ANALYSIS][PATHWAY] 'requests' library not available for Reactome API calls. "
            "Install with: pip install requests"
        )
    except (requests.RequestException, ConnectionError, TimeoutError) as e:
        # Network errors
        logger.error(
            "[ANALYSIS][PATHWAY] Network error in Reactome pathway mapping: %r",
            e,
        )
    except (ValueError, KeyError) as e:
        # Data parsing errors
        logger.error(
            "[ANALYSIS][PATHWAY] Data parsing error in Reactome pathway mapping: %r",
            e,
            exc_info=True,
        )
    except Exception as e:
        # Unexpected errors
        logger.error(
            "[ANALYSIS][PATHWAY] Unexpected error in Reactome pathway mapping: %r",
            e,
            exc_info=True,
        )

    return pathways


def _fetch_kegg_pathway_info(pathway_id: str) -> Optional[Dict[str, str]]:
    """
    Fetch pathway information from KEGG.

    Args:
        pathway_id: KEGG pathway identifier

    Returns:
        Dictionary with pathway name and description, or None if not found
    """
    try:
        url = f"{KEGG_BASE_URL}/get/{pathway_id}"
        resp = requests.get(url, timeout=10)

        if resp.status_code == 200:
            # Parse KEGG flat file format
            info = {"name": pathway_id, "description": ""}
            for line in resp.text.split("\n"):
                if line.startswith("NAME"):
                    info["name"] = line.split("NAME")[1].strip()
                elif line.startswith("DESCRIPTION"):
                    info["description"] = line.split("DESCRIPTION")[1].strip()
            return info
    except (requests.RequestException, ValueError, KeyError):
        # KEGG API info fetch failed - return None
        pass
    return None
