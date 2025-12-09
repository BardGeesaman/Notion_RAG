"""
Cross-Omics Pathway Analysis.

Maps features (genes, proteins, metabolites, lipids) to biological pathways
using KEGG and Reactome databases, and performs pathway enrichment analysis.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set


from amprenta_rag.analysis.id_mapping import (
    map_gene_to_kegg,
    map_gene_to_reactome,
    map_metabolite_to_kegg,
    map_protein_to_kegg,
    map_protein_to_reactome,
    map_protein_to_uniprot,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


@dataclass
class Pathway:
    """
    Represents a biological pathway.
    
    Attributes:
        pathway_id: Pathway identifier (e.g., "hsa00010" for KEGG, "R-HSA-71291" for Reactome)
        name: Pathway name
        source: Source database ("KEGG" or "Reactome")
        description: Pathway description
        features: Set of feature identifiers mapped to this pathway
        feature_types: Set of omics types (gene, protein, metabolite, lipid)
    """
    pathway_id: str
    name: str
    source: str  # "KEGG" or "Reactome"
    description: Optional[str] = None
    features: Set[str] = field(default_factory=set)
    feature_types: Set[str] = field(default_factory=set)


@dataclass
class PathwayEnrichmentResult:
    """
    Results from pathway enrichment analysis.
    
    Attributes:
        pathway: Pathway object
        input_features: Number of input features in pathway
        pathway_size: Total number of features in pathway
        background_size: Total number of features in background
        p_value: Statistical p-value (Fisher's exact test or hypergeometric)
        adjusted_p_value: FDR-adjusted p-value (Benjamini-Hochberg)
        enrichment_ratio: Ratio of observed/expected
        matched_features: List of feature names that matched
    """
    pathway: Pathway
    input_features: int
    pathway_size: int
    background_size: int
    p_value: float
    adjusted_p_value: float
    enrichment_ratio: float
    matched_features: List[str] = field(default_factory=list)


def map_features_to_kegg_pathways(
    features: Set[str],
    feature_type: str,
) -> Dict[str, Pathway]:
    """
    Map features to KEGG pathways.
    
    Args:
        features: Set of feature identifiers
        feature_type: One of "gene", "protein", "metabolite"
        
    Returns:
        Dictionary mapping pathway_id -> Pathway object
    """
    logger.info(
        "[ANALYSIS][PATHWAY] Mapping %d %s features to KEGG pathways",
        len(features),
        feature_type,
    )
    
    # KEGG API integration
    # For now, return empty dict - will implement actual API calls
    # KEGG REST API: https://www.kegg.jp/kegg/rest/keggapi.html
    
    pathways: Dict[str, Pathway] = {}
    
    try:
        import requests
        
        # KEGG REST API base URL
        KEGG_BASE_URL = "https://rest.kegg.jp"
        
        # Map feature type to KEGG database
        db_map = {
            "gene": "hsa",  # Human genes (KEGG organism code)
            "protein": "hsa",
            "metabolite": "cpd",  # KEGG compound database
        }
        
        if feature_type not in db_map:
            logger.warning(
                "[ANALYSIS][PATHWAY] KEGG mapping not supported for feature_type '%s'",
                feature_type,
            )
            return pathways
        
        db_map[feature_type]
        
        # For genes/proteins: use KEGG gene-to-pathway mapping
        # For metabolites: use KEGG compound-to-pathway mapping
        
        # Batch query KEGG (rate limit: ~1 request/second)
        for feature in features:
            try:
                # Convert feature to KEGG ID format using ID mapping service
                if feature_type == "gene":
                    kegg_id = map_gene_to_kegg(feature, organism="hsa")
                elif feature_type == "protein":
                    kegg_id = map_protein_to_kegg(feature, organism="hsa")
                elif feature_type == "metabolite":
                    kegg_id = map_metabolite_to_kegg(feature)
                else:
                    kegg_id = None
                
                if not kegg_id:
                    continue
                
                # Query pathway mapping
                url = f"{KEGG_BASE_URL}/link/pathway/{kegg_id}"
                resp = requests.get(url, timeout=10)
                
                if resp.status_code == 200:
                    # Parse response (format: pathway_id\tkegg_id)
                    for line in resp.text.strip().split("\n"):
                        if not line:
                            continue
                        parts = line.split("\t")
                        if len(parts) >= 2:
                            pathway_id = parts[0]
                            if pathway_id not in pathways:
                                # Fetch pathway info
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
                            if pathway_id in pathways:
                                pathways[pathway_id].features.add(feature)
                
                # Rate limiting
                import time
                time.sleep(0.5)  # KEGG rate limit
                
            except Exception as e:
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
            "[ANALYSIS][PATHWAY] 'requests' not available for KEGG API calls"
        )
    except Exception as e:
        logger.error(
            "[ANALYSIS][PATHWAY] Error in KEGG pathway mapping: %r",
            e,
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
        import requests
        
        # Reactome REST API base URL
        REACTOME_BASE_URL = "https://reactome.org/ContentService"
        
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
                import time
                time.sleep(0.1)  # Reactome rate limit
                
            except Exception as e:
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
            "[ANALYSIS][PATHWAY] 'requests' not available for Reactome API calls"
        )
    except Exception as e:
        logger.error(
            "[ANALYSIS][PATHWAY] Error in Reactome pathway mapping: %r",
            e,
        )
    
    return pathways


def _fetch_kegg_pathway_info(pathway_id: str) -> Optional[Dict[str, str]]:
    """
    Fetch pathway information from KEGG.
    """
    try:
        import requests
        
        url = f"https://rest.kegg.jp/get/{pathway_id}"
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
    except Exception:
        pass
    return None


def perform_pathway_enrichment(
    input_features: Set[str],
    input_feature_types: Set[str],
    background_features: Optional[Set[str]] = None,
    pathway_sources: List[str] = None,
    p_value_threshold: float = 0.05,
) -> List[PathwayEnrichmentResult]:
    """
    Perform pathway enrichment analysis.
    
    Args:
        input_features: Set of feature identifiers to test
        input_feature_types: Set of omics types (gene, protein, metabolite, lipid)
        background_features: Optional background set (if None, uses all known features)
        pathway_sources: List of sources to use (["KEGG", "Reactome"] or None for all)
        p_value_threshold: P-value threshold for significance
        
    Returns:
        List of PathwayEnrichmentResult objects, sorted by p-value
    """
    logger.info(
        "[ANALYSIS][PATHWAY] Performing pathway enrichment for %d features",
        len(input_features),
    )
    
    if pathway_sources is None:
        pathway_sources = ["KEGG", "Reactome"]
    
    # Map features to pathways
    all_pathways: Dict[str, Pathway] = {}
    
    for feature_type in input_feature_types:
        # Get features of this type
        type_features = input_features  # Simplified - would need to filter by type
        
        if "KEGG" in pathway_sources:
            kegg_pathways = map_features_to_kegg_pathways(type_features, feature_type)
            all_pathways.update(kegg_pathways)
        
        if "Reactome" in pathway_sources and feature_type in ["gene", "protein"]:
            reactome_pathways = map_features_to_reactome_pathways(type_features, feature_type)
            all_pathways.update(reactome_pathways)
    
    if not all_pathways:
        logger.warning(
            "[ANALYSIS][PATHWAY] No pathways found for input features"
        )
        return []
    
    # Perform enrichment analysis
    enrichment_results: List[PathwayEnrichmentResult] = []
    
    # Background size (total features in database)
    if background_features is None:
        # Use pathway sizes as proxy for background
        background_size = max(
            len(pathway.features) for pathway in all_pathways.values()
        ) if all_pathways else len(input_features)
    else:
        background_size = len(background_features)
    
    for pathway_id, pathway in all_pathways.items():
        # Count input features in pathway
        matched_features = [f for f in input_features if f in pathway.features]
        input_in_pathway = len(matched_features)
        pathway_size = len(pathway.features)
        
        if input_in_pathway == 0:
            continue
        
        # Contingency table:
        #                 In Pathway  |  Not in Pathway
        # Input features |    a      |      b
        # Background    |    c      |      d
        a = input_in_pathway
        b = len(input_features) - input_in_pathway
        c = pathway_size - input_in_pathway
        d = background_size - len(input_features) - c
        
        # Fisher's exact test (hypergeometric distribution)
        try:
            from scipy.stats import fisher_exact
        except ImportError:
            logger.warning(
                "[ANALYSIS][PATHWAY] scipy not available, using simplified p-value calculation"
            )
            # Simplified p-value approximation
            p_value = 1.0
            odds_ratio = 0.0
            if a + c > 0 and b + d > 0:
                # Approximate using hypergeometric distribution
                pass
                # Simplified: use ratio-based approximation
                p_value = min(1.0, (a / (a + b)) / ((a + c) / (a + b + c + d)) if (a + b + c + d) > 0 else 1.0)
                odds_ratio = (a * d) / (b * c) if (b * c) > 0 else 0.0
        
        if a + c == 0 or b + d == 0:
            continue
        
        try:
            odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
            
            # Enrichment ratio
            expected = (len(input_features) * pathway_size) / background_size if background_size > 0 else 0
            enrichment_ratio = input_in_pathway / expected if expected > 0 else 0
            
            enrichment_results.append(
                PathwayEnrichmentResult(
                    pathway=pathway,
                    input_features=input_in_pathway,
                    pathway_size=pathway_size,
                    background_size=background_size,
                    p_value=p_value,
                    adjusted_p_value=p_value,  # Will adjust below
                    enrichment_ratio=enrichment_ratio,
                    matched_features=matched_features,
                )
            )
        except Exception as e:
            logger.debug(
                "[ANALYSIS][PATHWAY] Error computing enrichment for pathway %s: %r",
                pathway_id,
                e,
            )
            continue
    
    # Multiple testing correction (Benjamini-Hochberg FDR)
    if enrichment_results:
        try:
            from statsmodels.stats.multitest import multipletests
            
            p_values = [r.p_value for r in enrichment_results]
            _, adjusted_p_values, _, _ = multipletests(
                p_values,
                method="fdr_bh",
            )
            
            for i, result in enumerate(enrichment_results):
                result.adjusted_p_value = adjusted_p_values[i]
        except ImportError:
            logger.warning(
                "[ANALYSIS][PATHWAY] statsmodels not available, skipping FDR correction"
            )
            # Use Bonferroni correction as fallback
            n_tests = len(enrichment_results)
            for result in enrichment_results:
                result.adjusted_p_value = min(1.0, result.p_value * n_tests)
    
    # Filter by threshold and sort
    significant_results = [
        r for r in enrichment_results
        if r.adjusted_p_value <= p_value_threshold
    ]
    significant_results.sort(key=lambda r: r.adjusted_p_value)
    
    logger.info(
        "[ANALYSIS][PATHWAY] Found %d significantly enriched pathways (p < %.3f)",
        len(significant_results),
        p_value_threshold,
    )
    
    return significant_results

