"""
Data models for pathway analysis.

This module defines the core data structures used in pathway analysis:
- Pathway: Represents a biological pathway from KEGG or Reactome
- PathwayEnrichmentResult: Statistical results from enrichment analysis

These models are used throughout the pathway analysis pipeline for
mapping features to pathways and performing enrichment calculations.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Set


@dataclass
class Pathway:
    """
    Represents a biological pathway from KEGG or Reactome databases.

    Pathways are collections of genes, proteins, or metabolites that work together
    in a biological process. This class stores pathway metadata and the features
    (genes/proteins/metabolites) that are associated with it.

    Attributes:
        pathway_id: Unique pathway identifier
            - KEGG format: "hsa00010" (organism code + pathway number)
            - Reactome format: "R-HSA-71291" (stable identifier)
        name: Human-readable pathway name (e.g., "Glycolysis / Gluconeogenesis")
        source: Source database, either "KEGG" or "Reactome"
        description: Optional detailed pathway description
        features: Set of feature identifiers (gene symbols, protein IDs, metabolite IDs)
            that are mapped to this pathway
        feature_types: Set of omics types present in this pathway
            (e.g., {"gene", "protein"} or {"metabolite"})

    Example:
        >>> pathway = Pathway(
        ...     pathway_id="hsa00010",
        ...     name="Glycolysis / Gluconeogenesis",
        ...     source="KEGG",
        ...     description="Central metabolic pathway",
        ...     features={"ALDOA", "GAPDH", "PKM"},
        ...     feature_types={"gene", "protein"}
        ... )
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

    Contains statistical results from testing whether a set of input features
    is significantly enriched in a particular pathway compared to a background.

    Attributes:
        pathway: The Pathway object being tested
        input_features: Number of input features found in this pathway
        pathway_size: Total number of features in the pathway database
        background_size: Total number of features in the background/reference set
        p_value: Raw statistical p-value from Fisher's exact test
            (lower values indicate stronger enrichment)
        adjusted_p_value: FDR-adjusted p-value using Benjamini-Hochberg correction
            (accounts for multiple testing; use this for significance testing)
        enrichment_ratio: Ratio of observed to expected features
            - > 1.0: Enriched (more features than expected)
            - < 1.0: Depleted (fewer features than expected)
            - = 1.0: No enrichment
        matched_features: List of feature names from input that matched this pathway

    Example:
        >>> result = PathwayEnrichmentResult(
        ...     pathway=pathway,
        ...     input_features=15,
        ...     pathway_size=100,
        ...     background_size=20000,
        ...     p_value=0.001,
        ...     adjusted_p_value=0.01,
        ...     enrichment_ratio=2.5,
        ...     matched_features=["ALDOA", "GAPDH", "PKM"]
        ... )
        >>> result.adjusted_p_value < 0.05  # Significant enrichment
        True
    """

    pathway: Pathway
    input_features: int
    pathway_size: int
    background_size: int
    p_value: float
    adjusted_p_value: float
    enrichment_ratio: float
    matched_features: List[str] = field(default_factory=list)
