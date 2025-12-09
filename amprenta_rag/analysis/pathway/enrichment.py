"""
Pathway enrichment analysis utilities.

This module provides statistical enrichment analysis to identify pathways
that are significantly over-represented in a set of input features compared
to a background set.

The analysis uses:
- Fisher's exact test for statistical significance
- Benjamini-Hochberg FDR correction for multiple testing
- Hypergeometric distribution for p-value calculation

Usage:
    >>> from amprenta_rag.analysis.pathway.enrichment import perform_pathway_enrichment
    >>> input_features = {"ALDOA", "GAPDH", "PKM", "ENO1"}
    >>> input_types = {"gene", "protein"}
    >>> results = perform_pathway_enrichment(
    ...     input_features=input_features,
    ...     input_feature_types=input_types,
    ...     p_value_threshold=0.05
    ... )
    >>> len(results) > 0
    True
"""

from __future__ import annotations

from typing import Dict, List, Optional, Set

from amprenta_rag.analysis.pathway.mapping import map_features_to_kegg_pathways, map_features_to_reactome_pathways
from amprenta_rag.analysis.pathway.models import Pathway, PathwayEnrichmentResult
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def perform_pathway_enrichment(
    input_features: Set[str],
    input_feature_types: Set[str],
    background_features: Optional[Set[str]] = None,
    pathway_sources: Optional[List[str]] = None,
    p_value_threshold: float = 0.05,
) -> List[PathwayEnrichmentResult]:
    """
    Perform pathway enrichment analysis.

    Tests whether input features are significantly enriched in biological pathways
    compared to a background set. Uses Fisher's exact test with FDR correction.

    The analysis:
    1. Maps input features to pathways (KEGG and/or Reactome)
    2. Calculates statistical significance for each pathway
    3. Applies multiple testing correction (Benjamini-Hochberg FDR)
    4. Filters results by p-value threshold
    5. Returns sorted list (most significant first)

    Args:
        input_features: Set of feature identifiers to test for enrichment
            (e.g., {"ALDOA", "GAPDH", "PKM"})
        input_feature_types: Set of omics types in input_features
            (e.g., {"gene", "protein"} or {"metabolite"})
        background_features: Optional background/reference set
            - If None: Uses pathway sizes as proxy for background
            - If provided: Uses this set as the reference population
        pathway_sources: List of pathway databases to query
            - None or ["KEGG", "Reactome"]: Query both databases
            - ["KEGG"]: Only KEGG pathways
            - ["Reactome"]: Only Reactome pathways (genes/proteins only)
        p_value_threshold: Adjusted p-value threshold for significance
            (default: 0.05, meaning 5% false discovery rate)

    Returns:
        List of PathwayEnrichmentResult objects, sorted by adjusted_p_value
        (most significant first). Only includes pathways with adjusted_p_value
        <= p_value_threshold.

    Note:
        - Reactome only supports genes and proteins (not metabolites)
        - Background size affects statistical power (larger = more conservative)
        - Multiple testing correction is critical when testing many pathways

    Example:
        >>> # Test if glycolysis genes are enriched
        >>> glycolysis_genes = {"ALDOA", "GAPDH", "PKM", "ENO1", "PFKM"}
        >>> results = perform_pathway_enrichment(
        ...     input_features=glycolysis_genes,
        ...     input_feature_types={"gene"},
        ...     pathway_sources=["KEGG"],
        ...     p_value_threshold=0.05
        ... )
        >>> # Should find glycolysis pathway as significantly enriched
        >>> any("glycolysis" in r.pathway.name.lower() for r in results)
        True
    """
    logger.info(
        "[ANALYSIS][PATHWAY] Performing pathway enrichment for %d features",
        len(input_features),
    )

    # Default to querying both KEGG and Reactome if not specified
    if pathway_sources is None:
        pathway_sources = ["KEGG", "Reactome"]

    # Step 1: Map input features to pathways
    # This queries KEGG and/or Reactome APIs to find which pathways contain our features
    all_pathways: Dict[str, Pathway] = {}

    for feature_type in input_feature_types:
        # Note: Currently using all input_features for each type
        # TODO: Filter features by type for more accurate mapping
        type_features = input_features

        # Query KEGG pathways (supports genes, proteins, metabolites)
        if "KEGG" in pathway_sources:
            kegg_pathways = map_features_to_kegg_pathways(type_features, feature_type)
            all_pathways.update(kegg_pathways)

        # Query Reactome pathways (only supports genes and proteins, not metabolites)
        if "Reactome" in pathway_sources and feature_type in ["gene", "protein"]:
            reactome_pathways = map_features_to_reactome_pathways(type_features, feature_type)
            all_pathways.update(reactome_pathways)

    # Early return if no pathways found
    if not all_pathways:
        logger.warning("[ANALYSIS][PATHWAY] No pathways found for input features")
        return []

    # Step 2: Perform statistical enrichment analysis
    enrichment_results: List[PathwayEnrichmentResult] = []

    # Determine background size for statistical testing
    # The background represents the total "universe" of features we're comparing against
    if background_features is None:
        # If no background provided, use the largest pathway size as a proxy
        # This is a heuristic - ideally should use actual background set
        background_size = (
            max(len(pathway.features) for pathway in all_pathways.values()) if all_pathways else len(input_features)
        )
    else:
        # Use provided background set size
        background_size = len(background_features)

    # Step 3: Calculate enrichment statistics for each pathway
    for pathway_id, pathway in all_pathways.items():
        # Count how many of our input features are in this pathway
        matched_features = [f for f in input_features if f in pathway.features]
        input_in_pathway = len(matched_features)
        pathway_size = len(pathway.features)

        # Skip pathways with no matches (can't be enriched)
        if input_in_pathway == 0:
            continue

        # Build 2x2 contingency table for Fisher's exact test
        # This table compares observed vs expected feature counts
        #
        #                 In Pathway  |  Not in Pathway
        # Input features |    a      |      b
        # Background    |    c      |      d
        #
        # Where:
        #   a = input features in pathway (observed)
        #   b = input features not in pathway
        #   c = background features in pathway (excluding input)
        #   d = background features not in pathway (excluding input)
        a = input_in_pathway
        b = len(input_features) - input_in_pathway
        c = pathway_size - input_in_pathway
        d = background_size - len(input_features) - c

        # Skip if contingency table is invalid (can't compute statistics)
        if a + c == 0 or b + d == 0:
            continue

        # Step 4: Calculate statistical significance using Fisher's exact test
        # Fisher's exact test uses the hypergeometric distribution to test if
        # the observed overlap (a) is significantly greater than expected by chance
        try:
            from scipy.stats import fisher_exact

            # Perform one-sided test (alternative="greater")
            # Tests if input features are OVER-represented in pathway
            odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")

            # Calculate enrichment ratio (observed / expected)
            # This tells us how many times more features we see than expected
            # - ratio > 1.0: Enriched (more than expected)
            # - ratio = 1.0: No enrichment
            # - ratio < 1.0: Depleted (fewer than expected)
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
        except ImportError:
            logger.warning("[ANALYSIS][PATHWAY] scipy not available, using simplified p-value calculation")
            # Simplified p-value approximation
            p_value = _calculate_simplified_p_value(a, b, c, d)
            (a * d) / (b * c) if (b * c) > 0 else 0.0

            if p_value < 1.0:  # Only add if significant
                expected = (len(input_features) * pathway_size) / background_size if background_size > 0 else 0
                enrichment_ratio = input_in_pathway / expected if expected > 0 else 0

                enrichment_results.append(
                    PathwayEnrichmentResult(
                        pathway=pathway,
                        input_features=input_in_pathway,
                        pathway_size=pathway_size,
                        background_size=background_size,
                        p_value=p_value,
                        adjusted_p_value=p_value,
                        enrichment_ratio=enrichment_ratio,
                        matched_features=matched_features,
                    )
                )
        except (ValueError, ZeroDivisionError, AttributeError) as e:
            logger.debug(
                "[ANALYSIS][PATHWAY] Error computing enrichment for pathway %s: %r",
                pathway_id,
                e,
            )
            continue

    # Step 5: Apply multiple testing correction
    # When testing many pathways, we need to correct for multiple comparisons
    # to control the false discovery rate (FDR)
    if enrichment_results:
        enrichment_results = _apply_multiple_testing_correction(enrichment_results)

    # Step 6: Filter and sort results
    # Keep only pathways that meet the significance threshold
    # Sort by adjusted p-value (most significant first)
    significant_results = [r for r in enrichment_results if r.adjusted_p_value <= p_value_threshold]
    significant_results.sort(key=lambda r: r.adjusted_p_value)

    logger.info(
        "[ANALYSIS][PATHWAY] Found %d significantly enriched pathways (p < %.3f)",
        len(significant_results),
        p_value_threshold,
    )

    return significant_results


def _calculate_simplified_p_value(a: int, b: int, c: int, d: int) -> float:
    """
    Calculate simplified p-value approximation when scipy is not available.

    This is a fallback method that provides a rough approximation of statistical
    significance. For accurate results, scipy.stats.fisher_exact should be used.

    The approximation compares observed vs expected ratios:
    - observed_ratio = proportion of input features in pathway
    - expected_ratio = proportion of all features in pathway
    - p-value ≈ difference between ratios

    Args:
        a: Input features in pathway (observed count)
        b: Input features not in pathway
        c: Background features in pathway (excluding input)
        d: Background features not in pathway (excluding input)

    Returns:
        Approximate p-value (0.0 to 1.0)
        - Lower values indicate stronger enrichment
        - This is a heuristic, not a true statistical test

    Note:
        This is a simplified approximation. For production use, ensure scipy
        is installed for accurate Fisher's exact test calculations.
    """
    if a + b + c + d == 0:
        return 1.0

    # Calculate observed ratio: what proportion of input features are in pathway?
    observed_ratio = a / (a + b) if (a + b) > 0 else 0.0

    # Calculate expected ratio: what proportion of all features are in pathway?
    expected_ratio = (a + c) / (a + b + c + d) if (a + b + c + d) > 0 else 0.0

    if expected_ratio == 0:
        return 1.0

    # Approximate p-value as normalized difference from expected
    # This is a heuristic - not a true statistical test
    p_value = min(1.0, abs(observed_ratio - expected_ratio) / expected_ratio)
    return p_value


def _apply_multiple_testing_correction(
    enrichment_results: List[PathwayEnrichmentResult],
) -> List[PathwayEnrichmentResult]:
    """
    Apply multiple testing correction to enrichment results.

    When testing multiple pathways simultaneously, we need to correct for
    multiple comparisons to control the false discovery rate (FDR).

    Methods:
    - Primary: Benjamini-Hochberg FDR correction (less conservative)
    - Fallback: Bonferroni correction (more conservative, if statsmodels unavailable)

    Args:
        enrichment_results: List of enrichment results with raw p-values

    Returns:
        List of enrichment results with adjusted_p_value field updated

    Note:
        - Benjamini-Hochberg is preferred as it's less conservative than Bonferroni
        - Bonferroni multiplies each p-value by number of tests (very conservative)
        - Adjusted p-values should be used for significance testing, not raw p-values
    """
    try:
        from statsmodels.stats.multitest import multipletests

        # Extract raw p-values for correction
        p_values = [r.p_value for r in enrichment_results]

        # Apply Benjamini-Hochberg FDR correction
        # This method controls the false discovery rate (proportion of false positives)
        # Less conservative than Bonferroni, more appropriate for exploratory analysis
        _, adjusted_p_values, _, _ = multipletests(
            p_values,
            method="fdr_bh",  # Benjamini-Hochberg method
        )

        # Update each result with its adjusted p-value
        for i, result in enumerate(enrichment_results):
            result.adjusted_p_value = adjusted_p_values[i]
    except ImportError:
        logger.warning("[ANALYSIS][PATHWAY] statsmodels not available, using Bonferroni correction")
        # Fallback: Use Bonferroni correction (multiply by number of tests)
        # This is more conservative but doesn't require external dependencies
        # Bonferroni: adjusted_p = p_value * n_tests
        n_tests = len(enrichment_results)
        for result in enrichment_results:
            result.adjusted_p_value = min(1.0, result.p_value * n_tests)

    return enrichment_results
