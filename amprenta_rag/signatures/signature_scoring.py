"""
Direction-aware, weight-aware signature scoring engine.

Computes signature scores based on:
- Species presence/absence
- Direction matching (↑/↓)
- Component weights
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set

from .signature_loader import Signature, SignatureComponent
from .species_matching import match_species, normalize_species_name


@dataclass
class ComponentMatch:
    """
    Match result for a single signature component.

    Attributes:
        signature_species: Original signature species name
        matched_dataset_species: Matched dataset species (or None)
        match_type: Type of match (exact, class_fallback, none)
        direction_match: Direction match status (match, conflict, neutral, unknown)
        weight: Component weight
    """

    signature_species: str
    matched_dataset_species: Optional[str] = None
    match_type: str = "none"  # exact, class_fallback, none
    direction_match: str = "unknown"  # match, conflict, neutral, unknown
    weight: float = 1.0


@dataclass
class SignatureScoreResult:
    """
    Complete signature scoring result.

    Attributes:
        total_score: Overall signature score (0-1)
        component_matches: List of component match results
        missing_species: Species in signature but not in dataset
        conflicting_species: Species with direction conflicts
        matched_species: Species successfully matched
    """

    total_score: float
    component_matches: List[ComponentMatch] = field(default_factory=list)
    missing_species: List[str] = field(default_factory=list)
    conflicting_species: List[str] = field(default_factory=list)
    matched_species: List[str] = field(default_factory=list)


def _get_direction_from_dataset(
    dataset_species: str,
    dataset_data: Dict[str, Dict[str, float]],
) -> Optional[str]:
    """
    Extract direction from dataset for a species.

    This is a placeholder for future direction inference from multi-condition datasets.
    For now, returns None (no direction info).

    Args:
        dataset_species: Species name
        dataset_data: Dataset data structure (future: condition → value mapping)

    Returns:
        Direction (↑, ↓, neutral) or None if unknown
    """
    # TODO: Implement direction inference from multi-condition data
    # For now, return None (no direction info available)
    return None


def score_signature(
    signature: Signature,
    dataset_species: Set[str],
    dataset_directions: Optional[Dict[str, str]] = None,
    refmet_map: Optional[Dict[str, str]] = None,
) -> SignatureScoreResult:
    """
    Score a signature against a dataset.

    Scoring logic:
    - match = +1.0    if dataset direction matches signature direction
    - match = -1.0    if opposite direction
    - match = +0.3    if present but no direction info
    - match = 0.0     if absent

    Final score = Σ (weight_i * match_i) / Σ(weight_i)

    Args:
        signature: Signature definition
        dataset_species: Set of species names in the dataset
        dataset_directions: Optional dict mapping species → direction (↑/↓)
        refmet_map: Optional RefMet mapping for species matching

    Returns:
        SignatureScoreResult with detailed scoring information
    """
    if dataset_directions is None:
        dataset_directions = {}

    # Match species
    signature_species_set = {comp.species for comp in signature.components}
    matches = match_species(
        dataset_species=dataset_species,
        signature_species=signature_species_set,
        refmet_map=refmet_map,
    )

    # Score each component
    component_matches: List[ComponentMatch] = []
    total_weighted_score = 0.0
    total_weight = 0.0

    for comp in signature.components:
        sig_species = comp.species
        matched_species = matches.get(sig_species)

        # Determine match type
        match_type = "none"
        if matched_species:
            if normalize_species_name(sig_species) == normalize_species_name(
                matched_species
            ):
                match_type = "exact"
            else:
                match_type = "class_fallback"

        # Determine direction match
        direction_match = "unknown"
        match_value = 0.0

        if matched_species:
            # Get dataset direction (if available)
            dataset_dir = dataset_directions.get(matched_species)
            sig_dir = comp.direction

            if sig_dir and dataset_dir:
                # Both have direction info
                if sig_dir == dataset_dir:
                    direction_match = "match"
                    match_value = 1.0
                elif (sig_dir == "↑" and dataset_dir == "↓") or (
                    sig_dir == "↓" and dataset_dir == "↑"
                ):
                    direction_match = "conflict"
                    match_value = -1.0
                else:
                    direction_match = "neutral"
                    match_value = 0.3
            elif sig_dir:
                # Only signature has direction (dataset present but no direction)
                direction_match = "neutral"
                match_value = 0.3
            else:
                # Species present but no direction info
                direction_match = "unknown"
                match_value = 0.3
        else:
            # Species not found
            direction_match = "missing"
            match_value = 0.0

        # Weight the match
        weight = comp.weight or 1.0
        weighted_score = weight * match_value
        total_weighted_score += weighted_score
        total_weight += weight

        component_matches.append(
            ComponentMatch(
                signature_species=sig_species,
                matched_dataset_species=matched_species,
                match_type=match_type,
                direction_match=direction_match,
                weight=weight,
            )
        )

    # Calculate final score
    if total_weight > 0:
        total_score = total_weighted_score / total_weight
        # Normalize to 0-1 range (since match_value can be -1 to +1)
        total_score = (total_score + 1.0) / 2.0
    else:
        total_score = 0.0

    # Categorize results
    missing_species = [
        comp.signature_species
        for comp in component_matches
        if comp.matched_dataset_species is None
    ]

    conflicting_species = [
        comp.signature_species
        for comp in component_matches
        if comp.direction_match == "conflict"
    ]

    matched_species = [
        comp.matched_dataset_species
        for comp in component_matches
        if comp.matched_dataset_species is not None
    ]

    return SignatureScoreResult(
        total_score=total_score,
        component_matches=component_matches,
        missing_species=missing_species,
        conflicting_species=conflicting_species,
        matched_species=matched_species,
    )
