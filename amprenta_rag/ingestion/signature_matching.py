# amprenta_rag/ingestion/signature_matching.py

"""
Automatic signature matching and scoring for datasets.

This module provides:
- Signature scoring against datasets
- Automatic signature detection during dataset ingestion
- Notion writebacks for signature matches
- Lipidomics-aware species mapping
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

import requests

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config
from amprenta_rag.ingestion.feature_extraction import \
    extract_features_from_mwtab
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.signatures.signature_loader import (Signature,
                                                      SignatureComponent,
                                                      load_signature_from_tsv)
from amprenta_rag.signatures.signature_scoring import (SignatureScoreResult,
                                                       score_signature)
from amprenta_rag.signatures.species_matching import (match_species,
                                                      normalize_species_name)

logger = get_logger(__name__)


@dataclass
class SignatureMatchResult:
    """Result of matching a signature against a dataset."""

    signature_page_id: str
    signature_name: str
    score: float
    overlap_fraction: float
    matched_components: List[str]
    missing_components: List[str]
    conflicting_components: List[str]
    score_result: SignatureScoreResult


def map_raw_lipid_to_canonical_species(raw_name: str) -> Optional[str]:
    """
    Map raw lipid name from dataset to canonical species name.

    Uses normalization and class-level matching to handle vendor formats.

    Examples:
        "CER 16:0" → "Cer(d18:1/16:0)"
        "Ceramide C24:1" → "Cer(d18:1/24:1)"
        "SM 18:0" → "SM(d18:1/18:0)"

    Args:
        raw_name: Raw lipid name from dataset

    Returns:
        Canonical species name or None if cannot be mapped
    """
    if not raw_name:
        return None

    normalized = normalize_species_name(raw_name)

    # Try to reconstruct canonical format from normalized name
    # This is a simplified mapping - can be enhanced with Lipid Species DB lookup
    raw_lower = raw_name.lower().strip()

    # Handle common vendor formats
    # CER 16:0 → Cer(d18:1/16:0)
    cer_match = re.match(r"cer\s*([cC]?)(\d+):(\d+)", raw_lower)
    if cer_match:
        chain = f"{cer_match.group(2)}:{cer_match.group(3)}"
        return f"Cer(d18:1/{chain})"

    # SM 16:0 → SM(d18:1/16:0)
    sm_match = re.match(r"sm\s*([cC]?)(\d+):(\d+)", raw_lower)
    if sm_match:
        chain = f"{sm_match.group(2)}:{sm_match.group(3)}"
        return f"SM(d18:1/{chain})"

    # If already in canonical format, return as-is
    if "(" in raw_name and ")" in raw_name:
        return raw_name

    return None


def fetch_all_signatures_from_notion() -> List[Dict[str, Any]]:
    """
    Fetch all signature pages from Notion Lipid Signatures database.

    Returns:
        List of signature page dictionaries
    """
    cfg = get_config()

    if not cfg.notion.signature_db_id:
        logger.warning("[INGEST][SIGNATURE-MATCH] Signature database ID not configured")
        return []

    try:
        url = f"{cfg.notion.base_url}/databases/{cfg.notion.signature_db_id}/query"
        all_pages = []
        has_more = True
        start_cursor = None

        while has_more:
            payload = {
                "page_size": 100,
            }
            if start_cursor:
                payload["start_cursor"] = start_cursor

            resp = requests.post(
                url,
                headers=notion_headers(),
                json=payload,
                timeout=30,
            )
            resp.raise_for_status()

            data = resp.json()
            all_pages.extend(data.get("results", []))
            has_more = data.get("has_more", False)
            start_cursor = data.get("next_cursor")

        logger.debug(
            "[INGEST][SIGNATURE-MATCH] Fetched %d signatures from Notion",
            len(all_pages),
        )
        return all_pages

    except Exception as e:
        logger.warning(
            "[INGEST][SIGNATURE-MATCH] Error fetching signatures from Notion: %r",
            e,
        )
        return []


def load_signature_from_notion_page(
    signature_page: Dict[str, Any],
) -> Optional[Signature]:
    """
    Load a Signature object from a Notion page.

    Fetches components from the signature page and builds a Signature object.

    Args:
        signature_page: Notion page dictionary for a signature

    Returns:
        Signature object or None if loading failed
    """
    cfg = get_config()
    signature_page_id = signature_page.get("id", "")

    if not signature_page_id:
        logger.warning(
            "[INGEST][SIGNATURE-MATCH] Signature page has no ID",
        )
        return None

    # Get signature name
    props = signature_page.get("properties", {}) or {}
    name_prop = props.get("Name", {}).get("title", []) or []
    signature_name = name_prop[0].get("plain_text", "") if name_prop else "Unknown"

    # Get description if available
    desc_prop = props.get("Description", {}) or {}
    if desc_prop.get("rich_text"):
        description = desc_prop["rich_text"][0].get("plain_text", "")
    else:
        description = None

    # Fetch components from Signature Components DB
    if (
        not hasattr(cfg.notion, "signature_component_db_id")
        or not cfg.notion.signature_component_db_id
    ):
        logger.debug(
            "[INGEST][SIGNATURE-MATCH] Component DB ID not configured",
        )
        return None

    components: List[SignatureComponent] = []

    try:
        # Query for components linked to this signature
        url = f"{cfg.notion.base_url}/databases/{cfg.notion.signature_component_db_id}/query"
        all_components = []
        has_more = True
        start_cursor = None

        while has_more:
            payload = {
                "filter": {
                    "property": "Signature",
                    "relation": {"contains": signature_page_id},
                },
                "page_size": 100,
            }
            if start_cursor:
                payload["start_cursor"] = start_cursor

            resp = requests.post(
                url,
                headers=notion_headers(),
                json=payload,
                timeout=30,
            )
            resp.raise_for_status()

            data = resp.json()
            all_components.extend(data.get("results", []))
            has_more = data.get("has_more", False)
            start_cursor = data.get("next_cursor")

        # Parse components
        for comp_page in all_components:
            comp_props = comp_page.get("properties", {}) or {}

            # Get species name (Component Name)
            name_prop = comp_props.get("Component Name", {}).get("title", []) or []
            species = name_prop[0].get("plain_text", "") if name_prop else ""

            if not species:
                continue

            # Get direction
            direction_prop = comp_props.get("Direction", {}).get("select")
            direction = None
            if direction_prop:
                direction_str = direction_prop.get("name", "")
                # Map Notion values to signature format
                direction_map = {
                    "Up": "↑",
                    "Down": "↓",
                    "NoChange": "neutral",
                    "Complex": "complex",
                    "Unknown": None,
                }
                direction = direction_map.get(direction_str, None)

            # Get weight
            weight_prop = comp_props.get("Weight", {}).get("number")
            weight = float(weight_prop) if weight_prop is not None else None

            components.append(
                SignatureComponent(
                    species=species,
                    direction=direction,
                    weight=weight,
                )
            )

        if not components:
            logger.debug(
                "[INGEST][SIGNATURE-MATCH] No components found for signature %s",
                signature_name,
            )
            return None

        signature = Signature(
            name=signature_name,
            components=components,
            description=description,
        )

        logger.debug(
            "[INGEST][SIGNATURE-MATCH] Loaded signature '%s' with %d components",
            signature_name,
            len(components),
        )

        return signature

    except Exception as e:
        logger.warning(
            "[INGEST][SIGNATURE-MATCH] Error loading signature %s: %r",
            signature_name,
            e,
        )
        return None


def score_signature_against_dataset(
    signature: Signature,
    dataset_species: Set[str],
    dataset_directions: Optional[Dict[str, str]] = None,
    metadata: Optional[Dict[str, Any]] = None,
) -> SignatureScoreResult:
    """
    Score a signature against a dataset.

    Computes overlap/similarity score between a lipid signature and a dataset.

    Args:
        signature: Signature definition
        dataset_species: Set of species names found in the dataset
        dataset_directions: Optional dict mapping species → direction (↑/↓)
        metadata: Optional dataset metadata (for future use)

    Returns:
        SignatureScoreResult with detailed scoring information
    """
    return score_signature(
        signature=signature,
        dataset_species=dataset_species,
        dataset_directions=dataset_directions,
    )


def find_matching_signatures_for_dataset(
    dataset_species: Set[str],
    dataset_directions: Optional[Dict[str, str]] = None,
    overlap_threshold: float = 0.3,
) -> List[SignatureMatchResult]:
    """
    Find all signatures that match a dataset above the overlap threshold.

    Args:
        dataset_species: Set of species names in the dataset
        dataset_directions: Optional dict mapping species → direction
        overlap_threshold: Minimum overlap fraction to consider a match (default: 0.3)

    Returns:
        List of SignatureMatchResult objects for matching signatures
    """
    matches: List[SignatureMatchResult] = []

    # Fetch all signatures from Notion
    signature_pages = fetch_all_signatures_from_notion()

    if not signature_pages:
        logger.debug(
            "[INGEST][SIGNATURE-MATCH] No signatures found in Notion",
        )
        return matches

    for sig_page in signature_pages:
        try:
            # Load signature from Notion page
            signature = load_signature_from_notion_page(sig_page)
            if not signature:
                continue

            # Score signature against dataset
            score_result = score_signature_against_dataset(
                signature=signature,
                dataset_species=dataset_species,
                dataset_directions=dataset_directions,
            )

            # Calculate overlap fraction
            total_components = len(signature.components)
            matched_count = len(score_result.matched_species)
            overlap_fraction = (
                matched_count / total_components if total_components > 0 else 0.0
            )

            # Check if above threshold
            if overlap_fraction >= overlap_threshold:
                props = sig_page.get("properties", {}) or {}
                name_prop = props.get("Name", {}).get("title", []) or []
                signature_name = (
                    name_prop[0].get("plain_text", "") if name_prop else "Unknown"
                )

                match_result = SignatureMatchResult(
                    signature_page_id=sig_page.get("id", ""),
                    signature_name=signature_name,
                    score=score_result.total_score,
                    overlap_fraction=overlap_fraction,
                    matched_components=score_result.matched_species,
                    missing_components=score_result.missing_species,
                    conflicting_components=score_result.conflicting_species,
                    score_result=score_result,
                )
                matches.append(match_result)

                logger.info(
                    "[INGEST][SIGNATURE-MATCH] Found match: %s (overlap: %.2f, score: %.3f)",
                    signature_name,
                    overlap_fraction,
                    score_result.total_score,
                )

        except Exception as e:
            logger.warning(
                "[INGEST][SIGNATURE-MATCH] Error processing signature %s: %r",
                sig_page.get("id", ""),
                e,
            )
            continue

    return matches


def update_dataset_with_signature_matches(
    dataset_page_id: str,
    matches: List[SignatureMatchResult],
) -> None:
    """
    Update a dataset's Notion page with signature match information.

    Adds:
    - Relations to matching signatures
    - Summary text of matches
    - Highest match score

    Args:
        dataset_page_id: Notion page ID of dataset (with dashes)
        matches: List of signature match results
    """
    if not matches:
        logger.debug(
            "[INGEST][SIGNATURE-MATCH] No matches to write back for dataset %s",
            dataset_page_id,
        )
        return

    cfg = get_config()

    try:
        # Build signature relation list
        signature_relations = [{"id": match.signature_page_id} for match in matches]

        # Build summary text
        summary_parts = [
            f"Found {len(matches)} matching signature(s):",
        ]

        for match in sorted(matches, key=lambda m: m.score, reverse=True):
            summary_parts.append(
                f"\n• {match.signature_name}: "
                f"{len(match.matched_components)}/{len(match.matched_components) + len(match.missing_components)} "
                f"components matched (score: {match.score:.3f})"
            )
            if match.missing_components:
                summary_parts.append(
                    f"  Missing: {', '.join(match.missing_components[:5])}"
                    + ("..." if len(match.missing_components) > 5 else "")
                )

        summary_text = "\n".join(summary_parts)

        # Get highest score
        highest_score = max(match.score for match in matches) if matches else 0.0

        # Update Notion page
        url = f"{cfg.notion.base_url}/pages/{dataset_page_id}"

        properties: Dict[str, Any] = {}

        # Add relation to matching signatures
        # Use correct property name from schema: "Related Signature(s)"
        if signature_relations:
            properties["Related Signature(s)"] = {
                "relation": signature_relations,
            }

        # Add summary to existing Summary property (append to existing content)
        # First fetch current summary to append
        try:
            page_url = f"{cfg.notion.base_url}/pages/{dataset_page_id}"
            page_resp = requests.get(page_url, headers=notion_headers(), timeout=30)
            page_resp.raise_for_status()
            page_data = page_resp.json()
            page_props = page_data.get("properties", {}) or {}
            summary_prop = page_props.get("Summary", {}) or {}
            existing_summary_rich = summary_prop.get("rich_text", []) or []
            existing_summary = "".join(
                rt.get("plain_text", "") for rt in existing_summary_rich
            )

            # Append signature match summary
            if existing_summary and not existing_summary.strip().endswith("."):
                existing_summary += "."
            combined_summary = (
                f"{existing_summary}\n\n## Signature Matches\n\n{summary_text}"
                if existing_summary
                else summary_text
            )

            properties["Summary"] = {
                "rich_text": [
                    {
                        "text": {
                            "content": combined_summary,
                        },
                    },
                ],
            }
        except Exception as e:
            logger.warning(
                "[INGEST][SIGNATURE-MATCH] Could not fetch existing summary for %s: %r",
                dataset_page_id,
                e,
            )
            # Fallback: just add summary
            properties["Summary"] = {
                "rich_text": [
                    {
                        "text": {
                            "content": summary_text,
                        },
                    },
                ],
            }

        payload = {
            "properties": properties,
        }

        resp = requests.patch(
            url,
            headers=notion_headers(),
            json=payload,
            timeout=30,
        )

        if resp.status_code >= 300:
            logger.warning(
                "[INGEST][SIGNATURE-MATCH] Failed to update dataset %s: %s",
                dataset_page_id,
                resp.text,
            )
        else:
            logger.info(
                "[INGEST][SIGNATURE-MATCH] Updated dataset %s with %d signature match(es)",
                dataset_page_id,
                len(matches),
            )

    except Exception as e:
        logger.warning(
            "[INGEST][SIGNATURE-MATCH] Error updating dataset %s with matches: %r",
            dataset_page_id,
            e,
        )
