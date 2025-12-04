"""
Notion writeback for signature matches.

Handles updating dataset pages in Notion with signature match information.
"""

from __future__ import annotations

from typing import Any, Dict, List

import requests

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config
from amprenta_rag.ingestion.signature_matching.models import SignatureMatchResult
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


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
        logger.info(
            "[INGEST][SIGNATURE-MATCH] No signature matches found for dataset %s — Signature Match Score unchanged",
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

        # Write Signature Match Score (numeric property)
        # Only write if we have matches (idempotent: don't clear if no matches)
        if matches and highest_score > 0:
            properties["Signature Match Score"] = {
                "number": highest_score,
            }
            logger.info(
                "[INGEST][SIGNATURE-MATCH] Writing Signature Match Score = %.3f to dataset %s",
                highest_score,
                dataset_page_id,
            )

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
            # Check if error is due to missing property
            error_text = resp.text
            if "Signature Match Score" in error_text and "not a property" in error_text:
                logger.warning(
                    "[INGEST][SIGNATURE-MATCH] Signature Match Score property not found on dataset %s. "
                    "Property may not exist in database schema. Continuing without score update.",
                    dataset_page_id,
                )
                # Retry without the Signature Match Score property
                if "Signature Match Score" in properties:
                    del properties["Signature Match Score"]
                    payload_retry = {"properties": properties}
                    retry_resp = requests.patch(
                        url,
                        headers=notion_headers(),
                        json=payload_retry,
                        timeout=30,
                    )
                    if retry_resp.status_code < 300:
                        logger.info(
                            "[INGEST][SIGNATURE-MATCH] Updated dataset %s with %d signature match(es) "
                            "(without Signature Match Score property)",
                            dataset_page_id,
                            len(matches),
                        )
                    else:
                        logger.warning(
                            "[INGEST][SIGNATURE-MATCH] Failed to update dataset %s even without score: %s",
                            dataset_page_id,
                            retry_resp.text,
                        )
            else:
                logger.warning(
                    "[INGEST][SIGNATURE-MATCH] Failed to update dataset %s: %s",
                    dataset_page_id,
                    error_text,
                )
        else:
            logger.info(
                "[INGEST][SIGNATURE-MATCH] Updated dataset %s with %d signature match(es) "
                "(Signature Match Score = %.3f)",
                dataset_page_id,
                len(matches),
                highest_score if matches else 0.0,
            )

    except Exception as e:
        logger.warning(
            "[INGEST][SIGNATURE-MATCH] Error updating dataset %s with matches: %r",
            dataset_page_id,
            e,
        )

