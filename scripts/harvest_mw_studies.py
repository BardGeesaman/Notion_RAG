#!/usr/bin/env python3
"""
Harvest Metabolomics Workbench studies and create/update Dataset pages in Postgres.

This script fetches study metadata and mwTab data from the Metabolomics Workbench
REST API and creates corresponding datasets in Postgres.

Notion support has been removed - Postgres is now the source of truth.
"""

import argparse
import json
from typing import Any, Dict, List, Optional

import requests

from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def notion_headers() -> Dict[str, str]:
    """Stub: Notion support removed. Returns empty headers."""
    return {}

MW_BASE_URL = "https://www.metabolomicsworkbench.org/rest"
MW_STUDY_SUMMARY_URL = f"{MW_BASE_URL}/study/study_id/ST/summary"
MW_REFMET_URL_TEMPLATE = f"{MW_BASE_URL}/study/refmet_name/{{refmet_name}}/data/txt"


def _fetch_all_study_summaries() -> List[Dict[str, Any]]:
    """
    Fetch all public study summaries from MW.

    Uses the documented endpoint /rest/study/study_id/ST/summary
    which returns a list of studies in JSON.

    Returns:
        List of study summary dictionaries
    """
    try:
        resp = requests.get(MW_STUDY_SUMMARY_URL, params={"format": "json"}, timeout=60)
        resp.raise_for_status()
        data = resp.json()

        # Expect data to be a list of dicts; if it's a dict keyed by study_id, normalize.
        if isinstance(data, dict):
            # Some MW endpoints return { "ST001111": {...}, "ST001112": {...} }
            studies = []
            for sid, summary in data.items():
                if isinstance(summary, dict):
                    summary.setdefault("study_id", sid)
                    studies.append(summary)
            return studies

        if isinstance(data, list):
            return data

        logger.warning("[MW] Unexpected response format from study summaries endpoint")
        return []

    except Exception as e:
        logger.error("[MW] Error fetching all study summaries: %r", e)
        raise


def fetch_mw_study_summary(study_id: str) -> Dict[str, Any]:
    """
    Fetch study summary metadata from Metabolomics Workbench REST API.

    Args:
        study_id: MW study ID (e.g., "ST001111")

    Returns:
        Dictionary with study metadata
    """
    url = f"{MW_BASE_URL}/study/study_id/{study_id}/summary"

    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()

        # Handle both JSON and text responses
        content_type = resp.headers.get("content-type", "")
        if "application/json" in content_type:
            data = resp.json()
            # MW API may return a list or wrapped dict
            if isinstance(data, list) and data:
                return data[0]
            elif isinstance(data, dict) and "summary" in data:
                return data["summary"]
            return data if isinstance(data, dict) else {}
        else:
            # Try to parse as JSON anyway
            try:
                return resp.json()
            except:
                logger.warning("[MW] Study summary for %s not in JSON format", study_id)
                return {}
    except Exception as e:
        logger.error("[MW] Error fetching study summary for %s: %r", study_id, e)
        raise


def fetch_mw_mwtab(study_id: str) -> str:
    """
    Fetch mwTab data from Metabolomics Workbench REST API.

    Args:
        study_id: MW study ID (e.g., "ST001111")

    Returns:
        mwTab text content
    """
    url = f"{MW_BASE_URL}/study/study_id/{study_id}/mwtab"

    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()

        # Always try to get as text first
        text_content = resp.text

        # If content type suggests JSON, try parsing but fallback to raw text
        content_type = resp.headers.get("content-type", "")
        if "application/json" in content_type:
            try:
                data = resp.json()
                if isinstance(data, dict):
                    return json.dumps(data, indent=2)
                elif isinstance(data, list):
                    return json.dumps(data, indent=2)
                return str(data)
            except (json.JSONDecodeError, ValueError):
                # If JSON parsing fails, return raw text
                logger.debug(
                    "[MW] mwTab response looks like JSON but parse failed, using raw text"
                )
                return text_content

        return text_content
    except Exception as e:
        logger.warning("[MW] Error fetching mwTab for %s: %r", study_id, e)
        # Try alternative endpoint
        try:
            url = f"{MW_BASE_URL}/study/study_id/{study_id}/data"
            resp = requests.get(url, timeout=30)
            resp.raise_for_status()
            return resp.text
        except Exception as e2:
            logger.warning(
                "[MW] Error fetching data from alternative endpoint for %s: %r",
                study_id,
                e2,
            )
            # Return empty string instead of raising - allows page creation to continue
            return ""


def extract_mw_metadata(summary: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extract relevant metadata from MW study summary.

    Args:
        summary: Raw study summary from MW API

    Returns:
        Dictionary with extracted metadata fields
    """
    # MW API structure may vary, so we'll try common fields
    result = {}

    # Try different possible field names
    result["title"] = (
        summary.get("study_title")
        or summary.get("title")
        or summary.get("study_name")
        or summary.get("name")
        or "(Untitled Study)"
    )

    result["organism"] = (
        summary.get("organism")
        or summary.get("species")
        or summary.get("sample_species")
        or ""
    )

    result["sample_type"] = (
        summary.get("sample_type")
        or summary.get("sample_source")
        or summary.get("matrix")
        or ""
    )

    result["disease"] = (
        summary.get("disease")
        or summary.get("condition")
        or summary.get("disease_condition")
        or ""
    )

    result["summary"] = (
        summary.get("study_summary")
        or summary.get("summary")
        or summary.get("description")
        or ""
    )

    result["doi"] = summary.get("doi") or summary.get("DOI") or ""
    result["pubmed_id"] = (
        summary.get("pubmed_id") or summary.get("pubmed") or summary.get("pmid") or ""
    )
    result["url"] = summary.get("url") or summary.get("link") or ""

    return result


def _get_database_schema() -> Optional[Dict[str, Any]]:
    """
    Fetch the Experimental Data Assets database schema to see available properties.

    Returns:
        Database schema dict or None if error
    """
    cfg = get_config().notion
    if not cfg.exp_data_db_id:
        return None

    try:
        url = f"{cfg.base_url}/databases/{cfg.exp_data_db_id}"
        resp = requests.get(url, headers=notion_headers(), timeout=30)
        resp.raise_for_status()
        return resp.json()
    except Exception as e:
        logger.error("[MW] Error fetching database schema: %r", e)
        return None


def find_existing_dataset_page(study_id: str) -> Optional[str]:
    """
    Find existing Dataset page by MW Study ID.

    Args:
        study_id: MW study ID to search for

    Returns:
        Page ID if found, None otherwise
    """
    cfg = get_config().notion
    if not cfg.exp_data_db_id:
        logger.error("[MW] Experimental Data Assets DB ID not configured")
        return None

    # Try multiple possible property names
    possible_property_names = [
        "MW Study ID",
        "MW Study",
        "Study ID",
        "Metabolomics Workbench Study ID",
    ]

    query_url = f"{cfg.base_url}/databases/{cfg.exp_data_db_id}/query"

    for prop_name in possible_property_names:
        payload = {
            "page_size": 1,
            "filter": {
                "property": prop_name,
                "rich_text": {"equals": study_id},
            },
        }

        try:
            resp = requests.post(
                query_url, headers=notion_headers(), json=payload, timeout=30
            )
            if resp.status_code == 200:
                results = resp.json().get("results", [])
                if results:
                    return results[0]["id"]
        except Exception:
            # Try next property name
            continue

    # Try searching in Summary field for MW Study ID
    payload = {
        "page_size": 100,  # Get more pages to search through
    }

    try:
        resp = requests.post(
            query_url, headers=notion_headers(), json=payload, timeout=30
        )
        resp.raise_for_status()
        results = resp.json().get("results", [])

        # Search through Summary field for the study ID
        for page in results:
            summary_prop = page.get("properties", {}).get("Summary", {})
            rich_text = summary_prop.get("rich_text", [])
            if rich_text:
                summary_content = "".join(rt.get("plain_text", "") for rt in rich_text)
                if f"MW Study ID: {study_id}" in summary_content:
                    return page["id"]
    except Exception as e:
        if hasattr(e, "response") and e.response is not None:
            logger.warning(
                "[MW] Error querying for existing page: %s - Response: %s",
                str(e),
                e.response.text,
            )
        else:
            logger.warning("[MW] Error querying for existing page: %r", e)

    return None


def create_dataset_page(study_id: str, mw_meta: Dict[str, Any], mwtab_text: str) -> str:
    """
    Create a new Dataset page in Notion.

    Args:
        study_id: MW study ID
        mw_meta: Extracted MW metadata
        mwtab_text: mwTab text content

    Returns:
        Created page ID
    """
    cfg = get_config().notion
    if not cfg.exp_data_db_id:
        raise RuntimeError("Experimental Data Assets DB ID not configured")

    # Build properties - match actual database schema
    props: Dict[str, Any] = {
        "Experiment Name": {"title": [{"text": {"content": mw_meta["title"]}}]},
        "Data Origin": {"select": {"name": "External – Published"}},
        "Dataset Source Type": {"select": {"name": "Processed table"}},
    }

    # Add MW Study ID to Summary if available (since MW Study ID property doesn't exist)
    if mw_meta.get("title"):
        summary_text = f"MW Study ID: {study_id}\n"
        if mw_meta.get("summary"):
            summary_text += f"\n{mw_meta['summary']}"
        props["Summary"] = {"rich_text": [{"text": {"content": summary_text}}]}

    # Source URL
    if mw_meta.get("doi"):
        props["Source URL / DOI"] = {"url": mw_meta["doi"]}
    elif mw_meta.get("url"):
        props["Source URL / DOI"] = {"url": mw_meta["url"]}
    else:
        props["Source URL / DOI"] = {
            "url": f"https://www.metabolomicsworkbench.org/study/index.php?study_id={study_id}"
        }

    # Disease (multi-select)
    if mw_meta.get("disease"):
        props["Disease"] = {"multi_select": [{"name": mw_meta["disease"]}]}

    # Matrix (multi-select)
    if mw_meta.get("sample_type"):
        props["Matrix"] = {"multi_select": [{"name": mw_meta["sample_type"]}]}

    # Model Systems (multi-select)
    if mw_meta.get("organism"):
        props["Model Systems"] = {"multi_select": [{"name": mw_meta["organism"]}]}

    # Create page
    payload = {
        "parent": {"database_id": cfg.exp_data_db_id},
        "properties": props,
    }

    try:
        resp = requests.post(
            f"{cfg.base_url}/pages",
            headers=notion_headers(),
            json=payload,
            timeout=30,
        )
        resp.raise_for_status()
        page_id = resp.json()["id"]
        logger.info("[MW] Created Dataset page %s for study %s", page_id, study_id)
        return page_id
    except Exception as e:
        if hasattr(e, "response") and e.response is not None:
            logger.error(
                "[MW] Error creating Dataset page for %s: %s - Response: %s",
                study_id,
                str(e),
                e.response.text,
            )
        else:
            logger.error("[MW] Error creating Dataset page for %s: %r", study_id, e)
        raise


def update_dataset_page(page_id: str, mw_meta: Dict[str, Any]) -> None:
    """
    Update existing Dataset page with new metadata.

    Args:
        page_id: Notion page ID to update
        mw_meta: Extracted MW metadata
    """
    props: Dict[str, Any] = {}

    # Update title - use "Experiment Name" to match schema
    if mw_meta.get("title"):
        props["Experiment Name"] = {"title": [{"text": {"content": mw_meta["title"]}}]}

    # Update other fields if needed (similar to create)
    if mw_meta.get("doi"):
        props["Source URL / DOI"] = {"url": mw_meta["doi"]}
    elif mw_meta.get("url"):
        props["Source URL / DOI"] = {"url": mw_meta["url"]}

    # Disease (multi-select)
    if mw_meta.get("disease"):
        props["Disease"] = {"multi_select": [{"name": mw_meta["disease"]}]}

    # Matrix (multi-select)
    if mw_meta.get("sample_type"):
        props["Matrix"] = {"multi_select": [{"name": mw_meta["sample_type"]}]}

    # Model Systems (multi-select)
    if mw_meta.get("organism"):
        props["Model Systems"] = {"multi_select": [{"name": mw_meta["organism"]}]}

    if props:
        cfg = get_config().notion
        payload = {"properties": props}

        try:
            resp = requests.patch(
                f"{cfg.base_url}/pages/{page_id}",
                headers=notion_headers(),
                json=payload,
                timeout=30,
            )
            resp.raise_for_status()
            logger.info("[MW] Updated Dataset page %s", page_id)
        except Exception as e:
            if hasattr(e, "response") and e.response is not None:
                logger.error(
                    "[MW] Error updating Dataset page %s: %s - Response: %s",
                    page_id,
                    str(e),
                    e.response.text,
                )
            else:
                logger.error("[MW] Error updating Dataset page %s: %r", page_id, e)
            raise


def add_mwtab_block(page_id: str, mwtab_text: str) -> None:
    """
    Add mwTab content as blocks on the Dataset page.

    Splits content into <=2000 char chunks per code block to ensure all content
    is preserved and readable by extract_page_content().

    Args:
        page_id: Notion page ID
        mwtab_text: mwTab text content
    """
    if not mwtab_text or not mwtab_text.strip():
        logger.warning("[MW] Empty mwTab content, skipping block creation")
        return

    cfg = get_config().notion
    max_block_length = 2000  # Notion code block limit

    # Create blocks: heading + code blocks
    blocks = [
        {
            "object": "block",
            "type": "heading_2",
            "heading_2": {
                "rich_text": [{"type": "text", "text": {"content": "mwTab Data"}}]
            },
        },
    ]

    # Split mwTab text into chunks of <=2000 characters
    # Try to split at newlines to avoid breaking lines in the middle
    chunks = []
    remaining = mwtab_text

    while len(remaining) > max_block_length:
        # Find the last newline within the chunk limit
        chunk = remaining[:max_block_length]
        last_newline = chunk.rfind("\n")

        if (
            last_newline > max_block_length * 0.8
        ):  # Only split at newline if we're not losing too much
            chunks.append(remaining[: last_newline + 1])
            remaining = remaining[last_newline + 1 :]
        else:
            # Force split if no good newline found
            chunks.append(chunk)
            remaining = remaining[max_block_length:]

    # Add the last chunk
    if remaining:
        chunks.append(remaining)

    logger.info(
        "[MW] Splitting mwTab into %d code block(s) for page %s", len(chunks), page_id
    )

    # Create a code block for each chunk
    for i, chunk in enumerate(chunks):
        blocks.append(
            {
                "object": "block",
                "type": "code",
                "code": {
                    "rich_text": [{"type": "text", "text": {"content": chunk}}],
                    "language": "plain text",
                },
            }
        )

    # Notion API allows adding up to 100 children at once
    # Split into batches if needed
    batch_size = 100
    for i in range(0, len(blocks), batch_size):
        batch = blocks[i : i + batch_size]

        try:
            resp = requests.patch(
                f"{cfg.base_url}/blocks/{page_id}/children",
                headers=notion_headers(),
                json={"children": batch},
                timeout=30,
            )
            resp.raise_for_status()

            if i == 0:
                logger.info(
                    "[MW] Added mwTab heading and first batch of blocks to page %s",
                    page_id,
                )
            else:
                logger.info(
                    "[MW] Added batch %d of mwTab blocks to page %s",
                    (i // batch_size) + 1,
                    page_id,
                )
        except Exception as e:
            if hasattr(e, "response") and e.response is not None:
                logger.error(
                    "[MW] Error adding mwTab blocks to page %s (batch %d): %s - Response: %s",
                    page_id,
                    (i // batch_size) + 1,
                    str(e),
                    e.response.text,
                )
            else:
                logger.error(
                    "[MW] Error adding mwTab blocks to page %s (batch %d): %r",
                    page_id,
                    (i // batch_size) + 1,
                    e,
                )
            # Continue with remaining batches even if one fails


def _search_mw_study_ids_by_keyword(
    keyword: str, summaries: Optional[List[Dict[str, Any]]] = None
) -> List[str]:
    """
    Filter all MW studies by keyword in title / summary / disease fields.

    This is a simple client-side filter over /study/study_id/ST/summary output.

    Args:
        keyword: Search keyword (e.g., "ALS", "amyotrophic")
        summaries: Optional pre-fetched study summaries. If None, fetches all summaries.

    Returns:
        List of MW study IDs matching the keyword
    """
    keyword_lower = keyword.lower()
    study_ids: List[str] = []

    if summaries is None:
        summaries = _fetch_all_study_summaries()

    for row in summaries:
        sid = row.get("study_id") or row.get("STUDY_ID") or row.get("Study_ID")
        if not sid:
            continue

        title = row.get("study_title") or row.get("Study_title") or ""
        summary_text = row.get("study_summary") or row.get("Study_summary") or ""
        disease = row.get("disease") or row.get("Disease") or ""

        combined = f"{title} {summary_text} {disease}".lower()
        if keyword_lower in combined:
            study_ids.append(sid)

    logger.info("[MW] Found %d studies matching keyword '%s'", len(study_ids), keyword)
    return study_ids


def _study_contains_lipids_matching_filters(
    study_id: str, lipid_filters: List[str]
) -> bool:
    """
    Use MW's /study/refmet_name/{name}/data/txt endpoint to check whether a given study_id
    is present among studies containing specific refmet_name(s).

    lipid_filters are strings like 'Ceramide', 'Sphingomyelin', etc.

    Args:
        study_id: MW study ID to check
        lipid_filters: List of RefMet name filter strings (e.g., ["Ceramide", "Sphingomyelin"])

    Returns:
        True if study contains metabolites matching any filter, False otherwise
    """
    study_id_upper = study_id.upper()

    for lf in lipid_filters:
        refmet_name = lf.strip()
        if not refmet_name:
            continue

        url = MW_REFMET_URL_TEMPLATE.format(refmet_name=refmet_name)
        try:
            resp = requests.get(url, timeout=60)
            if resp.status_code != 200:
                logger.debug(
                    "[MW] RefMet endpoint returned %d for %s",
                    resp.status_code,
                    refmet_name,
                )
                continue  # skip silently for now

            text = resp.text
            # The /data/txt endpoint returns tab-separated values with header row
            # Format: refmet_name\tkegg_id\tstudy_id
            # Parse and extract study IDs from the third column
            lines = text.strip().split("\n")
            if len(lines) < 2:
                continue  # Only header or empty

            # Skip header and check each line for the study ID
            for line in lines[1:]:  # Skip header row
                parts = line.split("\t")
                if len(parts) >= 3:
                    # Third column is study_id
                    found_study_id = parts[2].strip().upper()
                    if found_study_id == study_id_upper:
                        logger.debug(
                            "[MW] Study %s found in RefMet data for %s",
                            study_id,
                            refmet_name,
                        )
                        return True
        except Exception as e:
            logger.debug(
                "[MW] Error checking RefMet endpoint for %s: %r", refmet_name, e
            )
            continue

    return False


def discover_study_ids_by_keyword(
    keywords: List[str],
    max_results: int = 50,
    lipid_filters: Optional[List[str]] = None,
) -> List[str]:
    """
    Use MW REST API to search for studies by keyword(s), optionally filter by
    metabolite names containing specified lipid filter strings (e.g. 'ceramide').

    Args:
        keywords: List of keywords to search for (e.g., ["ALS", "amyotrophic"])
        max_results: Maximum number of studies to return
        lipid_filters: Optional list of lipid filter strings (e.g., ["ceramide"])

    Returns:
        List of MW study IDs matching the criteria
    """
    if lipid_filters is None:
        lipid_filters = []

    # 1) Fetch all study summaries once (optimization)
    logger.info("[MW] Fetching all study summaries for keyword search...")
    summaries = _fetch_all_study_summaries()
    logger.info("[MW] Loaded %d study summaries", len(summaries))

    # 2) Query MW for candidate studies using keywords
    candidate_ids: List[str] = []
    seen_ids = set()

    for kw in keywords:
        ids_for_kw = _search_mw_study_ids_by_keyword(kw, summaries=summaries)
        for sid in ids_for_kw:
            if sid not in seen_ids:
                candidate_ids.append(sid)
                seen_ids.add(sid)

        if len(candidate_ids) >= max_results:
            break

    # Truncate to max_results
    study_ids = candidate_ids[:max_results]

    logger.info("[MW] Found %d unique studies from keyword search", len(study_ids))

    if not lipid_filters:
        return study_ids

    # 2) Filter studies by metabolite names
    logger.info("[MW] Filtering studies by lipid filters: %s", lipid_filters)
    filtered_ids: List[str] = []

    for sid in study_ids:
        try:
            if _study_contains_lipids_matching_filters(sid, lipid_filters):
                filtered_ids.append(sid)
        except Exception as e:
            logger.warning("[MW] Failed to check metabolites for %s: %r", sid, e)
            continue

    logger.info("[MW] %d studies passed lipid filter", len(filtered_ids))
    return filtered_ids


def harvest_study(
    study_id: str,
    create_notion: bool = False,
    ingest: bool = False,
    dry_run: bool = False,
) -> Optional[str]:
    """
    Harvest a single MW study and optionally create/update Notion page.

    Args:
        study_id: MW study ID (e.g., "ST001111")
        create_notion: If True, create/update Notion page
        ingest: If True, trigger dataset ingestion after creating page
        dry_run: If True, only print what would be done

    Returns:
        Created/updated page ID if successful, None otherwise
    """
    logger.info("[MW] Harvesting study %s", study_id)

    try:
        # Fetch MW data
        summary = fetch_mw_study_summary(study_id)
        mw_meta = extract_mw_metadata(summary)

        if dry_run:
            print(f"\n📊 Study: {study_id}")
            print(f"  Title: {mw_meta['title']}")
            print(f"  Organism: {mw_meta.get('organism', 'N/A')}")
            print(f"  Sample Type: {mw_meta.get('sample_type', 'N/A')}")
            print(f"  Disease: {mw_meta.get('disease', 'N/A')}")
            print(f"  Would create/update Notion page")
            return None

        if not create_notion:
            logger.info("[MW] Skipping Notion creation (--create-notion not set)")
            return None

        # Check for existing page
        existing_page_id = find_existing_dataset_page(study_id)

        # Fetch mwTab (needed for both new and existing pages)
        try:
            mwtab_text = fetch_mw_mwtab(study_id)
        except Exception as e:
            logger.warning(
                "[MW] Could not fetch mwTab for %s, continuing without it: %r",
                study_id,
                e,
            )
            mwtab_text = ""

        if existing_page_id:
            logger.info(
                "[MW] Found existing page %s for study %s", existing_page_id, study_id
            )
            page_id = existing_page_id
            update_dataset_page(page_id, mw_meta)

            # Add mwTab blocks to existing page if available
            if mwtab_text:
                add_mwtab_block(page_id, mwtab_text)
        else:
            # Create new page
            page_id = create_dataset_page(study_id, mw_meta, mwtab_text)

            # Add mwTab blocks if available
            if mwtab_text:
                add_mwtab_block(page_id, mwtab_text)

        # Optional ingestion
        if ingest:
            logger.info("[MW] Triggering ingestion for page %s", page_id)
            try:
                from amprenta_rag.ingestion.dataset_ingestion import \
                    ingest_dataset

                ingest_dataset(page_id, force=True)
            except Exception as e:
                logger.error("[MW] Error ingesting dataset %s: %r", page_id, e)
                raise

        return page_id

    except Exception as e:
        logger.error("[MW] Error harvesting study %s: %r", study_id, e)
        raise


def main():
    parser = argparse.ArgumentParser(
        description="Harvest Metabolomics Workbench studies and create/update Dataset pages in Notion."
    )
    parser.add_argument(
        "--study-id",
        dest="study_ids",
        action="append",
        default=[],
        help="MW study ID (e.g., ST001111). Can be specified multiple times.",
    )
    parser.add_argument(
        "--search-keyword",
        action="append",
        default=[],
        help=(
            "Keyword to search in MW study title/summary/disease fields "
            "(e.g. ALS, amyotrophic, 'motor neuron'). "
            "Can be supplied multiple times."
        ),
    )
    parser.add_argument(
        "--max-search-results",
        type=int,
        default=50,
        help="Maximum number of studies to consider from keyword search.",
    )
    parser.add_argument(
        "--lipid-filter",
        action="append",
        default=[],
        help=(
            "Optional RefMet name filter(s) (e.g. Ceramide, Sphingomyelin). "
            "Uses /study/refmet_name/{name}/data/txt to restrict to studies "
            "containing those lipids."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be done without creating/updating Notion pages.",
    )
    parser.add_argument(
        "--create-notion",
        action="store_true",
        help="Create or update pages in Notion.",
    )
    parser.add_argument(
        "--ingest",
        action="store_true",
        help="Trigger dataset ingestion after creating/updating pages.",
    )

    args = parser.parse_args()

    if args.ingest and not args.create_notion:
        parser.error("--ingest requires --create-notion")

    # Determine study IDs: explicit mode or search mode
    if args.study_ids:
        # Explicit mode: use provided study IDs
        study_ids = args.study_ids
        logger.info("[MW] Using explicit study IDs: %s", study_ids)
    elif args.search_keyword:
        # Search mode: discover study IDs by keyword
        logger.info("[MW] Searching for studies with keywords: %s", args.search_keyword)
        study_ids = discover_study_ids_by_keyword(
            keywords=args.search_keyword,
            max_results=args.max_search_results,
            lipid_filters=args.lipid_filter if args.lipid_filter else None,
        )

        if not study_ids:
            print("No studies found matching the provided keywords and lipid filters.")
            return

        print(f"\n🔍 Found {len(study_ids)} candidate study ID(s):")
        for sid in study_ids:
            print(f"  - {sid}")
        print()

        if args.dry_run:
            print("(Dry-run mode: not creating/updating Notion pages)")
            return
    else:
        parser.error("You must provide either --study-id or --search-keyword")

    for study_id in study_ids:
        try:
            harvest_study(
                study_id=study_id,
                create_notion=args.create_notion,
                ingest=args.ingest,
                dry_run=args.dry_run,
            )
        except Exception as e:
            logger.error("[MW] Failed to harvest study %s: %r", study_id, e)
            if not args.dry_run:
                raise


if __name__ == "__main__":
    main()
