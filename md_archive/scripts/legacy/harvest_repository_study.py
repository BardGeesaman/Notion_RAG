#!/usr/bin/env python3
"""
Unified harvest script for repository studies.

Harvests studies from any repository (GEO, PRIDE, MetaboLights, MW)
and creates Postgres datasets.

Notion support has been removed - Postgres is now the source of truth.
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, Optional

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.config import get_config
from amprenta_rag.ingestion.postgres_dataset_ingestion import ingest_dataset_from_postgres
from amprenta_rag.ingestion.postgres_integration import create_or_update_dataset_in_postgres
from amprenta_rag.ingestion.repositories.discovery import fetch_study_metadata, get_repository
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.models.domain import OmicsType

logger = get_logger(__name__)

# Stub for removed Notion support
NOTION_EXP_DATA_DB_ID = None


def notion_headers() -> Dict[str, str]:
    """Stub: Notion support removed. Returns empty headers."""
    return {}


def find_existing_dataset_page(study_id: str, repository: str) -> str | None:
    """
    Find existing Dataset page by study ID.

    Args:
        study_id: Repository-specific study identifier
        repository: Repository name

    Returns:
        Notion page ID if found, None otherwise
    """
    exp_data_db_id = NOTION_EXP_DATA_DB_ID
    if not exp_data_db_id:
        logger.error("[HARVEST] Experimental Data Assets DB ID not configured")
        return None

    import requests

    # Try different property names for study ID
    property_candidates = [
        "MW Study ID",
        "GEO Study ID",
        "PRIDE Project ID",
        "MetaboLights Study ID",
        "Study ID",
    ]

    for prop_name in property_candidates:
        try:
            cfg = get_config()
            url = f"{cfg.notion.base_url}/databases/{exp_data_db_id}/query"
            payload = {
                "filter": {
                    "property": prop_name,
                    "rich_text": {"equals": study_id},
                },
                "page_size": 1,
            }

            resp = requests.post(
                url,
                headers=notion_headers(),
                json=payload,
                timeout=30,
            )
            resp.raise_for_status()

            results = resp.json().get("results", [])
            if results:
                page_id = results[0].get("id", "")
                logger.info(
                    "[HARVEST] Found existing page %s for %s study %s",
                    page_id,
                    repository,
                    study_id,
                )
                return page_id
        except Exception as e:
            logger.debug(
                "[HARVEST] Error querying for property %s: %r",
                prop_name,
                e,
            )
            continue

    # Also try searching in Summary field
    try:
        cfg = get_config()
        url = f"{cfg.notion.base_url}/databases/{exp_data_db_id}/query"
        payload = {
            "filter": {
                "property": "Summary",
                "rich_text": {"contains": study_id},
            },
            "page_size": 10,
        }

        resp = requests.post(
            url,
            headers=notion_headers(),
            json=payload,
            timeout=30,
        )
        resp.raise_for_status()

        results = resp.json().get("results", [])
        for result in results:
            props = result.get("properties", {})
            summary_prop = props.get("Summary", {})
            if summary_prop.get("type") == "rich_text":
                rich_text = summary_prop.get("rich_text", [])
                summary_content = "".join(
                    rt.get("plain_text", "") for rt in rich_text
                )
                if f"{repository} Study ID: {study_id}" in summary_content or f"Study ID: {study_id}" in summary_content:
                    page_id = result.get("id", "")
                    logger.info(
                        "[HARVEST] Found existing page %s for %s study %s (in Summary)",
                        page_id,
                        repository,
                        study_id,
                    )
                    return page_id
    except Exception as e:
        logger.debug("[HARVEST] Error searching Summary field: %r", e)

    return None


def create_dataset_page_from_metadata(
    metadata,
    repository: str,
) -> str | None:
    """
    Create a Notion Dataset page from StudyMetadata.

    Args:
        metadata: StudyMetadata object
        repository: Repository name

    Returns:
        Notion page ID if created, None otherwise
    """
    exp_data_db_id = NOTION_EXP_DATA_DB_ID
    if not exp_data_db_id:
        logger.error("[HARVEST] Experimental Data Assets DB ID not configured")
        return None

    import requests

    # Build page properties
    props = {
        "Experiment Name": {"title": [{"text": {"content": metadata.title}}]},
    }

    # Add study ID to Summary
    summary_text = f"{repository} Study ID: {metadata.study_id}\n"
    if metadata.summary:
        summary_text += f"\n{metadata.summary}"
    props["Summary"] = {"rich_text": [{"text": {"content": summary_text}}]}

    # Add DOI if available
    if metadata.doi:
        props["Source URL / DOI"] = {"url": metadata.doi}

    # Add disease
    if metadata.disease:
        props["Disease"] = {
            "multi_select": [{"name": d} for d in metadata.disease if d]
        }

    # Add matrix/sample type
    if metadata.sample_type:
        props["Matrix"] = {
            "multi_select": [{"name": st} for st in metadata.sample_type if st]
        }

    # Add model systems/organism
    if metadata.organism:
        props["Model Systems"] = {
            "multi_select": [{"name": o} for o in metadata.organism if o]
        }

    # Add omics type
    omics_type_map = {
        "transcriptomics": "Transcriptomics",
        "proteomics": "Proteomics",
        "metabolomics": "Metabolomics",
        "lipidomics": "Lipidomics",
    }
    omics_type_value = omics_type_map.get(metadata.omics_type, metadata.omics_type)
    props["Omics Type"] = {"select": {"name": omics_type_value}}

    # Create page
    try:
        cfg = get_config()
        url = f"{cfg.notion.base_url}/pages"
        payload = {
            "parent": {"database_id": exp_data_db_id},
            "properties": props,
        }

        resp = requests.post(
            url,
            headers=notion_headers(),
            json=payload,
            timeout=30,
        )
        resp.raise_for_status()

        page_id = resp.json().get("id", "")
        logger.info(
            "[HARVEST] Created Dataset page %s for %s study %s",
            page_id,
            repository,
            metadata.study_id,
        )
        return page_id

    except Exception as e:
        logger.error(
            "[HARVEST] Error creating Dataset page for %s study %s: %r",
            repository,
            metadata.study_id,
            e,
        )
        return None


def create_postgres_dataset_from_metadata(
    metadata,
    repository: str,
) -> tuple[Optional[str], Optional[str]]:
    """
    Create a Postgres dataset from StudyMetadata.

    Args:
        metadata: StudyMetadata object
        repository: Repository name

    Returns:
        Tuple of (dataset_id, notion_page_id) - dataset_id is Postgres UUID string
    """
    cfg = get_config()

    if not cfg.pipeline.use_postgres_as_sot:
        logger.debug("[HARVEST] Postgres not enabled, skipping Postgres dataset creation")
        return None, None

    try:
        # Map omics type string to OmicsType enum
        omics_type_map = {
            "transcriptomics": OmicsType.TRANSCRIPTOMICS,
            "proteomics": OmicsType.PROTEOMICS,
            "metabolomics": OmicsType.METABOLOMICS,
            "lipidomics": OmicsType.LIPIDOMICS,
        }
        omics_type = omics_type_map.get(metadata.omics_type)
        if not omics_type:
            logger.warning(
                "[HARVEST] Unknown omics type '%s', skipping Postgres creation",
                metadata.omics_type,
            )
            return None, None

        # Build external IDs dict with repository study ID
        external_ids = {
            f"{repository.lower()}_study_id": metadata.study_id,
        }
        if metadata.doi:
            external_ids["doi"] = metadata.doi
        if metadata.pubmed_id:
            external_ids["pubmed_id"] = metadata.pubmed_id

        # Build file URLs from data files
        file_urls = [df.download_url for df in metadata.data_files if df.download_url]

        # Create Postgres dataset
        dataset = create_or_update_dataset_in_postgres(
            name=metadata.title,
            omics_type=omics_type,
            description=metadata.summary,
            file_urls=file_urls if file_urls else None,
            organism=metadata.organism if metadata.organism else None,
            sample_type=metadata.sample_type if metadata.sample_type else None,
            disease=metadata.disease if metadata.disease else None,
            external_ids=external_ids,
            notion_page_id=None,  # Will be linked later if Notion page created
        )

        logger.info(
            "[HARVEST] Created Postgres dataset %s for %s study %s",
            dataset.id,
            repository,
            metadata.study_id,
        )

        return str(dataset.id), dataset.notion_page_id

    except Exception as e:
        logger.error(
            "[HARVEST] Error creating Postgres dataset for %s study %s: %r",
            repository,
            metadata.study_id,
            e,
        )
        return None, None


def harvest_study(
    study_id: str,
    repository: str,
    create_notion: bool = False,
    create_postgres: bool = True,
    ingest: bool = False,
    dry_run: bool = False,
) -> tuple[Optional[str], Optional[str]]:
    """
    Harvest a study from a repository.

    Args:
        study_id: Repository-specific study identifier
        repository: Repository name
        create_notion: If True, create/update Notion page
        create_postgres: If True, create Postgres dataset (default: True if Postgres is SoT)
        ingest: If True, trigger dataset ingestion after creating page
        dry_run: If True, only print what would be done

    Returns:
        Tuple of (dataset_id, notion_page_id) where dataset_id is Postgres UUID or None
    """
    cfg = get_config()

    # Auto-enable Postgres if it's configured as SoT
    if create_postgres is True and cfg.pipeline.use_postgres_as_sot:
        create_postgres = True
    else:
        create_postgres = False

    logger.info(
        "[HARVEST] Harvesting %s study %s (create_notion=%s, create_postgres=%s, ingest=%s)",
        repository,
        study_id,
        create_notion,
        create_postgres,
        ingest,
    )

    # Get repository instance
    # Load API keys and email from config if available
    geo_api_key = None
    ncbi_email = None
    if repository == "GEO":
        try:
            from amprenta_rag.config import GEO_API_KEY, NCBI_EMAIL
            geo_api_key = GEO_API_KEY or None
            ncbi_email = NCBI_EMAIL or None
        except Exception:
            import os
            geo_api_key = os.getenv("GEO_API_KEY", "") or None
            ncbi_email = os.getenv("NCBI_EMAIL", "") or None

    repo_instance = get_repository(repository, api_key=geo_api_key, email=ncbi_email)
    if not repo_instance:
        logger.error("[HARVEST] Repository '%s' not found", repository)
        return None, None

    # Fetch metadata
    metadata = fetch_study_metadata(study_id, repository)
    if not metadata:
        logger.error(
            "[HARVEST] Could not fetch metadata for %s study %s",
            repository,
            study_id,
        )
        return None, None

    if dry_run:
        print(f"\n{'=' * 80}")
        print(f"DRY RUN: Would harvest {repository} study {study_id}")
        print(f"{'=' * 80}")
        print(f"Title: {metadata.title}")
        print(f"Omics Type: {metadata.omics_type}")
        if metadata.disease:
            print(f"Disease: {', '.join(metadata.disease)}")
        if metadata.organism:
            print(f"Organism: {', '.join(metadata.organism)}")
        if metadata.sample_type:
            print(f"Sample Type: {', '.join(metadata.sample_type)}")
        print(f"{'=' * 80}\n")
        return None, None

    # Create Postgres dataset first (if enabled)
    dataset_id = None
    page_id = None

    if create_postgres:
        dataset_id, existing_notion_page_id = create_postgres_dataset_from_metadata(
            metadata,
            repository,
        )
        if existing_notion_page_id:
            page_id = existing_notion_page_id
            logger.info(
                "[HARVEST] Postgres dataset %s already linked to Notion page %s",
                dataset_id,
                page_id,
            )

        # Try to link to programs/experiments if metadata contains relevant information
        # Note: Repository metadata doesn't typically include program/experiment IDs,
        # but we can try to match based on disease or other metadata if needed
        if dataset_id:
            try:
                # For now, we don't link programs/experiments from repository metadata
                # as they typically don't include Notion page IDs
                # This can be enhanced later if repository metadata includes program/experiment references
                pass
            except Exception as e:
                logger.debug(
                    "[HARVEST] Program/experiment linking skipped for repository study: %r",
                    e,
                )

    # Find or create Notion page (only if explicitly requested)
    if create_notion:
        logger.info("[HARVEST] Notion sync enabled - querying Notion for existing pages")
        # Check if page already exists
        if not page_id:
            page_id = find_existing_dataset_page(study_id, repository)

        if not page_id:
            # Create new page
            page_id = create_dataset_page_from_metadata(metadata, repository)

            # Link Postgres dataset to Notion page if both exist
            if dataset_id and page_id:
                try:
                    from amprenta_rag.database.base import get_db
                    from amprenta_rag.database.models import Dataset as DatasetModel
                    from uuid import UUID

                    db = next(get_db())
                    dataset = db.query(DatasetModel).filter(
                        DatasetModel.id == UUID(dataset_id)
                    ).first()
                    if dataset:
                        dataset.notion_page_id = page_id
                        db.commit()
                        logger.info(
                            "[HARVEST] Linked Postgres dataset %s to Notion page %s",
                            dataset_id,
                            page_id,
                        )
                except Exception as e:
                    logger.warning(
                        "[HARVEST] Could not link Postgres dataset to Notion page: %r",
                        e,
                    )
        else:
            logger.info(
                "[HARVEST] Using existing page %s for %s study %s",
                page_id,
                repository,
                study_id,
            )

            # Link Notion page to Postgres dataset if dataset exists
            if dataset_id and page_id:
                try:
                    from amprenta_rag.database.base import get_db
                    from amprenta_rag.database.models import Dataset as DatasetModel
                    from uuid import UUID

                    db = next(get_db())
                    dataset = db.query(DatasetModel).filter(
                        DatasetModel.id == UUID(dataset_id)
                    ).first()
                    if dataset and not dataset.notion_page_id:
                        dataset.notion_page_id = page_id
                        db.commit()
                        logger.info(
                            "[HARVEST] Linked Notion page %s to Postgres dataset %s",
                            page_id,
                            dataset_id,
                        )
                except Exception as e:
                    logger.warning(
                        "[HARVEST] Could not link Notion page to Postgres dataset: %r",
                        e,
                    )

    # For MW studies, add mwTab data to the page
    if repository in ["MW", "MW_LIPIDOMICS", "MW_METABOLOMICS"] and page_id:
        try:
            from scripts.harvest_mw_studies import fetch_mw_mwtab, add_mwtab_block

            mwtab_text = fetch_mw_mwtab(study_id)
            if mwtab_text:
                add_mwtab_block(page_id, mwtab_text)
                logger.info(
                    "[HARVEST] Added mwTab data to page %s",
                    page_id,
                )
        except Exception as e:
            logger.warning(
                "[HARVEST] Could not add mwTab data to page %s: %r",
                page_id,
                e,
            )

    # Trigger ingestion if requested
    if ingest:
        # Use dataset_id if available, otherwise page_id
        identifier = dataset_id if dataset_id else page_id

        if identifier:
            logger.info(
                "[HARVEST] Triggering ingestion for dataset %s",
                identifier,
            )
            try:
                # Prefer Postgres-only ingestion (fast, no Notion API calls)
                if dataset_id:
                    from uuid import UUID
                    logger.info(
                        "[HARVEST] Using Postgres-only ingestion for dataset %s",
                        dataset_id,
                    )
                    ingest_dataset_from_postgres(
                        dataset_id=UUID(dataset_id),
                        force=False,
                        update_notion=False,  # Postgres-only mode
                    )
                elif page_id:
                    # Fallback to Notion-based ingestion if no Postgres dataset
                    logger.info(
                        "[HARVEST] Using Notion-based ingestion for page %s (no Postgres dataset)",
                        page_id,
                    )
                    ingest_dataset(page_id, force=False)
                else:
                    logger.warning(
                        "[HARVEST] Cannot trigger ingestion: no dataset ID or page ID available",
                    )
            except Exception as e:
                logger.error(
                    "[HARVEST] Error during ingestion for dataset %s: %r",
                    identifier,
                    e,
                )
        else:
            logger.warning(
                "[HARVEST] Cannot trigger ingestion: no dataset ID or page ID available",
            )

    return dataset_id, page_id


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Harvest a study from a public repository (GEO, PRIDE, MetaboLights, MW)"
    )
    parser.add_argument(
        "--study-id",
        required=True,
        help="Repository-specific study identifier (e.g., ST001111, GSE12345, PXD012345)",
    )
    parser.add_argument(
        "--repository",
        required=True,
        help="Repository name (GEO, PRIDE, MetaboLights, MW, MW_LIPIDOMICS, MW_METABOLOMICS)",
    )
    parser.add_argument(
        "--create-notion",
        action="store_true",
        help="Create/update Notion Dataset page (optional, Postgres is default)",
    )
    parser.add_argument(
        "--create-postgres",
        action="store_true",
        default=True,
        help="Create Postgres dataset (default: True if Postgres is SoT)",
    )
    parser.add_argument(
        "--no-postgres",
        action="store_true",
        help="Skip Postgres dataset creation",
    )
    parser.add_argument(
        "--ingest",
        action="store_true",
        help="Trigger dataset ingestion after creating page",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Only print what would be done, don't create pages",
    )

    args = parser.parse_args()

    # Allow Postgres-only operation (no Notion required)
    create_postgres = args.create_postgres and not args.no_postgres

    if not args.create_notion and not args.dry_run and not create_postgres:
        parser.error("Either --create-notion, --create-postgres, or --dry-run must be specified")

    try:
        dataset_id, page_id = harvest_study(
            study_id=args.study_id,
            repository=args.repository,
            create_notion=args.create_notion,
            create_postgres=create_postgres,
            ingest=args.ingest,
            dry_run=args.dry_run,
        )

        if dataset_id or page_id:
            print(f"\n✅ Harvest complete!")
            print(f"   Study ID: {args.study_id}")
            print(f"   Repository: {args.repository}")
            if dataset_id:
                print(f"   Postgres Dataset ID: {dataset_id}")
            if page_id:
                print(f"   Notion Page ID: {page_id}")
            if args.ingest:
                print(f"   Ingestion: Triggered")
        elif args.dry_run:
            print("\n✅ Dry run complete (no datasets created)")

    except Exception as e:
        logger.error("[HARVEST] Harvest failed: %r", e)
        print(f"\n❌ Harvest failed: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

