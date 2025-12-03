#!/usr/bin/env python3
"""
Create the "Metabolite Features" Notion database.

This script creates a new Notion database with the required schema for metabolite/feature tracking.

Usage:
    python scripts/create_metabolite_features_db.py --parent-page-id <parent_page_id>
    
After creation, the database ID will be printed. Add it to your .env file:
    NOTION_METABOLITE_FEATURES_DB_ID=<database_id>
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import requests
from amprenta_rag.config import get_config
from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def create_metabolite_features_database(parent_page_id: str) -> str:
    """
    Create the Metabolite Features database in Notion.
    
    Args:
        parent_page_id: Notion page ID (with or without dashes) where database will be created
        
    Returns:
        Database ID (without dashes)
    """
    cfg = get_config()
    
    # Normalize parent page ID (remove dashes)
    parent_page_id_clean = parent_page_id.replace("-", "")
    
    url = f"{cfg.notion.base_url}/databases"
    
    properties = {
        "Name": {
            "title": {}
        },
        "Class": {
            "select": {
                "options": [
                    {"name": "Amino Acid", "color": "blue"},
                    {"name": "Organic Acid", "color": "green"},
                    {"name": "Nucleotide", "color": "yellow"},
                    {"name": "Coenzyme", "color": "orange"},
                    {"name": "Sugar", "color": "pink"},
                    {"name": "Lipid", "color": "purple"},
                    {"name": "Other", "color": "gray"},
                ]
            }
        },
        "Synonyms": {
            "rich_text": {}
        },
        "Notes": {
            "rich_text": {}
        },
        "Pathways": {
            "multi_select": {
                "options": []
            }
        },
    }
    
    # Add relations to databases we have IDs for
    if cfg.notion.exp_data_db_id:
        properties["Datasets"] = {
            "relation": {
                "database_id": cfg.notion.exp_data_db_id,
                "type": "dual_property",
            }
        }
    
    if cfg.notion.lit_db_id:
        properties["Literature Mentions"] = {
            "relation": {
                "database_id": cfg.notion.lit_db_id,
                "type": "dual_property",
            }
        }
    
    if cfg.notion.email_db_id:
        properties["Emails / Notes"] = {
            "relation": {
                "database_id": cfg.notion.email_db_id,
                "type": "dual_property",
            }
        }
    
    # Note: Experiments, Lipid Species, and Lipid Signature Components relations
    # should be added manually in Notion UI or via API update once those DB IDs are available
    
    payload = {
        "parent": {
            "type": "page_id",
            "page_id": parent_page_id_clean,
        },
        "title": [
            {
                "type": "text",
                "text": {"content": "Metabolite Features"},
            }
        ],
        "properties": properties,
    }
    
    try:
        resp = requests.post(
            url,
            headers=notion_headers(),
            json=payload,
            timeout=30,
        )
        resp.raise_for_status()
        
        database = resp.json()
        db_id = database.get("id", "")
        
        # Remove dashes for consistency with other DB IDs
        db_id_clean = db_id.replace("-", "")
        
        logger.info(
            "[FEATURES] Created Metabolite Features database: %s (clean: %s)",
            db_id,
            db_id_clean,
        )
        
        return db_id_clean
    except Exception as e:
        logger.error(
            "[FEATURES] Error creating database: %r - Response: %s",
            e,
            e.response.text if hasattr(e, 'response') else "N/A",
        )
        raise


def main():
    parser = argparse.ArgumentParser(
        description="Create the Metabolite Features Notion database."
    )
    parser.add_argument(
        "--parent-page-id",
        required=True,
        help="Notion page ID (with or without dashes) where database will be created",
    )
    args = parser.parse_args()
    
    try:
        db_id = create_metabolite_features_database(args.parent_page_id)
        print("\n" + "=" * 80)
        print("✅ Successfully created Metabolite Features database!")
        print("=" * 80)
        print(f"\nDatabase ID: {db_id}")
        print("\nAdd this to your .env file:")
        print(f"NOTION_METABOLITE_FEATURES_DB_ID={db_id}")
        print("\n" + "=" * 80)
        return 0
    except Exception as e:
        print(f"\n❌ Error creating database: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())

