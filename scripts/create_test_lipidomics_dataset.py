#!/usr/bin/env python3
"""
Create a Notion Experimental Data Asset page for the synthetic lipidomics test dataset.

This creates a page that can be ingested using the dataset ingestion pipeline.
"""

import sys
from pathlib import Path
import json

sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.config import get_config
from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.logging_utils import get_logger
import requests

logger = get_logger(__name__)


def create_test_lipidomics_dataset_page() -> str:
    """
    Create a Notion Experimental Data Asset page for the synthetic lipidomics dataset.
    
    Returns:
        Created page ID
    """
    cfg = get_config()
    
    if not cfg.notion.exp_data_db_id:
        raise RuntimeError("Experimental Data Assets DB ID not configured")
    
    # Create minimal mwTab JSON structure for synthetic dataset
    # This allows the ingestion pipeline to extract species correctly
    mwtab_data = {
        "STUDY_ID": "SYNTHETIC-TEST",
        "STUDY": {
            "Study Title": "Synthetic Lipidomics Test Dataset",
            "Study Summary": "Synthetic lipidomics dataset created for signature scoring test. Contains 6 lipid species matching ALS-CSF-Core-6Ceramides signature.",
        },
        "SUBJECT_SAMPLE_FACTORS": {
            "Data": [
                {
                    "Sample ID": "Test1",
                    "Factor Value[Disease]": "ALS",
                    "Factor Value[Matrix]": "CSF",
                },
            ],
        },
        "MS_METABOLITE_DATA": {
            "Data": [
                {
                    "Metabolite": "Cer(d18:1/16:0)",
                    "Test1": 12345,
                },
                {
                    "Metabolite": "Cer(d18:1/18:0)",
                    "Test1": 18234,
                },
                {
                    "Metabolite": "Cer(d18:1/24:0)",
                    "Test1": 9812,
                },
                {
                    "Metabolite": "Cer(d18:1/24:1)",
                    "Test1": 15234,
                },
                {
                    "Metabolite": "SM(d18:1/16:0)",
                    "Test1": 23145,
                },
                {
                    "Metabolite": "SM(d18:1/18:0)",
                    "Test1": 11234,
                },
            ],
        },
    }
    
    mwtab_json_text = json.dumps(mwtab_data, indent=2)
    
    # Create page properties
    props = {
        "Experiment Name": {
            "title": [{"text": {"content": "Test Lipidomics Dataset (Synthetic)"}}],
        },
        "Summary": {
            "rich_text": [
                {
                    "text": {
                        "content": "Synthetic lipidomics dataset created for signature scoring test. Contains 6 lipid species matching ALS-CSF-Core-6Ceramides signature. Real lipidomics dataset ingested for signature scoring test.",
                    },
                },
            ],
        },
        "Dataset Source Type": {
            "select": {"name": "Processed table"},
        },
        "Data Origin": {
            "select": {"name": "External – Open Dataset"},
        },
        "Disease": {
            "multi_select": [{"name": "ALS"}],
        },
        "Matrix": {
            "multi_select": [{"name": "CSF"}],
        },
        "Model Systems": {
            "multi_select": [{"name": "Homo sapiens"}],
        },
    }
    
    # Create the page
    payload = {
        "parent": {"database_id": cfg.notion.exp_data_db_id},
        "properties": props,
    }
    
    try:
        url = f"{cfg.notion.base_url}/pages"
        resp = requests.post(
            url,
            headers=notion_headers(),
            json=payload,
            timeout=30,
        )
        resp.raise_for_status()
        page_id = resp.json()["id"]
        
        logger.info(
            "[DATASET] Created test lipidomics dataset page: %s",
            page_id,
        )
        
        # Add mwTab data as code blocks (matching the format ingestion expects)
        blocks = []
        # Add heading
        blocks.append({
            "object": "block",
            "type": "heading_2",
            "heading_2": {
                "rich_text": [{"type": "text", "text": {"content": "mwTab Data"}}],
            },
        })
        
        # Add mwTab JSON as code block
        # Split into chunks if needed (Notion code blocks have size limits)
        chunk_size = 1900
        chunks = [
            mwtab_json_text[i:i + chunk_size]
            for i in range(0, len(mwtab_json_text), chunk_size)
        ]
        
        for chunk in chunks:
            blocks.append({
                "object": "block",
                "type": "code",
                "code": {
                    "caption": [],
                    "rich_text": [{"type": "text", "text": {"content": chunk}}],
                    "language": "json",
                },
            })
        
        # Append blocks to page
        url = f"{cfg.notion.base_url}/blocks/{page_id}/children"
        resp = requests.patch(
            url,
            headers=notion_headers(),
            json={"children": blocks},
            timeout=30,
        )
        resp.raise_for_status()
        
        logger.info(
            "[DATASET] Added mwTab data to page %s (%d blocks)",
            page_id,
            len(blocks),
        )
        
        return page_id
    
    except Exception as e:
        logger.error(
            "[DATASET] Error creating test dataset page: %r",
            e,
        )
        raise


if __name__ == "__main__":
    try:
        page_id = create_test_lipidomics_dataset_page()
        print(f"\n✅ Created dataset page: {page_id}\n")
        print(f"To ingest, run:")
        print(f"  python scripts/ingest_dataset.py --dataset-page-id {page_id} --force")
    except Exception as e:
        print(f"\n❌ Error: {e}\n", file=sys.stderr)
        sys.exit(1)
