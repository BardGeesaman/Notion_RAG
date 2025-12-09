#!/usr/bin/env python3
"""
List existing lipidomics datasets in Notion to help identify datasets for re-ingestion.
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import requests
from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config

def list_lipidomics_datasets():
    """List lipidomics datasets from Experimental Data Assets database."""
    cfg = get_config()
    db_id = cfg.notion.exp_data_db_id

    if not db_id:
        print("❌ Experimental Data Assets DB ID not configured")
        return

    # Query for datasets
    url = f"{cfg.notion.base_url}/databases/{db_id}/query"
    resp = requests.post(
        url,
        headers=notion_headers(),
        json={"page_size": 50, "sorts": [{"property": "Created", "direction": "descending"}]},
        timeout=30,
    )
    resp.raise_for_status()
    results = resp.json().get("results", [])

    print(f"=== Found {len(results)} Experimental Data Asset pages ===\n")

    lipidomics_datasets = []

    for page in results:
        props = page.get("properties", {})

        # Get page title/name
        name_prop = props.get("Name", {}).get("title", [])
        page_name = name_prop[0]["plain_text"] if name_prop else "Untitled"

        # Check if it's a lipidomics dataset
        data_origin = props.get("Data Origin", {}).get("rich_text", [])
        origin_text = data_origin[0]["plain_text"] if data_origin else ""

        omics_type = props.get("Omics Type", {})
        omics_value = None
        if omics_type.get("type") == "select":
            omics_value = omics_type.get("select", {}).get("name")
        elif omics_type.get("type") == "multi_select":
            omics_values = omics_type.get("multi_select", [])
            omics_value = ", ".join([v.get("name") for v in omics_values])

        # Check summary for lipidomics keywords
        summary_prop = props.get("Summary", {}).get("rich_text", [])
        summary_text = "".join([rt.get("plain_text", "") for rt in summary_prop])

        is_lipidomics = (
            "lipid" in page_name.lower()
            or "lipidomics" in summary_text.lower()
            or (omics_value and "lipid" in str(omics_value).lower())
            or "Internal" in origin_text
        )

        page_id = page["id"]
        print(f"📄 {page_name}")
        print(f"   Page ID: {page_id}")
        print(f"   Data Origin: {origin_text or 'N/A'}")
        print(f"   Omics Type: {omics_value or 'N/A'}")
        
        if is_lipidomics:
            lipidomics_datasets.append({"name": page_name, "page_id": page_id, "origin": origin_text})
        print()

    if lipidomics_datasets:
        print(f"\n=== Found {len(lipidomics_datasets)} lipidomics datasets ===")
        print("\nTo re-ingest any of these, use:")
        print("\npython scripts/ingest_lipidomics.py \\")
        print("  --file <path_to_your_csv_file> \\")
        print("  --dataset-page-id <page_id_from_above>\n")


if __name__ == "__main__":
    list_lipidomics_datasets()

