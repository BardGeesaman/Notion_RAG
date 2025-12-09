#!/usr/bin/env python3
"""
Verify Notion database setup for the RAG system.

Checks that all required databases are configured and accessible.
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.config import get_config
from amprenta_rag.clients.notion_client import notion_headers
import requests

def check_database(db_id: str, db_name: str) -> tuple[bool, str]:
    """
    Check if a Notion database is accessible.
    
    Returns:
        (is_accessible, error_message)
    """
    if not db_id:
        return False, "Database ID not configured"
    
    cfg = get_config()
    url = f"{cfg.notion.base_url}/databases/{db_id}"
    
    try:
        resp = requests.get(url, headers=notion_headers(), timeout=10)
        if resp.status_code == 200:
            return True, "✅ Accessible"
        elif resp.status_code == 404:
            return False, "❌ Database not found (check ID)"
        else:
            return False, f"❌ Error: {resp.status_code}"
    except Exception as e:
        return False, f"❌ Error: {e}"


def main():
    print("\n" + "=" * 80)
    print("NOTION DATABASE SETUP VERIFICATION")
    print("=" * 80 + "\n")
    
    cfg = get_config()
    
    # Core databases
    core_databases = [
        ("Experimental Data Assets", cfg.notion.exp_data_db_id),
        ("Lipid Signatures", cfg.notion.signature_db_id),
        ("Lipid Signature Components", cfg.notion.signature_component_db_id),
        ("Lipid Species", cfg.notion.lipid_species_db_id),
        ("Programs", cfg.notion.programs_db_id),
        ("Experiments", cfg.notion.experiments_db_id),
    ]
    
    # Feature databases
    feature_databases = [
        ("Metabolite Features", cfg.notion.metabolite_features_db_id),
        ("Protein Features", cfg.notion.protein_features_db_id),
        ("Gene Features", cfg.notion.gene_features_db_id),
    ]
    
    # Chemistry databases
    chemistry_databases = [
        ("Compound Features", getattr(cfg.notion, "compound_features_db_id", None)),
        ("HTS Campaigns", getattr(cfg.notion, "hts_campaigns_db_id", None)),
        ("Biochemical Hits", getattr(cfg.notion, "biochemical_hits_db_id", None)),
    ]
    
    # Pathway database
    pathway_databases = [
        ("Pathways", getattr(cfg.notion, "pathways_db_id", None)),
    ]
    
    print("📊 CORE DATABASES")
    print("-" * 80)
    all_good = True
    for name, db_id in core_databases:
        accessible, message = check_database(db_id, name)
        status = "✅" if accessible else "⚠️ "
        print(f"{status} {name:40s} {message}")
        if not accessible:
            all_good = False
    
    print("\n🧬 FEATURE DATABASES")
    print("-" * 80)
    for name, db_id in feature_databases:
        accessible, message = check_database(db_id, name)
        status = "✅" if accessible else "⚠️ "
        print(f"{status} {name:40s} {message}")
        if not accessible and db_id:  # Only warn if configured but not accessible
            all_good = False
    
    print("\n🧪 CHEMISTRY DATABASES")
    print("-" * 80)
    for name, db_id in chemistry_databases:
        accessible, message = check_database(db_id, name)
        status = "✅" if accessible else "⚠️ "
        print(f"{status} {name:40s} {message}")
        if not accessible and db_id:  # Only warn if configured but not accessible
            all_good = False
    
    print("\n🧬 PATHWAY DATABASES")
    print("-" * 80)
    for name, db_id in pathway_databases:
        accessible, message = check_database(db_id, name)
        status = "✅" if accessible else "⚠️ "
        print(f"{status} {name:40s} {message}")
        if not accessible and db_id:  # Only warn if configured but not accessible
            all_good = False
    
    print("\n" + "=" * 80)
    if all_good:
        print("✅ All configured databases are accessible!")
    else:
        print("⚠️  Some databases are not configured or not accessible.")
        print("   See docs/NOTION_DATABASE_SETUP.md for setup instructions.")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()

