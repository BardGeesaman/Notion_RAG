#!/usr/bin/env python3
"""
Test chemistry Notion integration.

Tests creating a compound page in Notion to verify database access.
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from amprenta_rag.chemistry.notion_integration import create_compound_feature_page
from amprenta_rag.chemistry.schema import Compound
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def main():
    print("\n" + "=" * 80)
    print("CHEMISTRY NOTION INTEGRATION TEST")
    print("=" * 80 + "\n")
    
    # Create a test compound
    test_compound = Compound(
        compound_id="TEST-001",
        smiles="CCO",  # Ethanol
        inchi_key="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        canonical_smiles="CCO",
        molecular_formula="C2H6O",
        molecular_weight=46.07,
        logp=-0.31,
        hbd_count=1,
        hba_count=1,
        rotatable_bonds=0,
    )
    
    print(f"Test Compound: {test_compound.compound_id}")
    print(f"SMILES: {test_compound.smiles}")
    print(f"InChI Key: {test_compound.inchi_key}\n")
    
    print("Attempting to create Notion page...")
    page_id = create_compound_feature_page(test_compound)
    
    if page_id:
        print(f"✅ Successfully created Notion page: {page_id}")
        print(f"   View at: https://www.notion.so/{page_id.replace('-', '')}")
    else:
        print("❌ Failed to create Notion page")
        print("   Check:")
        print("   1. NOTION_COMPOUND_FEATURES_DB_ID is set in .env")
        print("   2. Database is shared with your integration")
        print("   3. Integration has write permissions")
    
    print("\n" + "=" * 80 + "\n")


if __name__ == "__main__":
    main()

