"""Chemistry procurement service for searching vendor catalogs."""

from __future__ import annotations

from typing import List, Dict

VENDORS = {
    "sigma_aldrich": {"name": "Sigma-Aldrich", "url": "https://www.sigmaaldrich.com"},
    "fisher": {"name": "Fisher Scientific", "url": "https://www.fishersci.com"},
    "enamine": {"name": "Enamine", "url": "https://enamine.net"},
    "molport": {"name": "MolPort", "url": "https://www.molport.com"},
    "mcule": {"name": "Mcule", "url": "https://mcule.com"},
    "chembridge": {"name": "ChemBridge", "url": "https://www.chembridge.com"},
    "combi_blocks": {"name": "Combi-Blocks", "url": "https://www.combi-blocks.com"},
}


def search_vendors(query: str) -> List[Dict[str, str]]:
    """
    Search vendors for compound. Returns mock data (real APIs need subscriptions).

    Args:
        query: Search query (SMILES, CAS number, compound name, etc.)

    Returns:
        List of vendor results with catalog_id, price, availability, url
    """
    # Return mock results showing structure
    # In production, this would call vendor APIs with proper authentication
    return [
        {
            "vendor": "Sigma-Aldrich",
            "catalog_id": "S1234",
            "price": "$45.00",
            "availability": "In Stock",
            "url": "https://sigmaaldrich.com/...",
        },
        {
            "vendor": "Enamine",
            "catalog_id": "EN300-12345",
            "price": "$25.00",
            "availability": "2-3 weeks",
            "url": "https://enamine.net/...",
        },
        {
            "vendor": "MolPort",
            "catalog_id": "MP-123456",
            "price": "$35.00",
            "availability": "In Stock",
            "url": "https://molport.com/...",
        },
    ]


def get_vendor_info() -> List[Dict[str, str]]:
    """
    Return list of supported vendors.

    Returns:
        List of vendor dictionaries with id, name, and url
    """
    return [{"id": k, **v} for k, v in VENDORS.items()]

