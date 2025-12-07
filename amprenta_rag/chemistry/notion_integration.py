"""
Notion integration for chemistry and HTS data.

Creates and updates Notion pages for promoted compounds, HTS campaigns,
and biochemical hits.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import requests

# DEPRECATED: Notion imports removed - Postgres is now source of truth
# from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config
from amprenta_rag.chemistry.database import get_chemistry_db_path
import sqlite3
from amprenta_rag.chemistry.schema import Compound, HTSCampaign, BiochemicalResult
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

def notion_headers() -> Dict[str, str]:
    """DEPRECATED: Notion support removed. Returns empty headers dict."""
    logger.debug("[CHEMISTRY][NOTION-INTEGRATION] notion_headers() deprecated - Notion support removed")
    return {}


def create_compound_feature_page(
    compound: Compound,
    program_id: Optional[str] = None,
) -> Optional[str]:
    """
    Create a Compound Features page in Notion.
    
    Args:
        compound: Compound object
        program_id: Optional program ID to link to
        
    Returns:
        Notion page ID if created, None otherwise
    """
    cfg = get_config()
    
    # Check if Compound Features DB is configured
    if not hasattr(cfg.notion, "compound_features_db_id") or not cfg.notion.compound_features_db_id:
        logger.warning(
            "[CHEMISTRY][NOTION] Compound Features DB ID not configured. "
            "Set NOTION_COMPOUND_FEATURES_DB_ID in .env"
        )
        return None
    
    db_id = cfg.notion.compound_features_db_id
    
    # Build page properties
    props: Dict[str, Any] = {
        "Name": {"title": [{"text": {"content": compound.compound_id}}]},
        "SMILES": {"rich_text": [{"text": {"content": compound.smiles}}]},
    }
    
    if compound.canonical_smiles:
        props["Canonical SMILES"] = {"rich_text": [{"text": {"content": compound.canonical_smiles}}]}
    
    if compound.inchi_key:
        props["InChI Key"] = {"rich_text": [{"text": {"content": compound.inchi_key}}]}
    
    if compound.molecular_formula:
        props["Molecular Formula"] = {"rich_text": [{"text": {"content": compound.molecular_formula}}]}
    
    if compound.molecular_weight:
        props["Molecular Weight"] = {"number": compound.molecular_weight}
    
    if compound.logp is not None:
        props["LogP"] = {"number": compound.logp}
    
    if compound.hbd_count is not None:
        props["HBD Count"] = {"number": compound.hbd_count}
    
    if compound.hba_count is not None:
        props["HBA Count"] = {"number": compound.hba_count}
    
    if compound.rotatable_bonds is not None:
        props["Rotatable Bonds"] = {"number": compound.rotatable_bonds}
    
    # Link to program if provided
    if program_id:
        props["Related Programs"] = {"relation": [{"id": program_id}]}
    
    try:
        url = f"{cfg.notion.base_url}/pages"
        payload = {
            "parent": {"database_id": db_id},
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
            "[CHEMISTRY][NOTION] Created Compound Features page %s for compound %s",
            page_id,
            compound.compound_id,
        )
        return page_id
        
    except Exception as e:
        logger.error(
            "[CHEMISTRY][NOTION] Error creating Compound Features page for %s: %r",
            compound.compound_id,
            e,
        )
        return None


def find_or_create_compound_page(
    compound: Compound,
    program_id: Optional[str] = None,
) -> Optional[str]:
    """
    Find existing or create new Compound Features page.
    
    Args:
        compound: Compound object
        program_id: Optional program ID to link to
        
    Returns:
        Notion page ID if found/created, None otherwise
    """
    cfg = get_config()
    
    if not hasattr(cfg.notion, "compound_features_db_id") or not cfg.notion.compound_features_db_id:
        return None
    
    db_id = cfg.notion.compound_features_db_id
    
    # Try to find existing page by InChI key or compound ID
    try:
        url = f"{cfg.notion.base_url}/databases/{db_id}/query"
        
        # Search by InChI key if available
        if compound.inchi_key:
            payload = {
                "filter": {
                    "property": "InChI Key",
                    "rich_text": {"equals": compound.inchi_key},
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
                logger.debug(
                    "[CHEMISTRY][NOTION] Found existing Compound Features page %s for compound %s",
                    page_id,
                    compound.compound_id,
                )
                return page_id
        
        # Search by compound ID
        payload = {
            "filter": {
                "property": "Name",
                "title": {"equals": compound.compound_id},
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
            logger.debug(
                "[CHEMISTRY][NOTION] Found existing Compound Features page %s for compound %s",
                page_id,
                compound.compound_id,
            )
            return page_id
        
    except Exception as e:
        logger.debug(
            "[CHEMISTRY][NOTION] Error searching for compound page: %r",
            e,
        )
    
    # Create new page
    return create_compound_feature_page(compound, program_id=program_id)


def create_hts_campaign_page(
    campaign: HTSCampaign,
    program_id: Optional[str] = None,
) -> Optional[str]:
    """
    Create an HTS Campaigns page in Notion.
    
    Args:
        campaign: HTSCampaign object
        program_id: Optional program ID to link to
        
    Returns:
        Notion page ID if created, None otherwise
    """
    cfg = get_config()
    
    # Check if HTS Campaigns DB is configured
    if not hasattr(cfg.notion, "hts_campaigns_db_id") or not cfg.notion.hts_campaigns_db_id:
        logger.warning(
            "[CHEMISTRY][NOTION] HTS Campaigns DB ID not configured. "
            "Set NOTION_HTS_CAMPAIGNS_DB_ID in .env"
        )
        return None
    
    db_id = cfg.notion.hts_campaigns_db_id
    
    # Build page properties
    props: Dict[str, Any] = {
        "Campaign Name": {"title": [{"text": {"content": campaign.campaign_name}}]},
        "Campaign ID": {"rich_text": [{"text": {"content": campaign.campaign_id}}]},
    }
    
    if campaign.description:
        props["Description"] = {"rich_text": [{"text": {"content": campaign.description}}]}
    
    if campaign.assay_type:
        props["Assay Type"] = {"select": {"name": campaign.assay_type}}
    
    if campaign.target:
        props["Target"] = {"rich_text": [{"text": {"content": campaign.target}}]}
    
    if campaign.total_wells:
        props["Total Wells"] = {"number": campaign.total_wells}
    
    if campaign.hit_count is not None:
        props["Hit Count"] = {"number": campaign.hit_count}
    
    if campaign.run_date:
        props["Run Date"] = {"date": {"start": campaign.run_date}}
    
    # Link to program if provided
    if program_id:
        props["Related Programs"] = {"relation": [{"id": program_id}]}
    
    try:
        url = f"{cfg.notion.base_url}/pages"
        payload = {
            "parent": {"database_id": db_id},
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
            "[CHEMISTRY][NOTION] Created HTS Campaigns page %s for campaign %s",
            page_id,
            campaign.campaign_id,
        )
        return page_id
        
    except Exception as e:
        logger.error(
            "[CHEMISTRY][NOTION] Error creating HTS Campaigns page for %s: %r",
            campaign.campaign_id,
            e,
        )
        return None


def find_or_create_hts_campaign_page(
    campaign: HTSCampaign,
    program_id: Optional[str] = None,
) -> Optional[str]:
    """
    Find existing or create new HTS Campaigns page.
    
    Args:
        campaign: HTSCampaign object
        program_id: Optional program ID to link to
        
    Returns:
        Notion page ID if found/created, None otherwise
    """
    cfg = get_config()
    
    if not hasattr(cfg.notion, "hts_campaigns_db_id") or not cfg.notion.hts_campaigns_db_id:
        return None
    
    db_id = cfg.notion.hts_campaigns_db_id
    
    # Try to find existing page
    try:
        url = f"{cfg.notion.base_url}/databases/{db_id}/query"
        payload = {
            "filter": {
                "property": "Campaign ID",
                "rich_text": {"equals": campaign.campaign_id},
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
            logger.debug(
                "[CHEMISTRY][NOTION] Found existing HTS Campaigns page %s for campaign %s",
                page_id,
                campaign.campaign_id,
            )
            return page_id
        
    except Exception as e:
        logger.debug(
            "[CHEMISTRY][NOTION] Error searching for campaign page: %r",
            e,
        )
    
    # Create new page
    return create_hts_campaign_page(campaign, program_id=program_id)


def create_biochemical_hit_page(
    result: BiochemicalResult,
    compound: Compound,
    program_id: Optional[str] = None,
) -> Optional[str]:
    """
    Create a Biochemical Hits page in Notion.
    
    Args:
        result: BiochemicalResult object
        compound: Associated Compound object
        program_id: Optional program ID to link to
        
    Returns:
        Notion page ID if created, None otherwise
    """
    cfg = get_config()
    
    # Check if Biochemical Hits DB is configured
    if not hasattr(cfg.notion, "biochemical_hits_db_id") or not cfg.notion.biochemical_hits_db_id:
        logger.warning(
            "[CHEMISTRY][NOTION] Biochemical Hits DB ID not configured. "
            "Set NOTION_BIOCHEMICAL_HITS_DB_ID in .env"
        )
        return None
    
    db_id = cfg.notion.biochemical_hits_db_id
    
    # Build page properties
    props: Dict[str, Any] = {
        "Assay Name": {"title": [{"text": {"content": result.assay_name}}]},
        "Result ID": {"rich_text": [{"text": {"content": result.result_id}}]},
    }
    
    # Link to compound if we have a Notion page for it
    compound_page_id = find_or_create_compound_page(compound, program_id=program_id)
    if compound_page_id:
        props["Compound"] = {"relation": [{"id": compound_page_id}]}
    
    if result.target:
        props["Target"] = {"rich_text": [{"text": {"content": result.target}}]}
    
    if result.ic50 is not None:
        props["IC50"] = {"number": result.ic50}
    
    if result.ec50 is not None:
        props["EC50"] = {"number": result.ec50}
    
    if result.ki is not None:
        props["Ki"] = {"number": result.ki}
    
    if result.kd is not None:
        props["Kd"] = {"number": result.kd}
    
    if result.activity_type:
        props["Activity Type"] = {"select": {"name": result.activity_type}}
    
    if result.units:
        props["Units"] = {"rich_text": [{"text": {"content": result.units}}]}
    
    if result.run_date:
        props["Run Date"] = {"date": {"start": result.run_date}}
    
    # Link to program if provided
    if program_id:
        props["Related Programs"] = {"relation": [{"id": program_id}]}
    
    try:
        url = f"{cfg.notion.base_url}/pages"
        payload = {
            "parent": {"database_id": db_id},
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
            "[CHEMISTRY][NOTION] Created Biochemical Hits page %s for result %s",
            page_id,
            result.result_id,
        )
        return page_id
        
    except Exception as e:
        logger.error(
            "[CHEMISTRY][NOTION] Error creating Biochemical Hits page for %s: %r",
            result.result_id,
            e,
        )
        return None

