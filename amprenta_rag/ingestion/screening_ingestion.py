"""
Screening data ingestion pipeline.

Ingests HTS campaign data, hit lists, and biochemical results into SQLite
and optionally creates/updates Notion pages for promoted compounds.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from amprenta_rag.chemistry.database import (
    get_hits_for_campaign,
    insert_compound,
    insert_hts_campaign,
    insert_hts_results,
)
from amprenta_rag.chemistry.normalization import (
    compute_molecular_descriptors,
    generate_compound_id,
    normalize_smiles,
)
from amprenta_rag.chemistry.schema import (
    BiochemicalResult,
    Compound,
    HTSCampaign,
    HTSResult,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def ingest_hts_campaign_metadata(
    metadata_file: Path,
    campaign_id: Optional[str] = None,
) -> HTSCampaign:
    """
    Ingest HTS campaign metadata from a file.
    
    Expected format: CSV/TSV with columns:
    - campaign_name, description, assay_type, target, library_id, total_wells, run_date
    
    Args:
        metadata_file: Path to metadata file
        campaign_id: Optional campaign ID (auto-generated if not provided)
        
    Returns:
        HTSCampaign object
    """
    logger.info(
        "[INGEST][SCREENING] Ingesting HTS campaign metadata from %s",
        metadata_file,
    )
    
    # Read metadata file
    try:
        df = pd.read_csv(metadata_file, sep="\t" if metadata_file.suffix == ".tsv" else ",")
    except Exception as e:
        raise ValueError(f"Error reading metadata file: {e}")
    
    # Extract metadata (assume single row or take first row)
    row = df.iloc[0].to_dict()
    
    if not campaign_id:
        campaign_id = row.get("campaign_id") or row.get("campaign_name", "UNKNOWN")
    
    campaign = HTSCampaign(
        campaign_id=campaign_id,
        campaign_name=row.get("campaign_name", campaign_id),
        description=row.get("description"),
        assay_type=row.get("assay_type"),
        target=row.get("target"),
        library_id=row.get("library_id"),
        total_wells=row.get("total_wells"),
        hit_count=row.get("hit_count", 0),
        run_date=row.get("run_date"),
    )
    
    # Insert into database
    insert_hts_campaign(campaign)
    
    logger.info(
        "[INGEST][SCREENING] Ingested HTS campaign %s: %s",
        campaign_id,
        campaign.campaign_name,
    )
    
    return campaign


def ingest_hts_hit_list(
    hit_list_file: Path,
    campaign_id: str,
    smiles_column: str = "SMILES",
    vendor_id_column: Optional[str] = None,
    value_column: str = "value",
    hit_threshold: Optional[float] = None,
) -> int:
    """
    Ingest HTS hit list from a file.
    
    Expected format: CSV/TSV with SMILES and assay values.
    
    Args:
        hit_list_file: Path to hit list file
        campaign_id: Campaign ID to associate results with
        smiles_column: Name of SMILES column
        vendor_id_column: Optional vendor ID column
        value_column: Name of value/activity column
        hit_threshold: Optional threshold for hit classification
        
    Returns:
        Number of compounds ingested
    """
    logger.info(
        "[INGEST][SCREENING] Ingesting HTS hit list from %s for campaign %s",
        hit_list_file,
        campaign_id,
    )
    
    # Read hit list file
    try:
        df = pd.read_csv(
            hit_list_file,
            sep="\t" if hit_list_file.suffix == ".tsv" else ",",
        )
    except Exception as e:
        raise ValueError(f"Error reading hit list file: {e}")
    
    if smiles_column not in df.columns:
        raise ValueError(f"SMILES column '{smiles_column}' not found in file")
    
    compounds: Dict[str, Compound] = {}
    results: List[HTSResult] = []
    
    for idx, row in df.iterrows():
        smiles = str(row[smiles_column]).strip()
        if not smiles or smiles == "nan":
            continue
        
        # Normalize compound
        canonical_smiles, inchi_key, molecular_formula = normalize_smiles(smiles)
        compound_id = generate_compound_id(smiles)
        
        # Compute descriptors
        descriptors = compute_molecular_descriptors(smiles)
        
        # Create compound object
        compound = Compound(
            compound_id=compound_id,
            smiles=smiles,
            inchi_key=inchi_key,
            canonical_smiles=canonical_smiles,
            molecular_formula=molecular_formula,
            **descriptors,
        )
        compounds[compound_id] = compound
        
        # Extract value
        raw_value = None
        if value_column in df.columns:
            try:
                raw_value = float(row[value_column])
            except (ValueError, TypeError):
                pass
        
        # Determine hit status
        hit_flag = 0
        hit_category = None
        if hit_threshold is not None and raw_value is not None:
            if raw_value >= hit_threshold:
                hit_flag = 1
                hit_category = "hit"
        elif raw_value is not None:
            # Auto-detect: assume higher values are hits
            # This is a simple heuristic - should be configurable
            hit_flag = 1 if raw_value > 0 else 0
            hit_category = "hit" if hit_flag else "inactive"
        
        # Create result
        result_id = f"{campaign_id}_{compound_id}_{idx}"
        vendor_id = row.get(vendor_id_column) if vendor_id_column else None
        
        result = HTSResult(
            result_id=result_id,
            campaign_id=campaign_id,
            compound_id=compound_id,
            well_position=None,  # Not available in hit lists typically
            raw_value=raw_value,
            normalized_value=raw_value,  # Could add normalization logic
            z_score=None,
            hit_flag=hit_flag,
            hit_category=hit_category,
        )
        results.append(result)
    
    # Batch insert compounds
    for compound in compounds.values():
        insert_compound(compound)
    
    # Batch insert results
    insert_hts_results(results)
    
    logger.info(
        "[INGEST][SCREENING] Ingested %d compounds and %d results for campaign %s",
        len(compounds),
        len(results),
        campaign_id,
    )
    
    return len(compounds)


def ingest_biochemical_results(
    results_file: Path,
    campaign_id: Optional[str] = None,
    smiles_column: str = "SMILES",
    assay_name: Optional[str] = None,
) -> int:
    """
    Ingest biochemical assay results from a file.
    
    Expected format: CSV/TSV with SMILES and IC50/EC50/Ki/Kd values.
    
    Args:
        results_file: Path to results file
        campaign_id: Optional campaign ID
        smiles_column: Name of SMILES column
        assay_name: Assay name (or auto-detect from filename)
        
    Returns:
        Number of results ingested
    """
    logger.info(
        "[INGEST][SCREENING] Ingesting biochemical results from %s",
        results_file,
    )
    
    if not assay_name:
        assay_name = results_file.stem
    
    # Read results file
    try:
        df = pd.read_csv(
            results_file,
            sep="\t" if results_file.suffix == ".tsv" else ",",
        )
    except Exception as e:
        raise ValueError(f"Error reading results file: {e}")
    
    if smiles_column not in df.columns:
        raise ValueError(f"SMILES column '{smiles_column}' not found in file")
    
    compounds: Dict[str, Compound] = {}
    biochemical_results: List[BiochemicalResult] = []
    
    for idx, row in df.iterrows():
        smiles = str(row[smiles_column]).strip()
        if not smiles or smiles == "nan":
            continue
        
        # Normalize compound
        canonical_smiles, inchi_key, molecular_formula = normalize_smiles(smiles)
        compound_id = generate_compound_id(smiles)
        
        # Compute descriptors
        descriptors = compute_molecular_descriptors(smiles)
        
        # Create compound object
        compound = Compound(
            compound_id=compound_id,
            smiles=smiles,
            inchi_key=inchi_key,
            canonical_smiles=canonical_smiles,
            molecular_formula=molecular_formula,
            **descriptors,
        )
        compounds[compound_id] = compound
        
        # Extract activity values
        ic50 = row.get("IC50") or row.get("ic50")
        ec50 = row.get("EC50") or row.get("ec50")
        ki = row.get("Ki") or row.get("ki")
        kd = row.get("Kd") or row.get("kd")
        
        # Convert to float if possible
        for val_name, val in [("ic50", ic50), ("ec50", ec50), ("ki", ki), ("kd", kd)]:
            if val is not None:
                try:
                    locals()[val_name] = float(val)
                except (ValueError, TypeError):
                    locals()[val_name] = None
        
        # Determine activity type
        activity_type = None
        if ic50:
            activity_type = "IC50"
        elif ec50:
            activity_type = "EC50"
        elif ki:
            activity_type = "Ki"
        elif kd:
            activity_type = "Kd"
        
        # Create biochemical result
        result_id = f"biochem_{compound_id}_{idx}"
        target = row.get("target") or row.get("Target")
        
        result = BiochemicalResult(
            result_id=result_id,
            compound_id=compound_id,
            assay_name=assay_name,
            target=target,
            ic50=ic50,
            ec50=ec50,
            ki=ki,
            kd=kd,
            activity_type=activity_type,
            units=row.get("units") or row.get("Units", "nM"),
            run_date=row.get("run_date") or row.get("date"),
        )
        biochemical_results.append(result)
    
    # Batch insert compounds
    for compound in compounds.values():
        insert_compound(compound)
    
    # Insert biochemical results
    from amprenta_rag.chemistry.database import insert_biochemical_results
    insert_biochemical_results(biochemical_results)
    
    logger.info(
        "[INGEST][SCREENING] Ingested %d compounds for biochemical assay %s",
        len(compounds),
        assay_name,
    )
    
    return len(biochemical_results)


def promote_compounds_to_notion(
    campaign_id: str,
    program_id: Optional[str] = None,
    min_hit_value: Optional[float] = None,
    max_compounds: Optional[int] = None,
) -> List[str]:
    """
    Promote selected compounds from an HTS campaign to Notion.
    
    Args:
        campaign_id: Campaign ID
        program_id: Optional program ID to link compounds to
        min_hit_value: Minimum value threshold for promotion
        max_compounds: Maximum number of compounds to promote
        
    Returns:
        List of promoted compound Notion page IDs
    """
    logger.info(
        "[INGEST][SCREENING] Promoting compounds from campaign %s to Notion",
        campaign_id,
    )
    
    # Get hits from campaign
    hits = get_hits_for_campaign(campaign_id)
    
    # Filter by threshold if provided
    if min_hit_value is not None:
        hits = [h for h in hits if h.normalized_value and h.normalized_value >= min_hit_value]
    
    # Sort by value (descending)
    hits.sort(key=lambda h: h.normalized_value or 0.0, reverse=True)
    
    # Limit if requested
    if max_compounds:
        hits = hits[:max_compounds]
    
    # Import Notion integration
    from amprenta_rag.chemistry.notion_integration import find_or_create_compound_page
    from amprenta_rag.chemistry.database import get_chemistry_db_path
    import sqlite3
    
    promoted_page_ids = []
    
    # Fetch compound details from database
    db_path = get_chemistry_db_path()
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    
    try:
        for hit in hits:
            # Fetch compound from database
            cursor.execute(
                "SELECT * FROM compounds WHERE compound_id = ?",
                (hit.compound_id,),
            )
            row = cursor.fetchone()
            
            if not row:
                logger.warning(
                    "[INGEST][SCREENING] Compound %s not found in database",
                    hit.compound_id,
                )
                continue
            
            # Create Compound object
            compound = Compound(
                compound_id=row[0],
                smiles=row[1],
                inchi_key=row[2],
                canonical_smiles=row[3],
                molecular_formula=row[4],
                molecular_weight=row[5],
                logp=row[6],
                hbd_count=row[7],
                hba_count=row[8],
                rotatable_bonds=row[9],
            )
            
            # Create/find Notion page
            page_id = find_or_create_compound_page(compound, program_id=program_id)
            if page_id:
                promoted_page_ids.append(page_id)
                logger.info(
                    "[INGEST][SCREENING] Promoted compound %s to Notion page %s (value: %s)",
                    hit.compound_id,
                    page_id,
                    hit.normalized_value,
                )
            else:
                logger.warning(
                    "[INGEST][SCREENING] Could not create Notion page for compound %s",
                    hit.compound_id,
                )
    
    finally:
        conn.close()
    
    logger.info(
        "[INGEST][SCREENING] Promoted %d/%d compounds to Notion",
        len(promoted_page_ids),
        len(hits),
    )
    
    return promoted_page_ids

