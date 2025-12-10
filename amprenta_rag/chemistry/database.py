"""
SQLite database operations for chemistry and HTS data.

Provides functions to initialize the database and perform CRUD operations.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import List, Optional

from amprenta_rag.chemistry.schema import (
    BIOCHEMICAL_RESULTS_TABLE,
    COMPOUND_PROGRAM_TABLE,
    COMPOUNDS_TABLE,
    COMPOUND_SIGNATURE_TABLE,
    HTS_CAMPAIGNS_TABLE,
    HTS_RESULTS_TABLE,
    INDEXES,
    LIBRARIES_TABLE,
    LIBRARY_COMPOUNDS_TABLE,
    BiochemicalResult,
    CompoundSignatureLink,
    Compound,
    HTSCampaign,
    HTSResult,
)
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def get_chemistry_db_path() -> Path:
    """
    Get the path to the chemistry SQLite database.
    
    Returns:
        Path to the database file
    """
    get_config()
    # Use a default path if not configured
    db_dir = Path.cwd() / "data" / "chemistry"
    db_dir.mkdir(parents=True, exist_ok=True)
    return db_dir / "chemistry.db"


def initialize_chemistry_database(db_path: Optional[Path] = None) -> None:
    """
    Initialize the chemistry SQLite database with all tables and indexes.
    
    Args:
        db_path: Optional path to database file (defaults to get_chemistry_db_path())
    """
    if db_path is None:
        db_path = get_chemistry_db_path()
    
    logger.info(
        "[CHEMISTRY][DB] Initializing chemistry database at %s",
        db_path,
    )
    
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    
    try:
        # Create tables
        cursor.execute(COMPOUNDS_TABLE)
        cursor.execute(LIBRARIES_TABLE)
        cursor.execute(LIBRARY_COMPOUNDS_TABLE)
        cursor.execute(HTS_CAMPAIGNS_TABLE)
        cursor.execute(HTS_RESULTS_TABLE)
        cursor.execute(BIOCHEMICAL_RESULTS_TABLE)
        cursor.execute(COMPOUND_PROGRAM_TABLE)
        cursor.execute(COMPOUND_SIGNATURE_TABLE)
        
        # Create indexes
        for index_sql in INDEXES:
            cursor.execute(index_sql)
        
        conn.commit()
        logger.info("[CHEMISTRY][DB] Database initialized successfully")
        
    except Exception as e:
        conn.rollback()
        logger.error("[CHEMISTRY][DB] Error initializing database: %r", e)
        raise
    finally:
        conn.close()


def insert_compound(compound: Compound, db_path: Optional[Path] = None) -> None:
    """
    Insert or update a compound in the database.
    
    Args:
        compound: Compound object
        db_path: Optional path to database file
    """
    if db_path is None:
        db_path = get_chemistry_db_path()
    
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    
    try:
        cursor.execute(
            """
            INSERT OR REPLACE INTO compounds (
                compound_id, smiles, inchi_key, canonical_smiles,
                molecular_formula, molecular_weight, logp,
                hbd_count, hba_count, rotatable_bonds, updated_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP)
            """,
            (
                compound.compound_id,
                compound.smiles,
                compound.inchi_key,
                compound.canonical_smiles,
                compound.molecular_formula,
                compound.molecular_weight,
                compound.logp,
                compound.hbd_count,
                compound.hba_count,
                compound.rotatable_bonds,
            ),
        )
        conn.commit()
        logger.debug(
            "[CHEMISTRY][DB] Inserted/updated compound %s",
            compound.compound_id,
        )
    except Exception as e:
        conn.rollback()
        logger.error(
            "[CHEMISTRY][DB] Error inserting compound %s: %r",
            compound.compound_id,
            e,
        )
        raise
    finally:
        conn.close()


def find_compound_by_inchi_key(
    inchi_key: str,
    db_path: Optional[Path] = None,
) -> Optional[Compound]:
    """
    Find a compound by InChI key.
    
    Args:
        inchi_key: InChI key to search for
        db_path: Optional path to database file
        
    Returns:
        Compound object or None if not found
    """
    if db_path is None:
        db_path = get_chemistry_db_path()
    
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    
    try:
        cursor.execute(
            "SELECT * FROM compounds WHERE inchi_key = ?",
            (inchi_key,),
        )
        row = cursor.fetchone()
        
        if row:
            return Compound(
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
        return None
    finally:
        conn.close()


def insert_hts_campaign(
    campaign: HTSCampaign,
    db_path: Optional[Path] = None,
) -> None:
    """
    Insert or update an HTS campaign.
    
    Args:
        campaign: HTSCampaign object
        db_path: Optional path to database file
    """
    if db_path is None:
        db_path = get_chemistry_db_path()
    
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    
    try:
        cursor.execute(
            """
            INSERT OR REPLACE INTO hts_campaigns (
                campaign_id, campaign_name, description, assay_type,
                target, library_id, total_wells, hit_count, run_date
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                campaign.campaign_id,
                campaign.campaign_name,
                campaign.description,
                campaign.assay_type,
                campaign.target,
                campaign.library_id,
                campaign.total_wells,
                campaign.hit_count,
                campaign.run_date,
            ),
        )
        conn.commit()
        logger.debug(
            "[CHEMISTRY][DB] Inserted/updated HTS campaign %s",
            campaign.campaign_id,
        )
    except Exception as e:
        conn.rollback()
        logger.error(
            "[CHEMISTRY][DB] Error inserting HTS campaign %s: %r",
            campaign.campaign_id,
            e,
        )
        raise
    finally:
        conn.close()


def insert_hts_results(
    results: List[HTSResult],
    db_path: Optional[Path] = None,
) -> None:
    """
    Insert HTS results in batch.
    
    Args:
        results: List of HTSResult objects
        db_path: Optional path to database file
    """
    if not results:
        return
    
    if db_path is None:
        db_path = get_chemistry_db_path()
    
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    
    try:
        cursor.executemany(
            """
            INSERT OR REPLACE INTO hts_results (
                result_id, campaign_id, compound_id, well_position,
                raw_value, normalized_value, z_score, hit_flag, hit_category
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            [
                (
                    r.result_id,
                    r.campaign_id,
                    r.compound_id,
                    r.well_position,
                    r.raw_value,
                    r.normalized_value,
                    r.z_score,
                    r.hit_flag,
                    r.hit_category,
                )
                for r in results
            ],
        )
        conn.commit()
        logger.info(
            "[CHEMISTRY][DB] Inserted %d HTS results",
            len(results),
        )
    except Exception as e:
        conn.rollback()
        logger.error(
            "[CHEMISTRY][DB] Error inserting HTS results: %r",
            e,
        )
        raise
    finally:
        conn.close()


def insert_biochemical_results(
    results: List[BiochemicalResult],
    db_path: Optional[Path] = None,
) -> None:
    """
    Insert biochemical results in batch.
    
    Args:
        results: List of BiochemicalResult objects
        db_path: Optional path to database file
    """
    if not results:
        return
    
    if db_path is None:
        db_path = get_chemistry_db_path()
    
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    
    try:
        cursor.executemany(
            """
            INSERT OR REPLACE INTO biochemical_results (
                result_id, compound_id, assay_name, target,
                ic50, ec50, ki, kd, activity_type, units, run_date
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            [
                (
                    r.result_id,
                    r.compound_id,
                    r.assay_name,
                    r.target,
                    r.ic50,
                    r.ec50,
                    r.ki,
                    r.kd,
                    r.activity_type,
                    r.units,
                    r.run_date,
                )
                for r in results
            ],
        )
        conn.commit()
        logger.info(
            "[CHEMISTRY][DB] Inserted %d biochemical results",
            len(results),
        )
    except Exception as e:
        conn.rollback()
        logger.error(
            "[CHEMISTRY][DB] Error inserting biochemical results: %r",
            e,
        )
        raise
    finally:
        conn.close()


def get_hits_for_campaign(
    campaign_id: str,
    db_path: Optional[Path] = None,
) -> List[HTSResult]:
    """
    Get all hits for an HTS campaign.
    
    Args:
        campaign_id: Campaign ID
        db_path: Optional path to database file
        
    Returns:
        List of HTSResult objects (hits only)
    """
    if db_path is None:
        db_path = get_chemistry_db_path()
    
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    
    try:
        cursor.execute(
            """
            SELECT * FROM hts_results
            WHERE campaign_id = ? AND hit_flag = 1
            ORDER BY normalized_value DESC
            """,
            (campaign_id,),
        )
        
        rows = cursor.fetchall()
        results = []
        for row in rows:
            results.append(
                HTSResult(
                    result_id=row[0],
                    campaign_id=row[1],
                    compound_id=row[2],
                    well_position=row[3],
                    raw_value=row[4],
                    normalized_value=row[5],
                    z_score=row[6],
                    hit_flag=row[7],
                    hit_category=row[8],
                )
            )
        return results
    finally:
        conn.close()


def insert_compound_signature_link(
    link: CompoundSignatureLink,
    db_path: Optional[Path] = None,
) -> None:
    """
    Insert or update a compound–signature link.

    Args:
        link: CompoundSignatureLink object
        db_path: Optional path to database file
    """
    if db_path is None:
        db_path = get_chemistry_db_path()

    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()

    try:
        cursor.execute(
            """
            INSERT OR REPLACE INTO compound_signature (
                compound_id,
                signature_id,
                effect_type,
                correlation,
                p_value,
                evidence_source
            ) VALUES (?, ?, ?, ?, ?, ?)
            """,
            (
                link.compound_id,
                link.signature_id,
                link.effect_type,
                link.correlation,
                link.p_value,
                link.evidence_source,
            ),
        )
        conn.commit()
        logger.debug(
            "[CHEMISTRY][DB] Inserted/updated compound_signature link %s → %s",
            link.compound_id,
            link.signature_id,
        )
    except Exception as e:
        conn.rollback()
        logger.error(
            "[CHEMISTRY][DB] Error inserting compound_signature link %s → %s: %r",
            link.compound_id,
            link.signature_id,
            e,
        )
        raise
    finally:
        conn.close()


def get_compounds_for_signature(
    signature_id: str,
    db_path: Optional[Path] = None,
) -> List[CompoundSignatureLink]:
    """
    Get all compound–signature links for a given signature_id.

    Args:
        signature_id: Postgres signature UUID as string
        db_path: Optional path to database file

    Returns:
        List of CompoundSignatureLink objects
    """
    if db_path is None:
        db_path = get_chemistry_db_path()

    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()

    try:
        cursor.execute(
            """
            SELECT compound_id,
                   signature_id,
                   effect_type,
                   correlation,
                   p_value,
                   evidence_source,
                   created_at
            FROM compound_signature
            WHERE signature_id = ?
            ORDER BY correlation DESC
            """,
            (signature_id,),
        )
        rows = cursor.fetchall()

        return [
            CompoundSignatureLink(
                compound_id=row[0],
                signature_id=row[1],
                effect_type=row[2],
                correlation=row[3],
                p_value=row[4],
                evidence_source=row[5],
                created_at=row[6],
            )
            for row in rows
        ]
    finally:
        conn.close()


def get_signatures_for_compound(
    compound_id: str,
    db_path: Optional[Path] = None,
) -> List[CompoundSignatureLink]:
    """
    Get all compound–signature links for a given compound_id.

    Args:
        compound_id: Compound identifier used in the chemistry DB
        db_path: Optional path to database file

    Returns:
        List of CompoundSignatureLink objects
    """
    if db_path is None:
        db_path = get_chemistry_db_path()

    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()

    try:
        cursor.execute(
            """
            SELECT compound_id,
                   signature_id,
                   effect_type,
                   correlation,
                   p_value,
                   evidence_source,
                   created_at
            FROM compound_signature
            WHERE compound_id = ?
            ORDER BY correlation DESC
            """,
            (compound_id,),
        )
        rows = cursor.fetchall()

        return [
            CompoundSignatureLink(
                compound_id=row[0],
                signature_id=row[1],
                effect_type=row[2],
                correlation=row[3],
                p_value=row[4],
                evidence_source=row[5],
                created_at=row[6],
            )
            for row in rows
        ]
    finally:
        conn.close()


