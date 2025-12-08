from __future__ import annotations

import sqlite3
from typing import List, Optional

from fastapi import HTTPException

from amprenta_rag.chemistry.database import get_chemistry_db_path
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def list_campaigns() -> List[dict]:
    db_path = get_chemistry_db_path()
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    try:
        cursor.execute(
            """
            SELECT campaign_id, campaign_name, description, assay_type,
                   target, library_id, total_wells, hit_count, run_date
            FROM hts_campaigns
            ORDER BY run_date DESC
            """
        )
        rows = cursor.fetchall()
        return [
            {
                "campaign_id": r[0],
                "campaign_name": r[1],
                "description": r[2],
                "assay_type": r[3],
                "target": r[4],
                "library_id": r[5],
                "total_wells": r[6],
                "hit_count": r[7],
                "run_date": r[8],
            }
            for r in rows
        ]
    except sqlite3.Error as e:
        logger.error("[API][SCREENING] Database error listing campaigns: %r", e)
        raise HTTPException(status_code=500, detail=f"Database error: {e}") from e
    finally:
        conn.close()


def get_campaign(campaign_id: str) -> Optional[dict]:
    db_path = get_chemistry_db_path()
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    try:
        cursor.execute(
            """
            SELECT campaign_id, campaign_name, description, assay_type,
                   target, library_id, total_wells, hit_count, run_date
            FROM hts_campaigns
            WHERE campaign_id = ?
            """,
            (campaign_id,),
        )
        r = cursor.fetchone()
        if not r:
            return None
        return {
            "campaign_id": r[0],
            "campaign_name": r[1],
            "description": r[2],
            "assay_type": r[3],
            "target": r[4],
            "library_id": r[5],
            "total_wells": r[6],
            "hit_count": r[7],
            "run_date": r[8],
        }
    except sqlite3.Error as e:
        logger.error("[API][SCREENING] Database error fetching campaign %s: %r", campaign_id, e)
        raise HTTPException(status_code=500, detail=f"Database error: {e}") from e
    finally:
        conn.close()


def get_campaign_hits(campaign_id: str) -> List[dict]:
    db_path = get_chemistry_db_path()
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    try:
        cursor.execute(
            """
            SELECT result_id, compound_id, well_position,
                   raw_value, normalized_value, z_score,
                   hit_flag, hit_category
            FROM hts_results
            WHERE campaign_id = ? AND hit_flag = 1
            ORDER BY normalized_value DESC
            """,
            (campaign_id,),
        )
        rows = cursor.fetchall()
        return [
            {
                "result_id": r[0],
                "compound_id": r[1],
                "well_position": r[2],
                "raw_value": r[3],
                "normalized_value": r[4],
                "z_score": r[5],
                "hit_flag": r[6],
                "hit_category": r[7],
            }
            for r in rows
        ]
    except sqlite3.Error as e:
        logger.error("[API][SCREENING] Database error fetching hits for %s: %r", campaign_id, e)
        raise HTTPException(status_code=500, detail=f"Database error: {e}") from e
    finally:
        conn.close()

