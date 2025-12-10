from __future__ import annotations

from typing import Dict, List, Optional

import sqlite3

from amprenta_rag.chemistry.database import get_chemistry_db_path
from amprenta_rag.chemistry.schema import Compound, HTSCampaign
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def _get_conn(db_path: Optional[str] = None) -> sqlite3.Connection:
    path = db_path or get_chemistry_db_path()
    return sqlite3.connect(str(path))


def _fetch_compound(compound_id: str, conn: sqlite3.Connection) -> Optional[Compound]:
    row = conn.execute(
        "SELECT compound_id, smiles, molecular_weight FROM compounds WHERE compound_id = ?",
        (compound_id,),
    ).fetchone()
    if not row:
        return None
    return Compound(
        compound_id=row[0],
        smiles=row[1],
        molecular_weight=row[2],
    )


def _fetch_program_links(compound_id: str, conn: sqlite3.Connection) -> List[str]:
    rows = conn.execute(
        "SELECT program_id FROM compound_program WHERE compound_id = ?",
        (compound_id,),
    ).fetchall()
    return [r[0] for r in rows]


def _fetch_signature_links(compound_id: str, conn: sqlite3.Connection):
    rows = conn.execute(
        """
        SELECT signature_id, effect_type, correlation, p_value
        FROM compound_signature
        WHERE compound_id = ?
        """,
        (compound_id,),
    ).fetchall()
    return rows


def generate_compound_summary(compound_id: str, db_path: Optional[str] = None) -> str:
    conn = _get_conn(db_path)
    try:
        compound = _fetch_compound(compound_id, conn)
        if not compound:
            return f"Compound {compound_id} not found."

        programs = _fetch_program_links(compound_id, conn)
        sig_links = _fetch_signature_links(compound_id, conn)

        summary_parts = [
            f"Compound {compound.compound_id} (SMILES: {compound.smiles}, MW: {compound.molecular_weight})"
        ]
        if programs:
            summary_parts.append(f"linked to programs: {', '.join(programs)}")
        if sig_links:
            sig_strs = []
            for sid, eff, corr, pval in sig_links:
                sig_desc = f"{sid}"
                if eff:
                    sig_desc += f" [{eff}]"
                if corr is not None:
                    sig_desc += f" corr={corr:.3f}"
                if pval is not None:
                    sig_desc += f" p={pval:.2e}"
                sig_strs.append(sig_desc)
            summary_parts.append("signature links: " + "; ".join(sig_strs))

        summary = " | ".join(summary_parts)
        logger.info("[CHEMISTRY][RAG] Generated compound summary for %s", compound_id)
        return summary
    finally:
        conn.close()


def _fetch_campaign(campaign_id: str, conn: sqlite3.Connection) -> Optional[HTSCampaign]:
    row = conn.execute(
        """
        SELECT campaign_id, campaign_name, description, assay_type, target, library_id, total_wells, hit_count, run_date
        FROM hts_campaigns
        WHERE campaign_id = ?
        """,
        (campaign_id,),
    ).fetchone()
    if not row:
        return None
    return HTSCampaign(
        campaign_id=row[0],
        campaign_name=row[1],
        description=row[2],
        assay_type=row[3],
        target=row[4],
        library_id=row[5],
        total_wells=row[6],
        hit_count=row[7],
        run_date=row[8],
    )


def generate_campaign_summary(campaign_id: str, db_path: Optional[str] = None) -> str:
    conn = _get_conn(db_path)
    try:
        campaign = _fetch_campaign(campaign_id, conn)
        if not campaign:
            return f"HTS Campaign {campaign_id} not found."

        total = campaign.total_wells or 0
        hits = campaign.hit_count or 0
        summary = (
            f"HTS Campaign {campaign.campaign_name} screened {total} wells/compounds "
            f"against {campaign.target or 'unknown target'}, found {hits} hits."
        )
        logger.info("[CHEMISTRY][RAG] Generated campaign summary for %s", campaign_id)
        return summary
    finally:
        conn.close()


def embed_compound_for_rag(compound_id: str, db_path: Optional[str] = None) -> Dict[str, str]:
    text = generate_compound_summary(compound_id, db_path=db_path)
    return {
        "text": text,
        "metadata": {"compound_id": compound_id},
    }


def embed_campaign_for_rag(campaign_id: str, db_path: Optional[str] = None) -> Dict[str, str]:
    text = generate_campaign_summary(campaign_id, db_path=db_path)
    return {
        "text": text,
        "metadata": {"campaign_id": campaign_id},
    }


def query_compounds_by_context(query_text: str, top_k: int = 10, db_path: Optional[str] = None):
    """
    Simple keyword search over compound summaries. Placeholder for future vector search.
    """
    conn = _get_conn(db_path)
    try:
        rows = conn.execute("SELECT compound_id FROM compounds").fetchall()
        matches = []
        q_lower = query_text.lower()
        for (cid,) in rows:
            summary = generate_compound_summary(cid, db_path=db_path)
            score = 1.0 if q_lower in summary.lower() else 0.0
            if score > 0:
                matches.append((cid, score))
        matches = sorted(matches, key=lambda x: x[1], reverse=True)[:top_k]
        logger.info("[CHEMISTRY][RAG] query_compounds_by_context matched %d compounds", len(matches))
        return matches
    finally:
        conn.close()

