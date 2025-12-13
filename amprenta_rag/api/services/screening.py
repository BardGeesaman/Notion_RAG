from __future__ import annotations

from typing import Dict, List, Optional

from amprenta_rag.database.session import db_session
from amprenta_rag.database.models import HTSCampaign, HTSResult


def list_campaigns() -> List[Dict]:
    """List HTS campaigns ordered by run date descending."""
    with db_session() as db:
        rows = db.query(HTSCampaign).order_by(HTSCampaign.run_date.desc()).all()
        return [_campaign_to_dict(c) for c in rows]


def get_campaign(campaign_id: str) -> Optional[Dict]:
    """Fetch a campaign by campaign_id."""
    with db_session() as db:
        obj = db.query(HTSCampaign).filter(HTSCampaign.campaign_id == campaign_id).first()
        return _campaign_to_dict(obj) if obj else None


def get_campaign_hits(campaign_id: str) -> List[Dict]:
    """Get hit results for a campaign."""
    with db_session() as db:
        rows = (
            db.query(HTSResult)
            .join(HTSCampaign, HTSResult.campaign_id == HTSCampaign.id)
            .filter(HTSCampaign.campaign_id == campaign_id, HTSResult.hit_flag.is_(True))
            .order_by(HTSResult.normalized_value.desc())
            .all()
        )
        return [_result_to_dict(r) for r in rows]


def _campaign_to_dict(c: HTSCampaign) -> Dict:
    if not c:
        return {}
    return {
        "id": str(getattr(c, "id", "")) if getattr(c, "id", None) else None,
        "campaign_id": getattr(c, "campaign_id", None),
        "campaign_name": getattr(c, "campaign_name", None),
        "description": getattr(c, "description", None),
        "assay_type": getattr(c, "assay_type", None),
        "target": getattr(c, "target", None),
        "library_id": getattr(c, "library_id", None),
        "total_wells": getattr(c, "total_wells", None),
        "hit_count": getattr(c, "hit_count", None),
        "run_date": c.run_date.isoformat() if getattr(c, "run_date", None) else None,
    }


def _result_to_dict(r: HTSResult) -> Dict:
    if not r:
        return {}
    return {
        "id": str(getattr(r, "id", "")) if getattr(r, "id", None) else None,
        "result_id": getattr(r, "result_id", None),
        "compound_id": getattr(r, "compound_id", None),
        "well_position": getattr(r, "well_position", None),
        "raw_value": getattr(r, "raw_value", None),
        "normalized_value": getattr(r, "normalized_value", None),
        "z_score": getattr(r, "z_score", None),
        "hit_flag": getattr(r, "hit_flag", None),
        "hit_category": getattr(r, "hit_category", None),
    }

