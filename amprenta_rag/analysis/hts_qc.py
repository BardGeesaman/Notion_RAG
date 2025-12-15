from __future__ import annotations

import math
import statistics
from dataclasses import dataclass, asdict
from typing import Any, Dict, List, Optional
from uuid import UUID

from amprenta_rag.database.models import HTSCampaign, HTSResult
from amprenta_rag.database.session import db_session


@dataclass
class PlateQCSummary:
    campaign_id: UUID
    total_wells: int
    hit_rate: float
    z_prime: Optional[float]
    hits: int
    pos_controls: int
    neg_controls: int

    def asdict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class WellData:
    well_position: Optional[str]
    normalized_value: Optional[float]
    z_score: Optional[float]
    hit_flag: Optional[bool]
    compound_id: Optional[UUID]
    result_id: Optional[str]

    def asdict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class HitCompound:
    result_id: str
    compound_id: UUID
    well_position: Optional[str]
    normalized_value: Optional[float]
    z_score: Optional[float]

    def asdict(self) -> Dict[str, Any]:
        return asdict(self)


def calculate_z_prime(pos_controls: List[float], neg_controls: List[float]) -> Optional[float]:
    """Calculate Z' factor."""
    if len(pos_controls) < 2 or len(neg_controls) < 2:
        return None
    try:
        mu_pos = statistics.fmean(pos_controls)
        mu_neg = statistics.fmean(neg_controls)
        sd_pos = statistics.pstdev(pos_controls)
        sd_neg = statistics.pstdev(neg_controls)
        denom = abs(mu_pos - mu_neg)
        if denom == 0:
            return None
        z_prime = 1 - (3 * (sd_pos + sd_neg) / denom)
        return round(z_prime, 3)
    except Exception:
        return None


def calculate_hit_rate(results: List[HTSResult]) -> float:
    if not results:
        return 0.0
    hits = sum(1 for r in results if getattr(r, "hit_flag", False))
    return round(hits / len(results) * 100.0, 2)


def _split_controls(results: List[HTSResult]) -> tuple[list[float], list[float]]:
    pos_controls: List[float] = []
    neg_controls: List[float] = []
    for r in results:
        cat = (getattr(r, "hit_category", None) or "").lower()
        if cat in ("pos", "positive", "control_pos", "pos_ctrl"):
            if r.normalized_value is not None:
                pos_controls.append(r.normalized_value)
        elif cat in ("neg", "negative", "control_neg", "neg_ctrl"):
            if r.normalized_value is not None:
                neg_controls.append(r.normalized_value)
    return pos_controls, neg_controls


def get_plate_qc_summary(campaign_id: UUID) -> PlateQCSummary:
    with db_session() as db:
        results: List[HTSResult] = (
            db.query(HTSResult).filter(HTSResult.campaign_id == campaign_id).all()
        )
    pos_controls, neg_controls = _split_controls(results)
    z_prime = calculate_z_prime(pos_controls, neg_controls)
    hit_rate = calculate_hit_rate(results)
    return PlateQCSummary(
        campaign_id=campaign_id,
        total_wells=len(results),
        hit_rate=hit_rate,
        z_prime=z_prime,
        hits=sum(1 for r in results if getattr(r, "hit_flag", False)),
        pos_controls=len(pos_controls),
        neg_controls=len(neg_controls),
    )


def get_plate_heatmap_data(campaign_id: UUID) -> List[WellData]:
    with db_session() as db:
        results: List[HTSResult] = (
            db.query(HTSResult).filter(HTSResult.campaign_id == campaign_id).all()
        )
    return [
        WellData(
            well_position=r.well_position,
            normalized_value=r.normalized_value,
            z_score=r.z_score,
            hit_flag=r.hit_flag,
            compound_id=r.compound_id,
            result_id=r.result_id,
        )
        for r in results
    ]


def get_hit_compounds(campaign_id: UUID) -> List[HitCompound]:
    with db_session() as db:
        hits: List[HTSResult] = (
            db.query(HTSResult)
            .filter(HTSResult.campaign_id == campaign_id, HTSResult.hit_flag.is_(True))
            .all()
        )
    return [
        HitCompound(
            result_id=h.result_id,
            compound_id=h.compound_id,
            well_position=h.well_position,
            normalized_value=h.normalized_value,
            z_score=h.z_score,
        )
        for h in hits
    ]

