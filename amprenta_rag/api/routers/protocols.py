from __future__ import annotations

from typing import List
from uuid import UUID

from fastapi import APIRouter, HTTPException

from amprenta_rag.analysis.protocol_diff import (
    audit_deviations,
    diff_protocols,
    get_protocol_history,
)
from amprenta_rag.api import schemas
from amprenta_rag.database.models import Protocol
from amprenta_rag.database.session import db_session

router = APIRouter()


def _get_protocol_or_404(db, protocol_id: UUID) -> Protocol:
    proto = db.query(Protocol).filter(Protocol.id == protocol_id).first()
    if not proto:
        raise HTTPException(status_code=404, detail="Protocol not found")
    return proto


@router.get(
    "/protocols/{protocol_id}/history",
    summary="Get protocol version history",
    response_model=List[schemas.ProtocolHistoryItem],
)
def protocol_history(protocol_id: UUID) -> List[schemas.ProtocolHistoryItem]:
    return [schemas.ProtocolHistoryItem(**item.asdict()) for item in get_protocol_history(protocol_id)]


@router.get(
    "/protocols/{protocol_id}/diff/{other_id}",
    summary="Diff two protocol versions",
    response_model=schemas.ProtocolDiff,
)
def protocol_diff(protocol_id: UUID, other_id: UUID) -> schemas.ProtocolDiff:
    with db_session() as db:
        p1 = _get_protocol_or_404(db, protocol_id)
        p2 = _get_protocol_or_404(db, other_id)
    diff = diff_protocols(p1, p2)
    return schemas.ProtocolDiff(**diff.asdict())


@router.get(
    "/experiments/{experiment_id}/deviations",
    summary="Protocol deviations for an experiment",
    response_model=List[schemas.DeviationReport],
)
def experiment_deviations(experiment_id: UUID) -> List[schemas.DeviationReport]:
    reports = audit_deviations(experiment_id)
    return [schemas.DeviationReport(**r.asdict()) for r in reports]

