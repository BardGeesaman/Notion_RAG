from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Any, Dict, List, Optional, cast
from uuid import UUID

from amprenta_rag.database.models import ExperimentProtocol, Protocol
from amprenta_rag.database.session import db_session


@dataclass
class ProtocolDiff:
    protocol_id: UUID
    other_id: UUID
    added_steps: List[Dict[str, Any]]
    removed_steps: List[Dict[str, Any]]
    changed_steps: List[Dict[str, Any]]
    materials_added: List[Dict[str, Any]]
    materials_removed: List[Dict[str, Any]]
    parameters_changed: List[Dict[str, Any]]

    def asdict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class ProtocolHistoryItem:
    protocol_id: UUID
    version: int
    parent_id: Optional[UUID]
    name: str

    def asdict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class DeviationReport:
    experiment_id: UUID
    protocol_id: UUID
    protocol_name: str
    protocol_version: int
    deviations: List[Any]

    def asdict(self) -> Dict[str, Any]:
        return asdict(self)


def _index_by_order(items: Optional[List[Dict[str, Any]]]) -> Dict[Any, Dict[str, Any]]:
    indexed: Dict[Any, Dict[str, Any]] = {}
    if not items:
        return indexed
    for item in items:
        key = item.get("order") or item.get("id") or item.get("name")
        indexed[key] = item
    return indexed


def diff_protocols(v1: Protocol, v2: Protocol) -> ProtocolDiff:
    """Compute diff between two protocol versions."""
    steps1 = _index_by_order(v1.steps)
    steps2 = _index_by_order(v2.steps)

    added_steps = [steps2[k] for k in steps2.keys() - steps1.keys()]
    removed_steps = [steps1[k] for k in steps1.keys() - steps2.keys()]
    changed_steps: List[Dict[str, Any]] = []
    for key in steps1.keys() & steps2.keys():
        if steps1[key] != steps2[key]:
            changed_steps.append({"from": steps1[key], "to": steps2[key]})

    mats1 = _index_by_order(v1.materials)
    mats2 = _index_by_order(v2.materials)
    materials_added = [mats2[k] for k in mats2.keys() - mats1.keys()]
    materials_removed = [mats1[k] for k in mats1.keys() - mats2.keys()]

    params1 = v1.parameters or {}
    params2 = v2.parameters or {}
    parameters_changed = []
    for key in set(params1.keys()) | set(params2.keys()):
        if params1.get(key) != params2.get(key):
            parameters_changed.append({"name": key, "from": params1.get(key), "to": params2.get(key)})

    return ProtocolDiff(
        protocol_id=cast(UUID, v1.id),
        other_id=cast(UUID, v2.id),
        added_steps=added_steps,
        removed_steps=removed_steps,
        changed_steps=changed_steps,
        materials_added=materials_added,
        materials_removed=materials_removed,
        parameters_changed=parameters_changed,
    )


def get_protocol_history(protocol_id: UUID) -> List[ProtocolHistoryItem]:
    """Return protocol version chain following parent_id."""
    history: List[ProtocolHistoryItem] = []
    with db_session() as db:
        current = db.query(Protocol).filter(Protocol.id == protocol_id).first()
        while current:
            history.append(
                ProtocolHistoryItem(
                    protocol_id=cast(UUID, current.id),
                    version=current.version or 1,
                    parent_id=cast(Optional[UUID], current.parent_id),
                    name=current.name or str(current.id),
                )
            )
            if current.parent_id:
                current = db.query(Protocol).filter(Protocol.id == current.parent_id).first()
            else:
                break
    # newest to oldest; reverse to chronological oldest->newest
    return list(reversed(history))


def audit_deviations(experiment_id: UUID) -> List[DeviationReport]:
    """Collect deviations for an experiment across linked protocols."""
    reports: List[DeviationReport] = []
    with db_session() as db:
        links: List[ExperimentProtocol] = (
            db.query(ExperimentProtocol)
            .filter(ExperimentProtocol.experiment_id == experiment_id)
            .all()
        )
        for link in links:
            deviations = link.deviations or []
            proto = link.protocol
            reports.append(
                DeviationReport(
                    experiment_id=experiment_id,
                    protocol_id=cast(UUID, proto.id) if proto else cast(UUID, link.protocol_id),
                    protocol_name=str(proto.name) if proto and proto.name is not None else str(link.protocol_id),
                    protocol_version=int(proto.version) if proto and proto.version is not None else 1,
                    deviations=deviations if isinstance(deviations, list) else [deviations],
                )
            )
    return reports

