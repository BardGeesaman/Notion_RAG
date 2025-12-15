from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Dict, List
from uuid import UUID

from amprenta_rag.analysis.quality_metrics import compute_quality_score
from amprenta_rag.database.models import Dataset
from amprenta_rag.database.session import db_session


@dataclass
class DatasetQualityReport:
    dataset_id: UUID
    dataset_name: str
    score: float
    status: str
    issues: List[str]
    metrics: Dict[str, float]

    def asdict(self) -> Dict[str, object]:
        return asdict(self)


def scan_all_datasets() -> List[DatasetQualityReport]:
    """Scan all datasets and compute quality scores."""
    reports: List[DatasetQualityReport] = []
    with db_session() as db:
        datasets: List[Dataset] = db.query(Dataset).all()
        for ds in datasets:
            quality = compute_quality_score(ds)
            reports.append(
                DatasetQualityReport(
                    dataset_id=ds.id,
                    dataset_name=ds.name or str(ds.id),
                    score=quality.get("score", 0.0),
                    status=quality.get("status", "low"),
                    issues=quality.get("issues", []),
                    metrics=quality.get("metrics", {}),
                )
            )
    return reports


def get_low_quality_datasets(threshold: float = 50.0) -> List[DatasetQualityReport]:
    """Return datasets with score below threshold."""
    return [r for r in scan_all_datasets() if r.score < threshold]


def get_quality_summary() -> Dict[str, int]:
    """Return counts by status."""
    reports = scan_all_datasets()
    summary = {"high": 0, "medium": 0, "low": 0, "total": len(reports)}
    for r in reports:
        if r.status == "high":
            summary["high"] += 1
        elif r.status == "medium":
            summary["medium"] += 1
        else:
            summary["low"] += 1
    return summary

